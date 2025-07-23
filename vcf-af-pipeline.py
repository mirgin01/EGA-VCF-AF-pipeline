from M1_preprocessing import convert_and_merge_vcfs, split_multiallelic, genotype_filtering, variant_filtering, sample_filtering
from M2_unrelated_samples import delete_related_samples
from M3_ancestry import subset_matrix, call_grafanc, annotate_ancestry
from M4_af_annotation import stats_by_sex, stats_by_ancestry, af_by_sex_ancestry, annotate_new_vcf, export_new_vcf
from utils import *
import hail as hl 

hl.init(
    spark_conf={
        'spark.driver.memory': config['spark_driver_memory'],  # Allocate sufficient memory for the driver
        'spark.executor.memory': config['spark_executor_memory'],  # Allocate memory for executors
        'spark.executor.instances': config['spark_executor_instances'],  # Use multiple executors
        'spark.executor.cores': config['spark_executor_cores'],  # Assign cores per executor
        'spark.rpc.askTimeout': config['spark_rpc_askTimeout'],  # Increase timeout for slow operations
        'spark.sql.shuffle.partitions': config['spark_sql_shuffle_partitions'],  # Reduce shuffle partitions for large data
        'spark.memory.fraction': config ['spark_memory_fraction'],  # Use most of the JVM heap for Spark execution
        'spark.local.dir': config["spark_local_dir"],  # Specify a temp directory for disk spill
        'spark.network.timeout': config["spark_network_timeout"],  # Avoid network timeouts
    },
    tmp_dir = config["tmp_dir"],
    local_tmpdir = config["local_tmpdir"]
)



def main():
    
    config = load_config()
    load_logging() # create proper logs 

    logging.info(f"+++ Running EGA standard VCF workflow v1 +++")

    csv_creator() 
    summary = []

    ## PREPROCESSING
    
    if config['preprocessing']:

        if config['convert_vcfs']:
            mt = convert_and_merge_vcfs(
                config['vcf_dir'],
                config['mt_from_vcf'],
                config['ref_gen']
            )

        mt = hl.read_matrix_table(config['mt_from_vcf'])

        original_size = mt.count() 
        logging.info(f"Original dataset size. Variants: {original_size[0]}, Samples: {original_size[1]}   ")

        if config['split_multiallelic']:
            mt = split_multiallelic(mt, original_size)
        if config['genotype_filtering']:
            mt = genotype_filtering(mt)
        if config['variant_filtering']:
            mt = variant_filtering(mt)

        # Impute sex for the samples 
        mt = impute_sex(mt)

        if config['sample_filtering']:
            mt = sample_filtering(mt, config['seq_type'])
                    

        final_size = mt.count()
        logging.info(f"After QC dataset size. Variants: {final_size[0]}, Samples: {final_size[1]}   ")

        mt.write(config['mt_afterQC'], overwrite=True)  # write matrix with QC 
        logging.info(f"MT with QC written in: {config['mt_afterQC']}")
        

    else:
        mt = hl.read_matrix_table(config['mt_from_vcf'])

        original_size = mt.count()
        
        logging.info(f"Original dataset size. Variants: {original_size[0]}, Samples: {original_size[1]}   ")
        logging.warning(f"No quality control applied to variants, samples and genotypes. Working with {config['mt_from_vcf']}")

        summary = []
        summary.append([])
        summary.append(["", "Variants", "Samples"])
        summary.append(["Original dataset size", original_size[0], original_size[1]])
        csv_writer(summary)

    ## DELETE RELATED SAMPLES
    
    if config["delete_related"]:
        logging.info(f"Deleting minimum number of samples for an unrelated dataset")
        mt = delete_related_samples(mt)
    
    ## ANCESTRY

    if config['ancestry']:
        ancestry_vcf = subset_matrix(mt)
        if config['preprocessing']:
            logging.info(f"Ancestry inference done using {config['mt_afterQC']}")
            ancestry_results = call_grafanc(ancestry_vcf, config['mt_afterQC']) # Annotate Anc to the correct matrix depending on conf.py
        else:
            logging.info(f"Ancestry inference done using {config['mt_from_vcf']}")
            ancestry_results = call_grafanc(ancestry_vcf, config['mt_from_vcf'])
        
        if ancestry_results is not None: # if ancestry was inferred annotate the results 
            mt = annotate_ancestry(ancestry_results, mt)
        

    ## AF RECALC

    if config['af_annotation']:
        if config['preprocessing'] == False:
            # Impute sex for the samples 
            mt = impute_sex(mt)
        mt, results_sex_agg = stats_by_sex(mt)
        if config['ancestry']:
            try:
                mt, results_ancestry_agg = stats_by_ancestry(mt)
            except: 
                logging.warning("AFs wont be grouped by ancestry")
                results_ancestry_agg = ""
        else:
            results_ancestry_agg = ""

        mt, AF_total, AC_total, AN_total, AC_hom_total, AF_sex, AF_ancestry = af_by_sex_ancestry(mt,
                                                                                             results_sex_agg,
                                                                                             results_ancestry_agg)

        mt = annotate_new_vcf(mt, AF_total, AC_total, AN_total, AC_hom_total, AF_sex, AF_ancestry)
        export_new_vcf(mt)




if __name__ == "__main__":
    main()