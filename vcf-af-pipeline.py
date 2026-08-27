from M1_preprocessing import convert_and_merge_vcfs, split_multiallelic, genotype_filtering, variant_filtering, sample_filtering_mad_thresholds, sample_filtering_hard_thresholds
from M2_unrelated_samples import delete_related_samples
from M3_ancestry import subset_matrix, call_grafanc, annotate_ancestry
from M4_af_annotation import stats_by_sex, stats_by_ancestry, af_by_sex_ancestry, annotate_new_vcf, export_new_vcf
from utils import *
import hail as hl 

if config['cluster']:
    hl.init(
        spark_conf={
            'spark.driver.memory': config['spark_driver_memory'],  # Allocate sufficient memory for the driver
            'spark.executor.memory': config['spark_executor_memory'],  # Allocate memory for executors
            'spark.executor.instances': config['spark_executor_instances'],  # Use multiple executors
            'spark.executor.cores': config['spark_executor_cores'],  # Assign cores per executor
            'spark.rpc.askTimeout': config['spark_rpc_askTimeout'],  # Increase timeout for slow operations
            'spark.sql.shuffle.partitions': config['spark_sql_shuffle_partitions'],  # Reduce shuffle partitions for large data
            'spark.local.dir': config["spark_local_dir"],  # Specify a temp directory for disk spill
            'spark.network.timeout': config["spark_network_timeout"],  # Avoid network timeouts
        },
        tmp_dir = config["tmp_dir"],
        local_tmpdir = config["local_tmpdir"]
    )
else:
    hl.init()



def main():
    
    config = load_config()

    logging.info(f"=== Pipeline settings ===")

    print_config(config)

    logging.info(f"=========================")
    logging.info(f"+++ RUNNING EGA STANDARD VCF WORKFLOW v1.1 +++")
    csv_creator() 
    summary = []
    basename = os.path.basename(config['mt_from_vcf'].rstrip('/'))

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

        mt = rename_chr(mt, config['ref_gen'])

        if config['split_multiallelic']:
            mt = split_multiallelic(mt, original_size)
        if config['genotype_filtering']:
            mt = genotype_filtering(mt)
        if config['variant_filtering']:
            mt = variant_filtering(mt)
        if config['sample_filtering']:
            mt = sample_filtering_hard_thresholds(mt, config['seq_type'], basename)
            mt = sample_filtering_mad_thresholds(mt, config['seq_type'], basename)
                    
        final_size = mt.count()
        logging.info(f"After QC dataset size. Variants: {final_size[0]}, Samples: {final_size[1]}   ")

        if final_size[0] == 0:
            logging.error("All the variants were deleted by the quality control - stopping the pipeline")
            raise ValueError(f"Empty dataset - ALL the variants were deleted by the quality control - Aborting execution")
        if final_size[1] == 0:
            logging.error("All the samples were deleted by the quality control - stopping the pipeline")
            raise ValueError(f"Empty dataset - ALL the samples were deleted by the quality control - Aborting execution")
        else:
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

        rename_chr(mt, config['ref_gen'])
 
    
    ## DELETE RELATED SAMPLES
    
    if config["delete_related"]:
        size = mt.count()
        if size[1] == 1: 
            logging.info (f"There's only one sample in the dataset - skipping trimming of related individuals")
        else:
            logging.info(f"Deleting minimum number of samples for an unrelated dataset")
            mt = delete_related_samples(mt)
    
    ## ANCESTRY

    if config['infer_ancestry']:
        ancestry_vcf, ancestry_snps_count = subset_matrix(mt)
        if ancestry_snps_count == 0: 
            logging.error("There are not enough ancestry SNPs to run GrafAnc - aborting ancestry annotation")
        else: 
            if config['preprocessing']:
                logging.info(f"Ancestry inference done using {config['mt_afterQC']}")
                ancestry_results = call_grafanc(ancestry_vcf, config['mt_afterQC']) # Annotate Anc to the correct matrix depending on conf.py
            else:
                logging.info(f"Ancestry inference done using {config['mt_from_vcf']}")
                ancestry_results = call_grafanc(ancestry_vcf, config['mt_from_vcf'])
            
            if ancestry_results is not None: # if ancestry was inferred annotate the results 
                mt = annotate_ancestry(ancestry_results, mt)

    if config['submit_ancestry']:
        logging.info(f"Ancestry/population read from: {config['ancestry_information']}")
        mt = annotate_ancestry(config['ancestry_information'], mt) 

        

    ## AF RECALC

    if config['af_annotation']:
        
        # Create sex groups
        if config['infer_sex']:
            mt = impute_sex(mt)
            mt, results_sex_agg = stats_by_sex(mt)
        else:            
            logging.warning("Imputed_sex not found. Skipping sex-stratified stats. This avoids the calculation of hemizygous calls. ")
            results_sex_agg = {} 
                
        # Create ancestry groups
        if config['infer_ancestry'] or config['submit_ancestry']:
            try:
                mt, results_ancestry_agg = stats_by_ancestry(mt)
            except: 
                logging.warning("AFs wont be grouped by ancestry")
                results_ancestry_agg = ""
        else:
            results_ancestry_agg = ""

        mt, AF_total, AC_total, AN_total, nhom_total, nhet_total, nhemi_total, AF_sex, AF_ancestry, AF_ancestry_sex = af_by_sex_ancestry(mt,
                                                                                             results_sex_agg,
                                                                                             results_ancestry_agg)

        mt = annotate_new_vcf(mt, AF_total, AC_total, AN_total, nhom_total, nhet_total, nhemi_total, AF_sex, AF_ancestry, AF_ancestry_sex)
        
        export_new_vcf(mt)




if __name__ == "__main__":
    load_logging() # create proper logs 
    main()