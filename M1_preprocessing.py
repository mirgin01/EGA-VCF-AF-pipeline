import hail as hl
import logging
import os
from utils import *
from hail.ggplot import ggplot, aes, labs


load_logging()
config = load_config()

def convert_and_merge_vcfs(vcf_dir, mt_from_vcf, reference_genome):

    # Get all VCF file paths in the directory
    vcfs = sorted([
        os.path.join(vcf_dir, f)
        for f in os.listdir(vcf_dir)
        if f.endswith('.vcf') or f.endswith('.vcf.bgz') or f.endswith('.vcf.gz')
    ])

    if not vcfs:
        csv_writer(["No VCF files found in the directory"])
        raise ValueError("No VCF files found in the directory.")

    logging.info(f"Found {len(vcfs)} VCF files to merge.")

    if len(vcfs) == 1:
        logging.info(f"Converting: {vcfs[0]}")
        mt = hl.import_vcf(vcfs[0], reference_genome=reference_genome, min_partitions=4)

    else:
        # Load the first VCF as the initial MatrixTable
        logging.info(f"Converting: {vcfs}")
        mt = hl.import_vcf(f"{config['vcf_dir']}*.vcf*", reference_genome=reference_genome, min_partitions=4)
        
    mt.write(mt_from_vcf, overwrite=True)
    
    logging.info(f"MatrixTable with {vcfs} written to {mt_from_vcf}   ")
    summary = [f"VCFs converted to {config['mt_from_vcf']}:", vcfs]
    csv_writer(summary)
    

    return mt


def impute_sex(mt):
    
    imputed_sex = hl.impute_sex(mt.GT, aaf_threshold=0.05, female_threshold=0.5, male_threshold=0.75)  # Imputed sex with suggested thresholds
    sex_expr = hl.if_else(hl.is_defined(imputed_sex.is_female), hl.if_else(imputed_sex.is_female, # rename imputed sex
                                                                        "female",
                                                                        "male"),
                                                                        "undefined")
    sex_ht = imputed_sex.annotate(imputed_sex=sex_expr)

    # annotate input (all chroms) mt with imputed sex
    sex_colnames = ["f_stat", "is_female", "imputed_sex"]
    mt = mt.annotate_cols(**sex_ht.select(*sex_colnames)[mt.col_key])

    # Save the sex imputation results
    basename = os.path.basename(config['mt_from_vcf'].rstrip('/'))
    sex_ht.select("imputed_sex", "f_stat", "is_female").export(f"{basename}_imputed_sex_results.tsv")

    return mt 

def split_multiallelic(mt, original_size):
    summary = []
    summary.append(["Split multiallelic variants"])
    summary.append(["", "Before", "After"])

    mt = hl.split_multi_hts(mt)

    if config['verbosity']: 
        after_splitting = mt.count()
        summary.append(["Splitting multiallelic variants", original_size[0], after_splitting[0]])
        logging.info(f"Number of multiallelic variants: {original_size[0] - after_splitting[0]}   ")
        csv_writer(summary)
    
    return mt

def genotype_filtering(mt):

    # Store filtering summary
    summary = []
    summary.append([])
    summary.append(["Genotype Filtering Steps"])
    
    # Genotype quality filtering
    if config['genotype_filters']['GQ_threshold'] is not None:
        if config['verbosity']:
            if "GT" in mt.entry:
                # General information about Genotype
                general_stats = mt.aggregate_entries(
                    hl.struct(
                        total=hl.agg.count(),
                        defined_GT=hl.agg.count_where(hl.is_defined(mt.GT)), # stats about samples without GT info
                    )
                )
                logging.info(f"Total number of entries: {general_stats.total}   ")
                logging.info(f"Total number of entries with GT information: {general_stats.defined_GT}  ")

        
        if "GQ" in mt.entry:
            if config['verbosity']:
                # Count total and passing genotypes in one go
                stats = mt.aggregate_entries(hl.struct(
                    passing=hl.agg.count_where(mt.GQ >= config['genotype_filters']['GQ_threshold']),
                    defined_GQ=hl.agg.count_where(hl.is_defined(mt.GQ)),
                    defined_both_GT_GQ=hl.agg.count_where(hl.is_defined(mt.GT) & hl.is_defined(mt.GQ)),
                    basic_stats=hl.agg.stats(mt.GQ)
                    )
                )
                # Stats for summary 
                clean_stats = {
                    'entries_with_GQ': stats.defined_GQ,
                    'entries_with_both_GT_GQ': stats.defined_both_GT_GQ,
                    'mean': stats.basic_stats.mean,
                    'stdev': stats.basic_stats.stdev,
                    'min': stats.basic_stats.min,
                    'max': stats.basic_stats.max,
                    'n': stats.basic_stats.n
                }
                
                # number of remaining genotypes
                genotypes_filtered = general_stats.defined_GT - stats.passing
                
                # create log
                logging.info(f"Total number of entries with GQ: {stats.defined_GQ}  ")
                logging.info(f"Total number of entries with GT and GQ: {stats.defined_both_GT_GQ}   ")
                logging.info(f"Stats for GQ field:{clean_stats}   ")
                logging.info(f"GQ filtering done - Genotypes removed: {genotypes_filtered}  ")
                
                # create summary
                summary.append(["Filter", "Before", "After", "Removed", "Stats"])
                summary.append(["GQ", general_stats.defined_GT, stats.passing, genotypes_filtered, clean_stats])

            mt = mt.annotate_entries(
                GT = hl.if_else(mt.GQ >= config['genotype_filters']['GQ_threshold'], mt.GT, hl.null(mt.GT.dtype))
                ) # set non passing GT to NULL 
        else:
            logging.info("GQ information not available")
            summary.append(["Filter", "Before", "After", "Removed", "Status"])
            summary.append(["GQ", "-", "-", "-", "Not available"])

    # AB filtering
    if config['genotype_filters']['AB_threshold'] is not None:
        if "GT" in mt.entry and "AD" in mt.entry:
            
            # Apply AB annotation and filtering
            mt = mt.annotate_entries(
                AB=hl.or_missing(
                    mt.GT.is_het(),
                    mt.AD[1] / hl.sum(mt.AD) # Calculate AB 
                )
            )

            if config['verbosity']:
                # Count genotypes before AB filtering
                pre_ab_stats = mt.aggregate_entries(
                    hl.struct(
                        defined_GT=hl.agg.count_where(hl.is_defined(mt.GT)),
                        het_genotypes=hl.agg.count_where(mt.GT.is_het()),
                        defined_GT_AD=hl.agg.count_where(hl.is_defined(mt.GT) & hl.is_defined(mt.AD)),
                        basic_stats=hl.agg.stats(mt.AB)
                    )
                )
                
                # Stats for summary 
                clean_stats = {
                    'entries_with_both_GT_AD': pre_ab_stats.defined_GT_AD, 
                    'het_genotypes': pre_ab_stats.het_genotypes,
                    'mean': pre_ab_stats.basic_stats.mean,
                    'stdev': pre_ab_stats.basic_stats.stdev,
                    'min': pre_ab_stats.basic_stats.min,
                    'max': pre_ab_stats.basic_stats.max,
                    'n': pre_ab_stats.basic_stats.n
                }
            
            mt = mt.annotate_entries( # set non passing GT to NULL 
                GT = hl.case()
                .when(mt.GT.is_het() & (mt.AB < config['genotype_filters']['AB_threshold']), hl.null(mt.GT.dtype)) # filter ONLY het genotypes by AB
                .default(mt.GT)
            )
            
            if config['verbosity']:
                # Count genotypes after AB filtering
                post_ab_stats = mt.aggregate_entries(
                    hl.struct(
                        defined_GT=hl.agg.count_where(hl.is_defined(mt.GT)),
                        remaining_genotypes=hl.agg.count_where(hl.is_defined(mt.GT))
                    )
                )
                
                genotypes_removed_ab = pre_ab_stats.defined_GT - post_ab_stats.defined_GT
                
                logging.info(f"Total heterozygous genotypes: {pre_ab_stats.het_genotypes}   ")
                logging.info(f"Total genotypes with GT and AD: {pre_ab_stats.defined_GT_AD} ")
                logging.info(f"General AB statistics: {pre_ab_stats.basic_stats   }")
                logging.info(f"Allele Balance filtering done - Genotypes removed: {genotypes_removed_ab}    ")
                
                summary.append(["AlleleBalance", pre_ab_stats.defined_GT, post_ab_stats.defined_GT, genotypes_removed_ab, clean_stats])
        else:
            logging.info("Allele Balance information not available")
            summary.append(["Filter", "Before", "After", "Removed", "Status"])
            summary.append(["AlleleBalance (Genotypes)", "-", "-", "-", "Not available"])
    
      
    csv_writer(summary)

    
    return mt


def variant_filtering(mt):
    
    mt = hl.variant_qc(mt)
   
    # Store filtering summary
    summary = []
    summary.append([])
    summary.append(["Variant Filtering Steps"])
    summary.append(["Filter", "Before", "After", "Removed", "Stats"])

      
    # Variant quality control
    
    if config['variant_filters']['QD_threshold'] is not None:
        if "QD" in mt.info:
            if config['verbosity']:
                pre_count = mt.count_rows()
                stats = create_stats(mt, "info.QD", "rows")
                
            mt = mt.filter_rows(mt.info.QD >= config['variant_filters']['QD_threshold'])

            if config['verbosity']:
                post_count = mt.count_rows()
                removed = pre_count - post_count
                logging.info(f"QD filtering done - Variants removed: {removed}  ")
                summary.append(["QD", pre_count, post_count, removed, stats])
            else:
                logging.info("QD filtering done")
        else:
            logging.info("QD information not available")
            summary.append(["QD", "-", "-", "-", "Not available"])
    
    if config['variant_filters']['DP_threshold'] is not None:
        if "DP" in mt.info:
            if config['verbosity']:
                pre_count = mt.count_rows()
                stats = create_stats(mt, "info.DP", "rows")

            mt = mt.filter_rows(mt.info.DP >= config['variant_filters']['DP_threshold'])
            
            if config['verbosity']:
                post_count = mt.count_rows()
                removed = pre_count - post_count
                logging.info(f"DP filtering done - Variants removed: {removed}  ")
                summary.append(["DP", pre_count, post_count, removed, stats])
            else:
                logging.info("DP filtering done")
        else:
            logging.info("DP information not available")
            summary.append(["DP", "-", "-", "-", "Not available"])

    if config['variant_filters']['QUAL_threshold'] is not None:
        try:
            if config['verbosity']:
                pre_count = mt.count_rows()
                stats = create_stats(mt, "qual", "rows")
            
            mt = mt.filter_rows(mt.qual >= config['variant_filters']['QUAL_threshold'])
            
            if config['verbosity']:
                post_count = mt.count_rows()
                removed = pre_count - post_count
                logging.info(f"QUAL filtering done - Variants removed: {removed}    ")
                summary.append(["QUAL", pre_count, post_count, removed, stats])
            else:
                logging.info("QUAL filtering done")
        except TypeError:
            logging.info("QUAL information not available")
            summary.append(["QUAL", "-", "-", "-", "Not available"])

    if config['variant_filters']['MQ_threshold'] is not None:
        if "MQ" in mt.info:
            if config['verbosity']:
                pre_count = mt.count_rows()
                stats = create_stats(mt, "info.MQ", "rows")
            
            mt = mt.filter_rows(mt.info.MQ >= config['variant_filters']['MQ_threshold'])
            
            if config['verbosity']:
                post_count = mt.count_rows()
                removed = pre_count - post_count
                logging.info(f"MQ filtering done - Variants removed: {removed}  ")
                summary.append(["MQ", pre_count, post_count, removed, stats])
            else:
                logging.info("MQ filtering done")
        else:
            logging.info("MQ information not available")
            summary.append(["MQ", "-", "-", "-", "Not available"])

    if config['variant_filters']['FS_threshold'] is not None:
        if "FS" in mt.info:
            if config['verbosity']:
                pre_count = mt.count_rows()
                stats = create_stats(mt, "info.FS", "rows")

            mt = mt.filter_rows(mt.info.FS <= config['variant_filters']['FS_threshold'])
            
            if config['verbosity']:
                post_count = mt.count_rows()
                removed = pre_count - post_count
                logging.info(f"FS filtering done - Variants removed: {removed}  ")
                summary.append(["FS", pre_count, post_count, removed, stats])
            else:
                logging.info("FS filtering done")
        else:
            logging.info("FS information not available")
            summary.append(["FS", "-", "-", "-", "Not available"])

    if config['variant_filters']['READPOSRANKSUM_threshold'] is not None:
        if "ReadPosRankSum" in mt.info:
            if config['verbosity']:
                pre_count = mt.count_rows()
                stats = create_stats(mt, "info.ReadPosRankSum", "rows")
            
            mt = mt.filter_rows(mt.info.ReadPosRankSum >= config['variant_filters']['READPOSRANKSUM_threshold'])
            
            if config['verbosity']:
                post_count = mt.count_rows()
                removed = pre_count - post_count
                logging.info(f"ReadPosRankSum filtering done - Variants removed: {removed}  ")
                summary.append(["ReadPosRankSum", pre_count, post_count, removed, stats])
            else:
                logging.info("ReadPosRankSum filtering done")
        else:
            logging.info("ReadPosRankSum information not available")
            summary.append(["ReadPosRankSum", "-", "-", "-", "Not available"])

    csv_writer(summary)


    return mt


def sample_filtering(mt, sequencingType): ## TODO print stats about qual fields as in sample and genotype

    # ensure the sequencing type entered in conf.py is valid
    assert sequencingType in ["WES", "WGS"], "sequencingType must be 'WES' or 'WGS'"

    # Compute sample QC
    mt = hl.sample_qc(mt)

    # Compute CHARR
    """if "AF" in mt.info:
        gnomad_sites = hl.read_table(config['gnomad_sites'])
        charr_result = hl.compute_charr(
            mt,
            ref_AF=(1 - gnomad_sites[mt.row_key].freq[0].AF)
        )
        #annotate charr results
        mt = mt.annotate_cols(charr=charr_result[mt.col_key].charr)"""


    config['verbosity'] = True

    if config["plots"]:
        logging.info("Creating Sample Filters plots")
        create_sample_plot(mt)

    if config['verbosity']:
        summary = []
        summary.append([])
        summary.append(["Sample Filtering Steps"])
        summary.append(["Filter", "Before", "After", "Removed", "Stats"])

    # Minimum coverage filtering
    if config['sample_filters']['DP_STATS.MEAN_WES_threshold'] is not None or config['sample_filters']['DP_STATS.MEAN_WGS_threshold'] is not None:
        if "dp_stats" in mt.sample_qc:
            if config['verbosity']:
                pre_count = mt.count()
                stats = create_stats(mt, "sample_qc.dp_stats.mean", "cols")
            
            min_coverage = config['sample_filters']['DP_STATS.MEAN_WES_threshold'] if sequencingType == "WES" else config['sample_filters']['DP_STATS.MEAN_WGS_threshold']
            mt = mt.filter_cols(mt.sample_qc.dp_stats.mean >= min_coverage)
            
            if config['verbosity']:
                post_count = mt.count()
                summary.append(["Minimum Coverage", pre_count[1], post_count[1], pre_count[1] - post_count[1], stats])
                logging.info(f"{sequencingType} coverage filtering done - Samples removed: {pre_count[1] - post_count[1]}   ")
            else:
                logging.info("DP filtering done")
        else:
            summary.append(["Minimum Coverage", "-", "-", "-", "Not available"])
            logging.info("DP information not available")

    # ti/tv ratio filtering
    if config['sample_filters']['R_TI_TV_WES_threshold'] is not None or config['sample_filters']['R_TI_TV_WGS_threshold'] is not None:
        if "r_ti_tv" in mt.sample_qc:
            if config['verbosity']:
                pre_count = mt.count()
                stats = create_stats(mt, "sample_qc.r_ti_tv", "cols")
            lower, upper = config['sample_filters']['R_TI_TV_WES_threshold'] if sequencingType == "WES" else config['sample_filters']['R_TI_TV_WGS_threshold']
            mt = mt.filter_cols((mt.sample_qc.r_ti_tv >= lower) & (mt.sample_qc.r_ti_tv <= upper))
            if config['verbosity']:
                post_count = mt.count()
                summary.append(["Ti/Tv Ratio", pre_count[1], post_count[1], pre_count[1] - post_count[1], stats])
                logging.info(f"{sequencingType} ti/tv filtering done — Samples removed: {pre_count[1] - post_count[1]}")
        else:
            summary.append(["Ti/Tv Ratio", "-", "-", "-", "Not available"])
            logging.info("Ti/Tv information not available")

    # Call rate filtering
    if config['sample_filters']['CALL_RATE_threshold'] is not None: 
        if "call_rate" in mt.sample_qc:
            if config['verbosity']:
                pre_count = mt.count()
                stats = create_stats(mt, "sample_qc.call_rate", "cols")

            mt = mt.filter_cols(mt.sample_qc.call_rate >= config['sample_filters']['CALL_RATE_threshold'])
            if config['verbosity']:    
                post_count = mt.count()
                summary.append(["Call Rate", pre_count[1], post_count[1], pre_count[1] - post_count[1], stats] )
                logging.info(f"Call rate filtering done - Samples removed: {pre_count[1] - post_count[1]}"  )
            else:
                logging.info("Call rate filtering done")
        else:
            summary.append(["Call Rate", "-", "-", "-", "Not available"])
            logging.info("Call rate information not available")

    # Singletons filtering
    if config['sample_filters']['N_SINGLETON_WES_threshold'] is not None or config['sample_filters']['R_TI_TV_WGS_threshold'] is not None:
        if "n_singleton" in mt.sample_qc: 
            if config['verbosity']:
                pre_count = mt.count()
                stats = create_stats(mt, "sample_qc.n_singleton", "cols")

            hard_cutoff = config['sample_filters']['N_SINGLETON_WES_threshold'] if sequencingType == "WES" else config['sample_filters']['N_SINGLETON_WGS_threshold']
            mt = mt.filter_cols(mt.sample_qc.n_singleton <= hard_cutoff)
            
            if config['verbosity']:
                post_count = mt.count()
                summary.append(["Singletons", pre_count[1], post_count[1], pre_count[1] - post_count[1], stats])
                logging.info(f"Singleton filtering done - Samples removed: {pre_count[1] - post_count[1]}   ")
            else:
                logging.info("Singleton filtering done")
        else:
            summary.append(["Singletons", "-", "-", "-", "Not available"])
            logging.info("Number of singletons not available")


    # CHARR filtering 
       
    """
    if config['sample_filters']['CHARR_threshold'] is not None:
        try:
            if config['verbosity']:
                pre_count = mt.count_rows()
                stats = create_stats(mt, "charr", "cols")
                
            mt = mt.filter_cols(mt.charr <= config['sample_filters']['CHARR_threshold']) # Filter samples with charr over the threshold
            
            if config['verbosity']:
                post_count = mt.count()
                logging.info(f"CHARR filtering done - Samples removed: {pre_count[1] - post_count[1]}   ")
                summary.append(["CHARR", pre_count[1], post_count[1], pre_count[1] - post_count[1], stats])  
            else: 
                logging.info("CHARR filtering done")
        except TypeError:
            logging.info("CHARR couldn't be calculated")
            summary.append(["CHARR", "-", "-", "-", "Not available"])"""


    # Het/Hom ratio filtering
    if config['sample_filters']['R_HET_HOM_VAR_WES_threshold'] is not None or config['sample_filters']['R_HET_HOM_VAR_WGS_threshold'] is not None:
        if "GT" in mt.entry:
            if config['verbosity']:
                pre_count = mt.count()
                stats = create_stats(mt, "sample_qc.r_het_hom_var", "cols")
            
            hard_cutoff = config['sample_filters']['R_HET_HOM_VAR_WES_threshold'] if sequencingType == "WES" else config['sample_filters']['R_HET_HOM_VAR_WGS_threshold']
            mt = mt.filter_cols(mt.sample_qc.r_het_hom_var <= hard_cutoff)
            
            if config['verbosity']:
                post_count = mt.count()
                summary.append(["Het/Hom Ratio", pre_count[1], post_count[1], pre_count[1] - post_count[1], stats])

                logging.info(f"{sequencingType} het/hom filtering done - Samples removed: {pre_count[1] - post_count[1]}    ")
            else:
                logging.info("Het/Hom ratio filtering done")
        else:
            summary.append(["Het/Hom Ratio", "-", "-", "-", "Not available"])
            logging.info("Het/Hom ratio information not available")

    if config['verbosity']:
        # Write CSV
        csv_writer(summary)
    
    return mt
