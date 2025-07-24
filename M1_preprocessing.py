import hail as hl
import logging
import os
from utils import *
from hail.ggplot import ggplot, aes, labs


config = load_config()

def convert_and_merge_vcfs(vcf_dir, mt_from_vcf, reference_genome):
    """
    Convert all the VCFs in a folder to a SINGLE hail matrix
    :params: directory where the VCFs are saved, path where the hail matrix will be saved, reference genome from the VCFs
    :return: Hail matrix with the genomic data converted
    """
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


def split_multiallelic(mt, original_size):
    """
    Splits multiallelic variants
    :params: mt without qc, original size to avoid recalculating
    :return: mt without multiallelic variants
    """
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
    """
    Applies genotype quality control based on the thresholds stated in config.yaml
    :params: mt without genotype qc
    :return: mt with genotypes QCed
    """
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
    """
    Applies variant quality control based on the thresholds stated in config.yaml
    :params: mt without variant qc
    :return: mt with variant QCed
    """   
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
                stats = create_stats(mt, "info.QD", "rows")            
                summary = verbosity_counts_variants(mt, "QD", mt.info.QD >= config['variant_filters']['QD_threshold'], summary, stats)                      

            # actually apply the filter to the mt
            mt = mt.filter_rows(mt.info.QD >= config['variant_filters']['QD_threshold'])
            logging.info("QD filtering done")
        
        else:
            logging.info("QD information not available")
            summary.append(["QD", "-", "-", "-", "Not available"])
    
    if config['variant_filters']['DP_threshold'] is not None:
        
        if "DP" in mt.info:
            
            if config['verbosity']:
                stats = create_stats(mt, "info.DP", "rows")
                summary = verbosity_counts_variants(mt, "DP", mt.info.DP >= config['variant_filters']['DP_threshold'], summary, stats)

            mt = mt.filter_rows(mt.info.DP >= config['variant_filters']['DP_threshold'])       
            logging.info("DP filtering done")   
            
        else:
            logging.info("DP information not available")
            summary.append(["DP", "-", "-", "-", "Not available"])

    if config['variant_filters']['QUAL_threshold'] is not None:
        try:
            if config['verbosity']:

                stats = create_stats(mt, "qual", "rows")
                summary = verbosity_counts_variants(mt, "QUAL", mt.qual >= config['variant_filters']['QUAL_threshold'], summary, stats)
            
            mt = mt.filter_rows(mt.qual >= config['variant_filters']['QUAL_threshold'])
            logging.info("QUAL filtering done")
        
        except TypeError:
            logging.info("QUAL information not available")
            summary.append(["QUAL", "-", "-", "-", "Not available"])

    if config['variant_filters']['MQ_threshold'] is not None:
        if "MQ" in mt.info:
            if config['verbosity']:
                stats = create_stats(mt, "info.MQ", "rows")
                summary = verbosity_counts_variants(mt, "MQ", mt.info.MQ >= config['variant_filters']['MQ_threshold'], summary, stats)
               
            mt = mt.filter_rows(mt.info.MQ >= config['variant_filters']['MQ_threshold'])
            logging.info("MQ filtering done")
        else:
            logging.info("MQ information not available")
            summary.append(["MQ", "-", "-", "-", "Not available"])

    if config['variant_filters']['FS_threshold'] is not None:
        if "FS" in mt.info:
            if config['verbosity']:
                stats = create_stats(mt, "info.FS", "rows")
                summary = verbosity_counts_variants(mt, "FS", mt.info.FS <= config['variant_filters']['FS_threshold'], summary, stats)

            mt = mt.filter_rows(mt.info.FS <= config['variant_filters']['FS_threshold'])
            logging.info("FS filtering done")
        else:
            logging.info("FS information not available")
            summary.append(["FS", "-", "-", "-", "Not available"])

    if config['variant_filters']['READPOSRANKSUM_threshold'] is not None:
        if "ReadPosRankSum" in mt.info:
            if config['verbosity']:
                stats = create_stats(mt, "info.ReadPosRankSum", "rows")
                summary = verbosity_counts_variants(mt, "ReadPosRankSum", mt.info.ReadPosRankSum >= config['variant_filters']['READPOSRANKSUM_threshold'], summary, stats)
            
            mt = mt.filter_rows(mt.info.ReadPosRankSum >= config['variant_filters']['READPOSRANKSUM_threshold'])
            logging.info("ReadPosRankSum filtering done")
        else:
            logging.info("ReadPosRankSum information not available")
            summary.append(["ReadPosRankSum", "-", "-", "-", "Not available"])

    csv_writer(summary)


    return mt


def sample_filtering(mt, sequencingType): 
    """
    Applies sample quality control based on the thresholds stated in config.yaml
    :params: mt without sample qc
    :return: mt with sample QCed
    """
    # ensure the sequencing type entered in conf.py is valid
    assert sequencingType in ["WES", "WGS"], "sequencingType must be 'WES' or 'WGS'"

    # Compute sample QC
    mt = hl.sample_qc(mt)

    # Compute CHARR
    if config['ref_gen'] == "GRCh37": 
        if not config.get('gnomad_sites_GRCh37') or config['gnomad_sites_GRCh37'].strip() == '':
            logging.error("gnomAD sites GRCh37 path is empty or not configured - check your conf.py")
        else:
            mt = run_charr(mt, config['gnomad_sites_GRCh37'])
    elif config['ref_gen'] == "GRCh38": 
        if not config.get('gnomad_sites_GRCh38') or config['gnomad_sites_GRCh38'].strip() == '':
            logging.error("gnomAD sites GRCh83 path is empty or not configured - check your conf.py")
        else:
            mt = run_charr(mt, config['gnomad_sites_GRCh38'])       
    else:
        logging.error("The value inserted as ref_gen is not valid. Please, choose between GRCh37 and GRCh38")


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
            # get correct threshold 
            min_coverage = config['sample_filters']['DP_STATS.MEAN_WES_threshold'] if sequencingType == "WES" else config['sample_filters']['DP_STATS.MEAN_WGS_threshold']
            
            if config['verbosity']:
                stats = create_stats(mt, "sample_qc.dp_stats.mean", "cols")
                repo_stats = mt.aggregate_cols(hl.struct(
                            total = hl.agg.count(),
                            passing=hl.agg.count_where(mt.sample_qc.dp_stats.mean >= min_coverage)
                            )
                )   
            logging.info(f"Minimum Coverage filtering done - Variants removed: {repo_stats.total - repo_stats.passing}  ")
            summary.append(["Minimum Coverage", repo_stats.total, repo_stats.total - (repo_stats.total-repo_stats.passing), repo_stats.total-repo_stats.passing, stats])
            
            mt = mt.filter_cols(mt.sample_qc.dp_stats.mean >= min_coverage)
            logging.info("DP filtering done")
        
        else:
            summary.append(["Minimum Coverage", "-", "-", "-", "Not available"])
            logging.info("DP information not available")

    # ti/tv ratio filtering
    if config['sample_filters']['R_TI_TV_WES_threshold'] is not None or config['sample_filters']['R_TI_TV_WGS_threshold'] is not None:
        
        lower, upper = config['sample_filters']['R_TI_TV_WES_threshold'] if sequencingType == "WES" else config['sample_filters']['R_TI_TV_WGS_threshold']
        
        if "r_ti_tv" in mt.sample_qc:
            if config['verbosity']:
                stats = create_stats(mt, "sample_qc.r_ti_tv", "cols")
                repo_stats = mt.aggregate_cols(hl.struct(
                            total = hl.agg.count(),
                            passing=hl.agg.count_where((mt.sample_qc.r_ti_tv >= lower) & (mt.sample_qc.r_ti_tv <= upper))
                            )
                )   
            logging.info(f"Ti/Tv ratio filtering done - Variants removed: {repo_stats.total - repo_stats.passing}  ")
            summary.append(["Ti/Tv ratio", repo_stats.total, repo_stats.total - (repo_stats.total-repo_stats.passing), repo_stats.total-repo_stats.passing, stats])
            
            mt = mt.filter_cols((mt.sample_qc.r_ti_tv >= lower) & (mt.sample_qc.r_ti_tv <= upper))
            logging.info("Ti/Tv ratio filtering done")
            
        else:
            summary.append(["Ti/Tv Ratio", "-", "-", "-", "Not available"])
            logging.info("Ti/Tv information not available")

    # Call rate filtering
    if config['sample_filters']['CALL_RATE_threshold'] is not None: 
        if "call_rate" in mt.sample_qc:
            if config['verbosity']:
                stats = create_stats(mt, "sample_qc.call_rate", "cols")
                repo_stats = mt.aggregate_cols(hl.struct(
                            total = hl.agg.count(),
                            passing=hl.agg.count_where(mt.sample_qc.call_rate >= config['sample_filters']['CALL_RATE_threshold'])
                            )
                )   
            logging.info(f"Call rate filtering done - Variants removed: {repo_stats.total - repo_stats.passing}  ")
            summary.append(["Call rate", repo_stats.total, repo_stats.total - (repo_stats.total-repo_stats.passing), repo_stats.total-repo_stats.passing, stats])

            mt = mt.filter_cols(mt.sample_qc.call_rate >= config['sample_filters']['CALL_RATE_threshold'])
            logging.info("Call rate filtering done")
        
        else:
            summary.append(["Call Rate", "-", "-", "-", "Not available"])
            logging.info("Call rate information not available")


    # Singletons filtering
    if config['sample_filters']['N_SINGLETON_WES_threshold'] is not None or config['sample_filters']['N_SINGLETON_WGS_threshold'] is not None:
        if "n_singleton" in mt.sample_qc: 
            hard_cutoff = config['sample_filters']['N_SINGLETON_WES_threshold'] if sequencingType == "WES" else config['sample_filters']['N_SINGLETON_WGS_threshold']
            if config['verbosity']:
                stats = create_stats(mt, "sample_qc.n_singleton", "cols")
                repo_stats = mt.aggregate_cols(hl.struct(
                            total = hl.agg.count(),
                            passing=hl.agg.count_where(mt.sample_qc.n_singleton <= hard_cutoff)
                            )
                )   
            logging.info(f"Singletons filtering done - Variants removed: {repo_stats.total - repo_stats.passing}  ")
            summary.append(["Singletons", repo_stats.total, repo_stats.total - (repo_stats.total-repo_stats.passing), repo_stats.total-repo_stats.passing, stats])

            mt = mt.filter_cols(mt.sample_qc.n_singleton <= hard_cutoff)
            logging.info("Singleton filtering done")
        else:
            summary.append(["Singletons", "-", "-", "-", "Not available"])
            logging.info("Number of singletons not available")


    # CHARR filtering 
    
    if config['sample_filters']['CHARR_threshold'] is not None:
        try:
            if config['verbosity']:
                stats = create_stats(mt, "charr", "cols")
                repo_stats = mt.aggregate_cols(hl.struct(
                            total = hl.agg.count(),
                            passing=hl.agg.count_where(mt.charr <= config['sample_filters']['CHARR_threshold'])
                            )
                )   
            logging.info(f"CHARR filtering done - Variants removed: {repo_stats.total - repo_stats.passing}  ")
            summary.append(["CHARR", repo_stats.total, repo_stats.total - (repo_stats.total-repo_stats.passing), repo_stats.total-repo_stats.passing, stats])
                
            mt = mt.filter_cols(mt.charr <= config['sample_filters']['CHARR_threshold']) # Filter samples with charr over the threshold
            logging.info("CHARR filtering done")
        except TypeError:
            logging.info("CHARR couldn't be calculated")
            summary.append(["CHARR", "-", "-", "-", "Not available"])
        except LookupError:
            logging.info("CHARR couldn't be calculated")
            summary.append(["CHARR", "-", "-", "-", "Not available"])


    # Het/Hom ratio filtering
    if config['sample_filters']['R_HET_HOM_VAR_WES_threshold'] is not None or config['sample_filters']['R_HET_HOM_VAR_WGS_threshold'] is not None:
        if "GT" in mt.entry:
            hard_cutoff = config['sample_filters']['R_HET_HOM_VAR_WES_threshold'] if sequencingType == "WES" else config['sample_filters']['R_HET_HOM_VAR_WGS_threshold']
            if config['verbosity']:
                stats = create_stats(mt, "sample_qc.r_het_hom_var", "cols")
                repo_stats = mt.aggregate_cols(hl.struct(
                            total = hl.agg.count(),
                            passing=hl.agg.count_where(mt.sample_qc.r_het_hom_var <= hard_cutoff)
                            )
                )   
                logging.info(f"Het/Hom ratio filtering done - Variants removed: {repo_stats.total - repo_stats.passing}  ")
                summary.append(["Het/Hom ratio", repo_stats.total, repo_stats.total - (repo_stats.total-repo_stats.passing), repo_stats.total-repo_stats.passing, stats])
                
            mt = mt.filter_cols(mt.sample_qc.r_het_hom_var <= hard_cutoff)
            logging.info("Het/Hom ratio filtering done")
        else:
            summary.append(["Het/Hom Ratio", "-", "-", "-", "Not available"])
            logging.info("Het/Hom ratio information not available")

    if config['verbosity']:
        csv_writer(summary)
    
    return mt
