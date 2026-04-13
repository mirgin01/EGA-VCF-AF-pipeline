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
        mt = hl.import_vcf(f"{config['vcf_dir']}/*.vcf*", reference_genome=reference_genome, min_partitions=4)
        
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
        logging.info(f"Number of multiallelic variants: {after_splitting[0] - original_size[0]}   ")
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
            logging.error("GQ information not available - terminating pipeline")
            summary.append(["Filter", "Before", "After", "Removed", "Status"])
            summary.append(["GQ", "-", "-", "-", "Not available"])
            raise ValueError("GQ information is not available. To proceed with filtering, add a hash sign to the threshold value to inactivate the metric and perform quality control without this step.")

    else: 
        logging.warning("GQ threshold not set - GQ filtering not performed ")

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
            logging.error("Allele Balance information not available - terminating pipeline")
            summary.append(["Filter", "Before", "After", "Removed", "Status"])
            summary.append(["AlleleBalance (Genotypes)", "-", "-", "-", "Not available"])
            raise ValueError("AB information is not available. To proceed with filtering, add a hash sign to the threshold value to inactivate the metric and perform quality control without this step.")
    else: 
        logging.warning("AB threshold not set - AB filtering not performed ")
      
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
            logging.error("QD information not available - terminating pipeline")
            summary.append(["QD", "-", "-", "-", "Not available"])
            raise ValueError("QD information is not available. To proceed with filtering, add a hash sign to the threshold value to inactivate the metric and perform quality control without this step.")
    else:
        logging.warning("QD threshold not set - QD filtering not performed ")
    
    if config['variant_filters']['DP_threshold'] is not None:
        
        if "DP" in mt.info:
            
            if config['verbosity']:
                stats = create_stats(mt, "info.DP", "rows")
                summary = verbosity_counts_variants(mt, "DP", 
                                                    (hl.is_missing(mt.info.DP)) | (mt.info.DP >= config['variant_filters']['DP_threshold']), 
                                                    summary, 
                                                    stats)

            mt = mt.filter_rows(hl.is_missing(mt.info.DP) | (mt.info.DP >= config['variant_filters']['DP_threshold']))     

            logging.info("DP filtering done")   
            
        else:
            logging.error("DP information not available - terminating pipeline")
            summary.append(["DP", "-", "-", "-", "Not available"])
            raise ValueError("DP information is not available. To proceed with filtering, guradd a hash sign to the threshold value to inactivate the metrication and perform quality control without this step.")

    else: 
        logging.warning("DP threshold not set - DP filtering not performed ")

    if config['variant_filters']['QUAL_threshold'] is not None:

        try:
            if config['verbosity']:

                stats = create_stats(mt, "qual", "rows")
                summary = verbosity_counts_variants(mt, "QUAL", mt.qual >= config['variant_filters']['QUAL_threshold'], summary, stats)
            
            
            mt = mt.filter_rows(mt.qual >= config['variant_filters']['QUAL_threshold'])
            logging.info("QUAL filtering done")
        
        except TypeError:
            logging.error("QUAL information not available - terminating pipeline")
            summary.append(["QUAL", "-", "-", "-", "Not available"])
            raise ValueError("QUAL information is not available. To proceed with filtering, add a hash sign to the threshold value to inactivate the metric and perform quality control without this step.")
    else: 
        logging.warning("QUAL threshold not set - QUAL filtering not performed ")

    if config['variant_filters']['MQ_threshold'] is not None:
        if "MQ" in mt.info:
            if config['verbosity']:
                stats = create_stats(mt, "info.MQ", "rows")
                summary = verbosity_counts_variants(mt, "MQ", mt.info.MQ >= config['variant_filters']['MQ_threshold'], summary, stats)
               
            mt = mt.filter_rows(mt.info.MQ >= config['variant_filters']['MQ_threshold'])
            logging.info("MQ filtering done")
        else:
            logging.error("MQ information not available - terminating pipeline")
            summary.append(["MQ", "-", "-", "-", "Not available"])
            raise ValueError("MQ information is not available. To proceed with filtering, add a hash sign to the threshold value to inactivate the metric and perform quality control without this step.")
    else: 
        logging.warning("MQ threshold not set - MQ filtering not performed ")

    if config['variant_filters']['FS_threshold'] is not None:
        if "FS" in mt.info:
            if config['verbosity']:
                stats = create_stats(mt, "info.FS", "rows")
                summary = verbosity_counts_variants(mt, "FS", mt.info.FS <= config['variant_filters']['FS_threshold'], summary, stats)

            mt = mt.filter_rows(mt.info.FS <= config['variant_filters']['FS_threshold'])
            logging.info("FS filtering done")
        else:
            logging.error("FS information not available - terminating pipeline")
            summary.append(["FS", "-", "-", "-", "Not available"])
            raise ValueError("FS information is not available. To proceed with filtering, add a hash sign to the threshold value to inactivate the metric and perform quality control without this step.")
    else: 
        logging.warning("FS threshold not set - FS filtering not performed ")

    if config['variant_filters']['READPOSRANKSUM_threshold'] is not None:
        
        filter_expr = (hl.is_missing(mt.info.ReadPosRankSum)) | (mt.info.ReadPosRankSum >= config['variant_filters']['READPOSRANKSUM_threshold'])

        logging.info(mt.info.ReadPosRankSum.dtype)
        logging.info(mt.info.describe())

        if "ReadPosRankSum" in mt.info:
            if config['verbosity']:
                stats = create_stats(mt, "info.ReadPosRankSum", "rows")
                summary = verbosity_counts_variants(mt, "ReadPosRankSum", mt.info.ReadPosRankSum >= config['variant_filters']['READPOSRANKSUM_threshold'], summary, stats)
            
            mt = mt.filter_rows(filter_expr)

            logging.info("ReadPosRankSum filtering done")
        else:
            logging.error("ReadPosRankSum information not available - terminating pipeline")
            summary.append(["ReadPosRankSum", "-", "-", "-", "Not available"])
            raise ValueError("ReadPosRankSum information is not available. To proceed with filtering, add a hash sign to the threshold value to inactivate the metric and perform quality control without this step.")
    
    else: 
        logging.warning("READPOSRANKSUM threshold not set - READPOSRANKSUM filtering not performed ")

    csv_writer(summary)


    return mt


def sample_filtering_hard_thresholds(mt, sequencingType, basename): 
    """
    Applies sample quality control based on the thresholds stated in config.yaml
    :params: mt without sample qc
    :return: mt with sample QCed
    """
    
    # ensure the sequencing type entered in conf.py is valid
    assert sequencingType in ["WES", "WGS"], "sequencingType must be 'WES' or 'WGS'"

    # Compute sample QC
    mt = hl.sample_qc(mt)

    
    # Open summary
    summary = [] # always creates to avoid errors in case verbosity is on but no sample filters are applied
    if config['verbosity']:
        summary.append([])
        summary.append(["Sample Filtering Steps"])
        summary.append(["Filter", "Before", "After", "Removed", "Stats"])

    # Compute CHARR
    try:
        if config['ref_gen'] == "GRCh37": 
            if not config.get('gnomad_sites_GRCh37') or config['gnomad_sites_GRCh37'].strip() == '':
                logging.error("gnomAD sites GRCh37 path is empty or not configured - check your conf.py")
            else:
                mt = run_charr(mt, config['gnomad_sites_GRCh37'])
        elif config['ref_gen'] == "GRCh38": 
            if not config.get('gnomad_sites_GRCh38') or config['gnomad_sites_GRCh38'].strip() == '':
                logging.error("gnomAD sites GRCh38 path is empty or not configured - check your conf.py")
            else:
                mt = run_charr(mt, config['gnomad_sites_GRCh38'])       
        else:
            logging.error("The value inserted as ref_gen is not valid. Please, choose between GRCh37 and GRCh38")
    except ValueError as e: 
        logging.error("CHARR couldn't be calculated:")
        logging.error("%r", e)
        summary.append("CHARR couldn't be calculated.", e)
        raise ValueError("CHARR information is not available. To proceed with filtering, add a hash sign to the threshold value to inactivate the metric and perform quality control without this step.")
    
    config['verbosity'] = True

    if config["plots"]:
        logging.info("Creating Sample Filters plots")
        create_distribution_plots(mt, sequencingType, basename)


    # Minimum coverage filtering
    if config['sample_filters']['DP_STATS.MEAN_WES_threshold'] is not None or config['sample_filters']['DP_STATS.MEAN_WGS_threshold'] is not None:
        if "dp_stats" in mt.sample_qc:
            # get correct threshold 
            min_coverage = config['sample_filters']['DP_STATS.MEAN_WES_threshold'] if sequencingType == "WES" else config['sample_filters']['DP_STATS.MEAN_WGS_threshold']
            mt = _apply_filter(mt,
                                condition=mt.sample_qc.dp_stats.mean >= min_coverage,
                                metric_path="sample_qc.dp_stats.mean",
                                metric_name="Minimum Coverage",
                                summary=summary)
        
        else:
            summary.append(["Minimum Coverage", "-", "-", "-", "Not available"])
            logging.error("DP information not available - terminating pipeline")
            raise ValueError("DP information is not available. To proceed with filtering, add a hash sign to the threshold value to inactivate the metric and perform quality control without this step.")
    else:
        logging.warning("DP threshold not set - DP filtering not performed ")

    
    # ti/tv ratio filtering
    """if config['sample_filters']['R_TI_TV_WES_threshold'] is not None or config['sample_filters']['R_TI_TV_WGS_threshold'] is not None:
        
        lower, upper = config['sample_filters']['R_TI_TV_WES_threshold'] if sequencingType == "WES" else config['sample_filters']['R_TI_TV_WGS_threshold']
        
        if "r_ti_tv" in mt.sample_qc:
            mt = _apply_filter(mt,
                condition=(mt.sample_qc.r_ti_tv >= lower) & (mt.sample_qc.r_ti_tv <= upper),
                metric_path="sample_qc.r_ti_tv",
                metric_name="Ti/Tv ratio",
                summary=summary)        
        else:
            summary.append(["Ti/Tv Ratio", "-", "-", "-", "Not available"])
            logging.error("Ti/Tv information not available - terminating pipeline")
            raise ValueError("Ti/Tv ratio information is not available. To proceed with filtering, add a hash sign to the threshold value to inactivate the metric and perform quality control without this step.")
    else: 
        logging.warning("Ti/Tv ratio threshold not set - Ti/Tv filtering not performed ")"""

    # Call rate filtering
    if config['sample_filters']['CALL_RATE_threshold'] is not None: 
        if "call_rate" in mt.sample_qc:
           mt = _apply_filter(mt,
                            condition=mt.sample_qc.call_rate >= config['sample_filters']['CALL_RATE_threshold'],
                            metric_path="sample_qc.call_rate",
                            metric_name="Call Rate",
                            summary=summary)
        
        else:
            summary.append(["Call Rate", "-", "-", "-", "Not available"])
            logging.error("Call rate information not available - terminating pipeline")
            raise ValueError("Call rate information is not available. To proceed with filtering, add a hash sign to the threshold value to inactivate the metric and perform quality control without this step.")

    else:
        logging.warning("Call rate threshold not set - Call rate filtering not performed ")


    # Singletons filtering
    if config['sample_filters']['N_SINGLETON_WES_threshold'] is not None or config['sample_filters']['N_SINGLETON_WGS_threshold'] is not None:
        if "n_singleton" in mt.sample_qc: 
            hard_cutoff = config['sample_filters']['N_SINGLETON_WES_threshold'] if sequencingType == "WES" else config['sample_filters']['N_SINGLETON_WGS_threshold']
            mt = _apply_filter(mt,
                                condition=mt.sample_qc.n_singleton <= hard_cutoff,
                                metric_path="sample_qc.n_singleton",
                                metric_name="Singletons",
                                summary=summary)
        else:
            summary.append(["Singletons", "-", "-", "-", "Not available"])
            logging.error("Number of singletons not available - terminating pipeline")
            raise ValueError("Singletons information is not available. To proceed with filtering, add a hash sign to the threshold value to inactivate the metric and perform quality control without this step.")
    else: 
        logging.warning("Number of singletons threshold not set - Singleton filtering not performed ")


    # CHARR filtering 
    
    if config['sample_filters']['CHARR_threshold'] is not None:
        try:
            mt = _apply_filter(mt,
                                condition=mt.charr <= config['sample_filters']['CHARR_threshold'],
                                metric_path="charr",
                                metric_name="CHARR",
                                summary=summary)
        except (TypeError, LookupError, ValueError) as e:
            logging.error(f"CHARR couldn't be calculated - terminating pipeline: {e}")
            if config['verbosity']:
                summary.append(["CHARR", "-", "-", "-", "Not available"])
            raise ValueError("CHARR information is not available. To proceed with filtering, add a hash sign to the threshold value to inactivate the metric and perform quality control without this step.")
    else:
        logging.warning("CHARR threshold not set - CHARR filtering not performed ")


    # Het/Hom ratio filtering
    """if config['sample_filters']['R_HET_HOM_VAR_WES_threshold'] is not None or config['sample_filters']['R_HET_HOM_VAR_WGS_threshold'] is not None:
        if "r_het_hom_var" in mt.sample_qc:
            hard_cutoff = config['sample_filters']['R_HET_HOM_VAR_WES_threshold'] if sequencingType == "WES" else config['sample_filters']['R_HET_HOM_VAR_WGS_threshold']
            
            mt = _apply_filter(mt,
                                condition=mt.sample_qc.r_het_hom_var <= hard_cutoff,
                                metric_path="sample_qc.r_het_hom_var",
                                metric_name="Het/Hom ratio",
                                summary=summary)
        else:
            summary.append(["Het/Hom Ratio", "-", "-", "-", "Not available"])
            logging.error("Het/Hom ratio information not available - terminating pipeline")
            raise ValueError("Het/Hom information is not available. To proceed with filtering, add a hash sign to the threshold value to inactivate the metric and perform quality control without this step.")
    else:
        logging.warning("Het/hom ratio threshold not set - Het/hom ratio filtering not performed")"""

    if config['verbosity']:
        csv_writer(summary)
    
    return mt


def sample_filtering_mad_thresholds(mt, sequencingType, basename):
    """
    Applies sample quality control using MAD-based thresholds following the gnomAD strategy.
    Thresholds are computed from the full dataset before any filtering is applied.
    :params: mt without sample qc
    :return: mt with sample QCed
    """

    assert sequencingType in ["WES", "WGS"], "sequencingType must be 'WES' or 'WGS'"

    # Compute sample QC
    mt = hl.sample_qc(mt)

    # Open summary
    summary = []
    if config['verbosity']:
        summary.append([])
        summary.append(["Sample Filtering Steps — MAD-based thresholds"])
        summary.append(["Filter", "Before", "After", "Removed", "Stats"])

    # Compute CHARR
    """try:
        if config['ref_gen'] == "GRCh37":
            if not config.get('gnomad_sites_GRCh37') or config['gnomad_sites_GRCh37'].strip() == '':
                logging.error("gnomAD sites GRCh37 path is empty or not configured - check your conf.py")
            else:
                mt = run_charr(mt, config['gnomad_sites_GRCh37'])
        elif config['ref_gen'] == "GRCh38":
            if not config.get('gnomad_sites_GRCh38') or config['gnomad_sites_GRCh38'].strip() == '':
                logging.error("gnomAD sites GRCh38 path is empty or not configured - check your conf.py")
            else:
                mt = run_charr(mt, config['gnomad_sites_GRCh38'])
        else:
            logging.error("The value inserted as ref_gen is not valid. Please, choose between GRCh37 and GRCh38")
    except (TypeError, LookupError, ValueError) as e:
        logging.error(f"CHARR couldn't be calculated: {e}")
        if config['verbosity']:
            summary.append(["CHARR", "-", "-", "-", "Not available"])
        raise ValueError("CHARR information is not available. To proceed with filtering, add a hash sign to the threshold value to inactivate the metric and perform quality control without this step.")
    """

    # ── collect all metrics in one Hail action ────────────────────────────────
    # must happen after sample_qc and charr, before any filtering
    # so all thresholds reflect the full dataset
    logging.info("start collecting metrics")

    collected = _collect_metrics(mt)
    
    logging.info("finished collecting metrics", collected)

    if config["plots"]:
        logging.info("Creating Sample Filters plots")
        sample_ids = create_distribution_plots(mt, sequencingType, basename)
        #create_outlier_plots(mt, z_thresh=3.0, name_prefix=get_name(),collected=collected, sample_ids=sample_ids)

    # ── compute all MAD thresholds upfront ────────────────────────────────────
    thresholds = {
        'dp_stats_mean': _safe_mad_threshold(collected, 'dp_stats_mean', k=4.0, one_sided='lower'),
        'r_ti_tv':       _safe_mad_threshold(collected, 'r_ti_tv',       k=4.0, one_sided='lower'),
        'call_rate':     _safe_mad_threshold(collected, 'call_rate',     k=4.0, one_sided='lower'),
        'n_singleton':   _safe_mad_threshold(collected, 'n_singleton',   k=4.0, one_sided='upper'),
        'charr':         _safe_mad_threshold(collected, 'charr',         k=4.0, one_sided='upper'),
        'r_het_hom_var': _safe_mad_threshold(collected, 'r_het_hom_var', k=4.0, one_sided=None),
    }

    # ── apply filters ─────────────────────────────────────────────────────────

    # Minimum coverage filtering
    """if "dp_stats" in mt.sample_qc:
        mt = _apply_mad_filter(mt,
                               metric_key='dp_stats_mean',
                               hail_expr=mt.sample_qc.dp_stats.mean,
                               metric_path="sample_qc.dp_stats.mean",
                               metric_name="Minimum Coverage",
                               summary=summary,
                               thresholds=thresholds)
    else:
        logging.warning("DP information not available - DP filtering not performed")"""

    # Ti/Tv ratio filtering
    if "r_ti_tv" in mt.sample_qc:
        mt = _apply_mad_filter(mt,
                               metric_key='r_ti_tv',
                               hail_expr=mt.sample_qc.r_ti_tv,
                               metric_path="sample_qc.r_ti_tv",
                               metric_name="Ti/Tv ratio",
                               summary=summary,
                               thresholds=thresholds)
    else:
        logging.warning("Ti/Tv information not available - Ti/Tv filtering not performed")

    # Call rate filtering
    """if "call_rate" in mt.sample_qc:
        mt = _apply_mad_filter(mt,
                               metric_key='call_rate',
                               hail_expr=mt.sample_qc.call_rate,
                               metric_path="sample_qc.call_rate",
                               metric_name="Call Rate",
                               summary=summary,
                               thresholds=thresholds)
    else:
        logging.warning("Call rate information not available - Call rate filtering not performed")

    # Singletons filtering
    if "n_singleton" in mt.sample_qc:
        mt = _apply_mad_filter(mt,
                               metric_key='n_singleton',
                               hail_expr=mt.sample_qc.n_singleton,
                               metric_path="sample_qc.n_singleton",
                               metric_name="Singletons",
                               summary=summary,
                               thresholds=thresholds)
    else:
        logging.warning("Singletons information not available - Singletons filtering not performed")

    # CHARR filtering
    if 'charr' in mt.col:
        mt = _apply_mad_filter(mt,
                               metric_key='charr',
                               hail_expr=mt.charr,
                               metric_path="charr",
                               metric_name="CHARR",
                               summary=summary,
                               thresholds=thresholds)
    else:
        logging.warning("CHARR not available - CHARR filtering not performed")"""

    # Het/Hom ratio filtering
    if "r_het_hom_var" in mt.sample_qc:
        mt = _apply_mad_filter(mt,
                               metric_key='r_het_hom_var',
                               hail_expr=mt.sample_qc.r_het_hom_var,
                               metric_path="sample_qc.r_het_hom_var",
                               metric_name="Het/Hom ratio",
                               summary=summary,
                               thresholds=thresholds)
    else:
        logging.warning("Het/Hom information not available - Het/Hom filtering not performed")

    if config['verbosity']:
        csv_writer(summary)

    return mt


def _apply_filter(mt, condition, metric_path, metric_name, summary):
    """
    Applies a single sample QC filter step.
    :param mt: MatrixTable to filter
    :param condition: Hail boolean expression for filtering
    :param metric_path: string path to metric for create_stats e.g. "sample_qc.call_rate"
    :param metric_name: display name for logging and summary e.g. "Call Rate"
    :param summary: summary list to append results to
    :return: filtered MatrixTable
    """
    repo_stats = mt.aggregate_cols(hl.struct(
        total=hl.agg.count(),
        passing=hl.agg.count_where(condition)
    ))

    if config['verbosity']:
        stats = create_stats(mt, metric_path, "cols")
        summary.append([metric_name, repo_stats.total, repo_stats.passing,
                        repo_stats.total - repo_stats.passing, stats])

    logging.info(f"{metric_name} filtering done - Samples removed: {repo_stats.total - repo_stats.passing}")
    mt = mt.filter_cols(condition)
    return mt


def _compute_mad_threshold(values, k=4.0, one_sided=None):
    """
    Computes MAD-based thresholds for a metric.
    :param values: numpy array of metric values
    :param k: MAD multiplier (gnomAD default = 4.0)
    :param one_sided: None (two-sided), 'lower' (flag below only), 'upper' (flag above only)
    :return: (lower, upper, median, mad)
    """
    median = float(np.median(values))
    mad    = float(np.median(np.abs(values - median)))

    # floor to avoid threshold collapse when MAD = 0
    if mad == 0:
        mad = 0.01 * abs(median) if median != 0 else 0.01

    spread = k * mad
    lower  = median - spread if one_sided != 'upper' else -np.inf
    upper  = median + spread if one_sided != 'lower' else  np.inf

    logging.info(f"MAD threshold | median={median:.4f}, MAD={mad:.4f}, "
                 f"lower={lower:.4f}, upper={upper:.4f}")

    return lower, upper, median, mad


def _collect_metrics(mt):
    """
    Collects all QC metrics from the MatrixTable to numpy arrays in one Hail action.
    :param mt: MatrixTable after hl.sample_qc() and run_charr()
    :return: dictionary of {metric_name: numpy array}
    """
    has_charr = 'charr' in mt.col
    extra     = ['charr'] if has_charr else []
    rows      = mt.cols().key_by().select('s', 'sample_qc', *extra).collect()

    extractors = {
        'call_rate':     lambda r: r.sample_qc.call_rate,
        'r_ti_tv':       lambda r: r.sample_qc.r_ti_tv,
        'r_het_hom_var': lambda r: r.sample_qc.r_het_hom_var,
        'n_singleton':   lambda r: r.sample_qc.n_singleton,
        'dp_stats_mean': lambda r: r.sample_qc.dp_stats.mean,
        'charr':         lambda r: r.charr if has_charr else None,
    }

    collected = {}
    for key, extractor in extractors.items():
        try:
            vals = [float(extractor(r)) if extractor(r) is not None else np.nan for r in rows]
            collected[key] = np.array(vals, dtype=float)
            n_valid = int(np.sum(~np.isnan(collected[key])))
            logging.info(f"Collected '{key}': {n_valid}/{len(rows)} valid values")
        except Exception as e:
            logging.warning(f"Could not collect metric '{key}': {e}")
            collected[key] = np.full(len(rows), np.nan)

    return collected

def _apply_mad_filter(mt, metric_key, hail_expr, metric_path, metric_name,
                      summary, thresholds):
    """
    Applies a pre-computed MAD-based filter to the MatrixTable.
    :param mt: MatrixTable to filter
    :param metric_key: key to look up in thresholds dict e.g. 'call_rate'
    :param hail_expr: Hail expression for the metric e.g. mt.sample_qc.call_rate
    :param metric_path: string path for create_stats e.g. "sample_qc.call_rate"
    :param metric_name: display name for logging and summary
    :param summary: summary list to append results to
    :param thresholds: dict of pre-computed thresholds from _compute_mad_threshold
    :return: filtered MatrixTable
    """
    # Avoids breaking when metric is not in matrix 
    if thresholds[metric_key] is None:
        logging.warning(f"{metric_name} MAD filter skipped — threshold could not be computed")
        return mt                  # return mt unchanged
    
    lower, upper, median, mad = thresholds[metric_key]

    condition = (hail_expr >= lower) & (hail_expr <= upper)

    logging.info(f"{metric_name} MAD filter | lower={lower:.4f}, upper={upper:.4f}")

    mt = _apply_filter(mt,
                       condition=condition,
                       metric_path=metric_path,
                       metric_name=metric_name,
                       summary=summary)
    return mt


def _safe_mad_threshold(collected, key, k=4.0, one_sided=None):
    """
    Wraps _compute_mad_threshold with a NaN check.
    Returns None if the metric is unavailable, threshold tuple otherwise.
    """
    values = collected.get(key)
    if values is None:
        logging.warning(f"MAD | '{key}' not found in collected metrics — skipping")
        return None
    
    valid = values[~np.isnan(values)]
    if len(valid) < 2:
        logging.warning(f"MAD | '{key}' has fewer than 2 valid values — skipping")
        return None

    return _compute_mad_threshold(valid, k=k, one_sided=one_sided)