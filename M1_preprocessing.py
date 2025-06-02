import hail as hl
import logging
import os
from csv_utils import *
import yaml

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s]: %(message)s",
    handlers=[
        logging.FileHandler("output.log", mode='w'),
        logging.StreamHandler()
    ]
)

def load_config(path="config.yaml"):
    with open(path, "r") as file:
        return yaml.safe_load(file)

config = load_config()

def convert_and_merge_vcfs(vcf_dir, mt_from_vcf, reference_genome):

    # Get all VCF file paths in the directory
    vcfs = sorted([
        os.path.join(vcf_dir, f)
        for f in os.listdir(vcf_dir)
        if f.endswith('.vcf') or f.endswith('.vcf.bgz') or f.endswith('.vcf.gz')
    ])

    if not vcfs:
        raise ValueError("No VCF files found in the directory.")

    logging.info(f"Found {len(vcfs)} VCF files to merge.")

    if len(vcfs) == 1:
        mt = hl.import_vcf(vcfs[0], reference_genome=reference_genome, min_partitions=4)
        mt.write(mt_from_vcf, overwrite=True)

    else:
        # Load the first VCF as the initial MatrixTable
        mt = hl.import_vcf(vcfs[0], reference_genome=reference_genome, min_partitions=4)
        logging.info(f"Loaded: {vcfs[0]} -> {mt.count()}")

        # Merge in the rest
        for vcf_path in vcfs[1:]:
            new_mt = hl.import_vcf(vcf_path, reference_genome=reference_genome, min_partitions=4)
            logging.info(f"Loaded: {vcf_path} -> {new_mt.count()}")
            mt = mt.union_cols(new_mt, row_join_type='outer') # row keys which exist in only one input dataset are
            # also included

        mt.write(mt_from_vcf, overwrite=True)

    final_count = mt.count()
    logging.info(f"+++ Running EGA standard VCF workflow v1 +++") ## TODO make it work even if convert_and_merge_vcfs
    # is not run --> add it to the main()
    logging.info(f"Final MatrixTable written to {mt_from_vcf}")
    logging.info(f"Original dataset dimensions: {final_count[0]} variants and {final_count[1]} samples")

    return mt

def split_multiallelic(mt):
    summary = []
    summary.append(["", "Variants", "Samples"])
    initial_count = mt.count()
    summary.append(["Nº Variants before splitting", initial_count[0], initial_count[1]])
    logging.info(f"Before splitting multiallelic variants: {initial_count[0]} variants and {initial_count[1]} samples")

    mt = hl.split_multi_hts(mt)

    after_splitting = mt.count()
    summary.append(["Nº Variants after splitting", after_splitting[0], after_splitting[1]])
    logging.info(f"After splitting multiallelic variants:  {after_splitting[0]} variants and {after_splitting[1]} samples")
    csv_writer(summary, "QCsummary.csv")
    return mt

def gnomad_variant_filtering(mt, output_csv="filtering_summary.csv"):
    mt = hl.variant_qc(mt)
    initial_count = mt.count()
    logging.info(f"Original number of variants: {initial_count[0]}")
    # Store filtering summary
    summary = []
    summary.append([])
    summary.append(["Variant Filtering Steps"])
    summary.append(["Filter", "Before", "After", "Removed"])

    # Filtering steps
    if "QD" in mt.info:
        pre_count = mt.count_rows()
        mt = mt.filter_rows(mt.info.QD >= config['QD_threshold'])
        post_count = mt.count_rows()
        removed = pre_count - post_count
        logging.info(f"QD filtering done - Variants removed: {removed}")
        summary.append(["QD", pre_count, post_count, removed])
    else:
        logging.info("QD information not available")
        summary.append(["QD", "-", "-", "-", "Not available"])

    if "DP" in mt.info:
        pre_count = mt.count_rows()
        mt = mt.filter_rows(mt.info.DP >= config['DP_threshold'])
        post_count = mt.count_rows()
        removed = pre_count - post_count
        logging.info(f"DP filtering done - Variants removed: {removed}")
        summary.append(["DP", pre_count, post_count, removed])
    else:
        logging.info("DP information not available")
        summary.append(["DP", "-", "-", "-", "Not available"])

    if mt.qual.dtype == hl.tfloat64:
        pre_count = mt.count_rows()
        mt = mt.filter_rows(mt.qual >= config['QUAL_threshold'])
        post_count = mt.count_rows()
        removed = pre_count - post_count
        logging.info(f"QUAL filtering done - Variants removed: {removed}")
        summary.append(["QUAL", pre_count, post_count, removed])
    else:
        logging.info("QUAL information not available")
        summary.append(["QUAL", "-", "-", "-", "Not available"])

    if "MQ" in mt.info:
        pre_count = mt.count_rows()
        mt = mt.filter_rows(mt.info.MQ >= config['MQ_threshold'])
        post_count = mt.count_rows()
        removed = pre_count - post_count
        logging.info(f"MQ filtering done - Variants removed: {removed}")
        summary.append(["MQ", pre_count, post_count, removed])
    else:
        logging.info("MQ information not available")
        summary.append(["MQ", "-", "-", "-", "Not available"])

    if "FS" in mt.info:
        pre_count = mt.count_rows()
        mt = mt.filter_rows(mt.info.FS <= config['FS_threshold'])
        post_count = mt.count_rows()
        removed = pre_count - post_count
        logging.info(f"FS filtering done - Variants removed: {removed}")
        summary.append(["FS", pre_count, post_count, removed])
    else:
        logging.info("FS information not available")
        summary.append(["FS", "-", "-", "-", "Not available"])

    if "ReadPosRankSum" in mt.info:
        pre_count = mt.count_rows()
        mt = mt.filter_rows(mt.info.ReadPosRankSum >= config['ReadPosRankSum_threshold'])
        post_count = mt.count_rows()
        removed = pre_count - post_count
        logging.info(f"ReadPosRankSum filtering done - Variants removed: {removed}")
        summary.append(["ReadPosRankSum", pre_count, post_count, removed])
    else:
        logging.info("ReadPosRankSum information not available")
        summary.append(["ReadPosRankSum", "-", "-", "-", "Not available"])

    if "GQ" in mt.entry:
        pre_count = mt.count()
        mt = mt.filter_entries(mt.GQ >= config['GQ_threshold'])
        post_count = mt.count()
        removed = (pre_count[0] * pre_count[1]) - (post_count[0] * post_count[1]) ## TODO check genotype count is
        # correct
        logging.info(f"GQ filtering done - Genotypes removed: {removed}")
        summary.append(["GQ (Genotypes)", pre_count[0] * pre_count[1], post_count[0] * post_count[1], removed])
    else:
        logging.info("GQ information not available")
        summary.append(["GQ (Genotypes)", "-", "-", "-", "Not available"])

    if "GT" in mt.entry and "AD" in mt.entry:
        pre_count = mt.count()
        mt = mt.annotate_entries(
            AB=hl.or_missing(
                mt.GT.is_het(),
                mt.AD[1] / hl.sum(mt.AD)
            )
        )
        mt = mt.filter_entries(
            hl.if_else(
                mt.GT.is_het(),
                mt.AB >= config['AB_threshold'],
                True
            )
        )
        post_count = mt.count()
        removed = (pre_count[0] * pre_count[1]) - (post_count[0] * post_count[1])
        logging.info(f"Allele Balance filtering done - Genotypes removed: {removed}")
        summary.append(
            ["AlleleBalance (Genotypes)", pre_count[0] * pre_count[1], post_count[0] * post_count[1], removed])
    else:
        logging.info("Allele Balance information not available")
        summary.append(["AlleleBalance (Genotypes)", "-", "-", "-", "Not available"])

    final_count = mt.count()
    logging.info(f"Final number of variants after QC: {final_count[0]}")

    csv_writer(summary, output_csv="QCsummary.csv")

    return mt


def gnomad_sample_filtering(mt, sequencingType):
    assert sequencingType in ["WES", "WGS"], "sequencingType must be 'WES' or 'WGS'"

    mt = hl.sample_qc(mt)
    initial_count = mt.count()
    summary = []
    summary.append([])
    summary.append(["Sample Filtering Steps"])
    summary.append(["Filter", "Before", "After", "Removed"])

    # Minimum coverage filtering
    if "dp_stats" in mt.sample_qc:
        pre_count = mt.count()
        min_coverage = config['DP_WES_threshold'] if sequencingType == "WES" else config['DP_WGS_threshold']
        mt = mt.filter_cols(mt.sample_qc.dp_stats.mean >= min_coverage)
        post_count = mt.count()
        summary.append(["Minimum Coverage", pre_count[1], post_count[1], pre_count[1] - post_count[1]])
        logging.info(f"{sequencingType} coverage filtering done — Samples removed: {pre_count[1] - post_count[1]}")
    else:
        summary.append(["Minimum Coverage", "-", "-", "-", "Not available"])

    # ti/tv ratio filtering
    if "r_ti_tv" in mt.sample_qc:
        pre_count = mt.count()
        lower, upper = config['TITV_WES_threshold'] if sequencingType == "WES" else config['TITV_WGS_threshold']
        mt = mt.filter_cols((mt.sample_qc.r_ti_tv >= lower) & (mt.sample_qc.r_ti_tv <= upper))
        post_count = mt.count()
        summary.append(["Ti/Tv Ratio", pre_count[1], post_count[1], pre_count[1] - post_count[1]])
        logging.info(f"{sequencingType} ti/tv filtering done — Samples removed: {pre_count[1] - post_count[1]}")
    else:
        summary.append(["Ti/Tv Ratio", "-", "-", "-", "Not available"])

    # Call rate filtering
    if "call_rate" in mt.sample_qc:
        pre_count = mt.count()
        mt = mt.filter_cols(mt.sample_qc.call_rate >= config['callRate_threshold'])
        post_count = mt.count()
        summary.append(["Call Rate", pre_count[1], post_count[1], pre_count[1] - post_count[1]])
        logging.info(f"Call rate filtering done — Samples removed: {pre_count[1] - post_count[1]}")
    else:
        summary.append(["Call Rate", "-", "-", "-", "Not available"])

    # Singletons filtering
    if "n_singleton" in mt.sample_qc: ## TODO make it work with relative and hard thresholds -- now only samples with
        # both thresholds are retained, check with abeer
        pre_count = mt.count()
        singleton_stats = mt.aggregate_cols(hl.agg.stats(mt.sample_qc.n_singleton))
        max_singletons_stdev = 2 * singleton_stats.stdev # relative threshold
        hard_cutoff = config['singletons_WES_threshold'] if sequencingType == "WES" else config['singletons_WGS_threshold']
        mt = mt.filter_cols(
            (mt.sample_qc.n_singleton <= max_singletons_stdev) & (mt.sample_qc.n_singleton <= hard_cutoff))
        post_count = mt.count()
        summary.append(["Singletons", pre_count[1], post_count[1], pre_count[1] - post_count[1]])
        logging.info(f"Singleton filtering done — Samples removed: {pre_count[1] - post_count[1]}")
    else:
        summary.append(["Singletons", "-", "-", "-", "Not available"])

    # CHARR filtering — placeholder
    summary.append(["CHARR", "-", "-", "-", "TODO"])

    # Het/Hom ratio filtering
    if "GT" in mt.entry:
        pre_count = mt.count()
        mt = mt.annotate_cols(
            n_het=hl.agg.count_where(mt.GT.is_het()),
            n_hom=hl.agg.count_where(mt.GT.is_hom_var()),
        )
        mt = mt.annotate_cols(
            het_hom_ratio=hl.float(mt.n_het) / hl.if_else(mt.n_hom == 0, 1.0, hl.float(mt.n_hom))
        )
        max_het_hom = config['hethom_WES_threshold'] if sequencingType == "WES" else config['hethom_WGS_threshold']
        mt = mt.filter_cols(mt.het_hom_ratio <= max_het_hom)
        post_count = mt.count()
        summary.append(["Het/Hom Ratio", pre_count[1], post_count[1], pre_count[1] - post_count[1]])

        logging.info(f"{sequencingType} het/hom filtering done — Samples removed: {pre_count[1] - post_count[1]}")
    else:
        summary.append(["Het/Hom Ratio", "-", "-", "-", "Not available"])

    # Final summary
    final_count = mt.count()
    logging.info(f"Final number of samples: {final_count[1]} (removed {initial_count[1] - final_count[1]})")

    # Write CSV
    csv_writer(summary, output_csv="QCsummary.csv")

    return mt







