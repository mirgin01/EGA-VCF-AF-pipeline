import logging
import hail as hl
import pandas as pd
from csv_utils import *
import yaml
import csv
import os

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

def subset_matrix(mt):
    # load matrix
    mt = mt.key_rows_by('locus', 'alleles')

    # load ancestry SNPs depending on the ref_gen
    ht_variants = hl.import_table(f"{config['ancestrySNPs']}ancestry_SNPs_GrafAnc_{config['ref_gen']}",
                                  delimiter='\t', impute=True, skip_blank_lines = True)

    ht_variants = ht_variants.key_by(**hl.parse_variant(ht_variants.ancestrySNPs)) #Convert variants in string
    # format to separate locus and allele hail fields
    ht_variants = ht_variants.key_by('locus', 'alleles')

    # filter variants
    ancestry_mt = mt.filter_rows(hl.is_defined(ht_variants[mt.locus, mt.alleles])) # filter_rows(to_keep)
    ancestry_counts = ancestry_mt.count()

    # apply HWE filter
    ancestry_mt = hl.variant_qc(ancestry_mt)
    pre_count = ancestry_mt.count_rows()
    ancestry_mt = ancestry_mt.filter_rows(ancestry_mt.variant_qc.p_value_hwe > 1e-6)
    post_count = ancestry_mt.count_rows()
    logging.info(f"HWE filtering done - Variants removed: {pre_count - post_count}")

    # export VCF
    print(ancestry_mt.count())
    mt_path = config['mt_afterQC'] if config['preprocessing'] else config['mt_from_vcf']
    logging.info(f"There are {pre_count} ancestry SNPs in {mt_path}")
    ancestry_vcf = os.path.join(os.path.dirname(mt_path), "ancestrySNPs.vcf")
    logging.info(f"Exporting ancestry matrix to {ancestry_vcf}")
    hl.export_vcf(ancestry_mt, ancestry_vcf)

    summary = []
    summary.append([])
    summary.append(["Ancestry filtering"])
    summary.append(["Nº of ancestry SNPs", ancestry_counts[0]])
    summary.append([])
    summary.append(["", "Before", "After", "Removed"])
    summary.append(["HWE equilbrium", pre_count, post_count, pre_count-post_count])
    csv_writer(summary, output_csv="QCsummary.csv")

    return ancestry_vcf

def call_grafanc(ancestry_vcf, mt_path):
    # run GrafAnc
    basename = os.path.basename(config['mt_afterQC']) # create name for results file
    name_only = os.path.splitext(basename)[0]
    ancestry_results = os.path.join(os.path.dirname(mt_path), f"{name_only}-GrafAnc_results")
    logging.info(f"Running GrafAnc: $ grafanc {ancestry_vcf} {ancestry_results}")
    os.system(f"grafanc {ancestry_vcf} {ancestry_results} --threads 4'")

    # update ancestry groups to super populations
    df = pd.read_csv(ancestry_results, sep='\t')

    ancestry_map = {  # Super population mapping based on the first digit
        '1': 'African',
        '2': 'Middle East and North Africa',
        '3': 'European',
        '4': 'South Asian',
        '5': 'East Asian',
        '6': 'American',
        '7': 'Oceania',
        '8': 'Multi-ancestry',
    }

    df['AncGroupNAME'] = df['AncGroupID'].astype(str).str[0].map(ancestry_map)
    df.to_csv(ancestry_results, sep="\t", index=False) # save new column to GrafAnc results

    return ancestry_results

def annotate_ancestry(ancestry_results, mt):
    ancestry_table = hl.import_table(ancestry_results, delimiter="\t", impute=True)
    ancestry_table = ancestry_table.select('Sample', 'AncGroupNAME')
    ancestry_table = ancestry_table.key_by('Sample')
    mt = mt.annotate_cols(ancestry=ancestry_table[mt.s].AncGroupNAME)

    return mt

