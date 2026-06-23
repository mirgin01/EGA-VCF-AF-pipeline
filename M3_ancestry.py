import logging
import hail as hl
import pandas as pd
from utils import *
import os

#load_logging()
config = load_config()

def subset_matrix(mt):
    """
    Finds all the ancestry informative SNPs in our matrix. Creates a subset matrix
    with ancestry SNPs. Export the subset matrix as a VCF.
    """

    mt = mt.key_rows_by('locus', 'alleles')

    # Load ancestry SNPs and parse variant string into locus+alleles
    ht_variants = hl.import_table(f"GrafAnc_SNPs/ancestry_SNPs_GrafAnc_{config['ref_gen']}",
                                  delimiter='\t', impute=True, skip_blank_lines=True)

    ht_variants = ht_variants.annotate(
        parsed=hl.parse_variant(ht_variants.ancestrySNPs)
    )
    ht_variants = ht_variants.annotate(
        locus=ht_variants.parsed.locus,
        alleles=ht_variants.parsed.alleles
    )

    # Match on locus only (mirrors GrafAnc's ExtractAncSnpsFromVcfGz.pl logic — no allele check)
    ht_locus_only = ht_variants.key_by('locus')
    mt_locus_only = mt.key_rows_by('locus')

    ancestry_mt = mt_locus_only.filter_rows(
        hl.is_defined(ht_locus_only[mt_locus_only.locus])
    )
    ancestry_mt = ancestry_mt.key_rows_by('locus', 'alleles')  # restore key for VCF export

    ancestry_counts = ancestry_mt.count_rows()
    logging.info(f"SNPs matched on locus only: {ancestry_counts}")

    # Export VCF
    mt_path = config['mt_afterQC'] if config['preprocessing'] else config['mt_from_vcf']
    ancestry_vcf = os.path.join(os.path.dirname(mt_path), "ancestrySNPs.vcf")
    logging.info(f"Exporting {ancestry_counts} ancestry SNPs to {ancestry_vcf}")
    hl.export_vcf(ancestry_mt, ancestry_vcf)

    summary = []
    summary.append([])
    summary.append(["Ancestry filtering"])
    summary.append(["", "Count"])
    summary.append(["Nº of ancestry SNPs", ancestry_counts])
    csv_writer(summary)

    return ancestry_vcf, ancestry_counts

def call_grafanc(ancestry_vcf, mt_path):
    """
    Calls grafanc and transform results to super populations
    :params: VCF with ancestry informative SNPs and path to the original matrix 
    :return: Ancestry information per sample
    """
    import os, logging
    logging.info(f"cwd before grafanc: {os.getcwd()}")
    logging.info(f"/app/data exists? {os.path.exists('/app/data/AncSnpPopAFs.txt')}")
    logging.info(f"./data exists? {os.path.exists('data/AncSnpPopAFs.txt')}")

    # run GrafAnc
    basename = os.path.basename(mt_path.rstrip('/'))  # create name for results file
    name_only = os.path.splitext(basename)[0]   
    ancestry_results = os.path.join(os.path.dirname(mt_path), f"{name_only}-GrafAnc_results")
    logging.info(f"Running GrafAnc: $ grafanc {ancestry_vcf} {ancestry_results}")
    os.system(f"./grafanc {ancestry_vcf} {ancestry_results} --threads 4")


    try:
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
    
    except FileNotFoundError: 
        logging.error("There are not enough ancestry SNPs to run GrafAnc - aborting ancestry annotation")
        return None

def annotate_ancestry(ancestry_results, mt):
    """
    Annotates GrafAnc results to original matrix.
    :params: ancestry results per sample from GrafAnc and original mt (all the SNPs, ancestry informative and non-informative)
    :return: matrix with ancestry information
    """
    if config['infer_ancestry']:
        ancestry_table = hl.import_table(ancestry_results, delimiter="\t", impute=True)
        ancestry_table = ancestry_table.select('Sample', 'AncGroupNAME')
        ancestry_table = ancestry_table.key_by('Sample')
        mt = mt.annotate_cols(ancestry=ancestry_table[mt.s].AncGroupNAME)
    if config['submit_ancestry']:
        ancestry_table = hl.import_table(ancestry_results, delimiter=",", impute=True)
        ancestry_table = ancestry_table.select('SampleID', 'Population')
        ancestry_table = ancestry_table.key_by('SampleID')
        mt = mt.annotate_cols(ancestry=ancestry_table[mt.s].Population)

    return mt

