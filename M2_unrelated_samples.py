import hail as hl
import yaml
import os 
from utils import *

def quality_control(mt):
    """
    Applies quality control relevant to relatedness estimation. The variants deleted with this QC are NOT deleted from the regular mt. Just deleted for the relatedness estimation. 
    :params: mt without relatedness qc
    :return: mt with relatedness QC
    """

    mt = mt.filter_rows( # these filters are not in the config because they are necessary for a correct relatedness estimation
    (hl.len(mt.alleles) == 2) &  # only biallelic SNPs
    hl.is_snp(mt.alleles[0], mt.alleles[1]) &  # only SNPs
    (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > 0.05) &  # only common SNPs
    (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 < 0.95) & # not present in more than 95% of samples
    (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.95)  # high call rate
    )
    
    n_var = mt.count()

    logging.info(f"Number of variants passing relatedness QC: {n_var[0]} Note: all variants will be included in downstream analysis after relatedness estimation.")

    return mt, n_var

def subset_100000(mt, n_var): 
    """
    Arround 100.000 variants are enough to estimate relatedness with PCA. Subset 100.000 from variants that passed relatedness QC. 
    :params: mt with relatedness qc, number of variants that pass the relatedness qc  
    :return: subset of arround 100.000 variants with relatedness qc 
    """

    logging.info("Creating Subset of arround 100.000 variants passing relatedness QC")
    
    # Calculate the fraction needed to get approximately 100,000 variants
    target_variants = 100000
    if n_var <= target_variants:
        logging.info(f"Dataset already less than {target_variants} variants, keeping all")
        fraction = 1.0
    else:
        fraction = target_variants / n_var
    
    logging.info(f"Sampling fraction: {fraction}    ")
    
    # Checkpoint to materialize any pending operations before sampling
    mt = mt.checkpoint(hl.utils.new_temp_file(extension='mt'))
    # Sample rows 
    mt_filtered = mt.sample_rows(fraction, seed=12345)

    # Verify the final count
    final_n_rows, final_n_cols = mt_filtered.count()
    logging.info(f"Number of variants and samples for PCA: {final_n_rows}, {final_n_cols}")

    
    return mt_filtered

def get_relatedness_output_path():
    basename = os.path.basename(config['mt_from_vcf'].rstrip('/'))
    return f"./{basename}_related_samples_deleted.txt"


def initialize_relatedness_report():
    """
    Create the report file with a clean header.
    """
    out_file = get_relatedness_output_path()
    with open(out_file, "w") as f:
        f.write("sampleId_1\tsampleId_2\tscore\tmethod\n")
    return out_file

def find_high_kin_score(mt):
    """
    Find samples with KING phi > 0.45 (duplicates / twins / near-identical samples).
    Returns a Hail table with one row per pair.
    """
    logging.info("Running KING")
    king_ht = hl.king(mt.GT)

    # Keep only very high kinship pairs
    high_kinship = king_ht.filter(king_ht.phi > 0.45)

    # Convert to row table
    high_kinship_table = high_kinship.entries()

    # Remove self-comparisons
    high_kinship_table = high_kinship_table.filter(high_kinship_table.s_1 != high_kinship_table.s)

    # Remove duplicated symmetric pairs (A-B and B-A), keep only one
    high_kinship_table = high_kinship_table.filter(high_kinship_table.s_1 < high_kinship_table.s)

    return high_kinship_table

def append_king_results_to_file(twins_table):
    """
    Append KING results to the report file with method column.
    """
    out_file = get_relatedness_output_path()

    rows = twins_table.collect()
    with open(out_file, "a") as f:
        for row in rows:
            f.write(f"{row.s_1}\t{row.s}\t{row.phi}\tKING\n")

    logging.info(f"KING-related samples saved in: {out_file}")

def delete_related_indv(mt, mt_filtered):
    """
    Find the samples which relatedness score is > 0.1 (up to about 2nd-degree relatedness). Estimation done by PCA  
    :params: subset of 100.000 variants with relatedness qc
    :return: mt (with all the original variants without relatedness qc) with maximal independent set of related samples
    """
   
    logging.info("hwe_normalized_pca running with subset QC relatedness matrix")
    pca_eigenvalues, pca_scores, _ = hl.hwe_normalized_pca(
        mt_filtered.GT,
        k=5,
        compute_loadings=False
    )

    logging.info("pc_relate running with subset QC relatedness matrix")
    relatedness_ht = hl.pc_relate(
        mt_filtered.GT,
        min_individual_maf=0.01,
        scores_expr=pca_scores[mt_filtered.col_key].scores,
        block_size=1000,
        min_kinship=0.1,
        statistics='all'
    )

    pairs = relatedness_ht.filter(relatedness_ht.kin > 0.1)

    # maximal independent set of samples to remove
    related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, False)

    # collect before filtering the MT
    samples_to_remove = related_samples_to_remove.collect()

    out_file = get_relatedness_output_path()
    with open(out_file, "a") as f:
        for sample in samples_to_remove:
            # sample['node'] is a struct like Struct(s='100c')
            sample_id = sample.node.s
            f.write(f"{sample_id}\t-\t-\tPC_relate\n")

    logging.info(f"PC-relate-related samples saved in: {out_file}")

    mt = mt.filter_cols(hl.is_defined(related_samples_to_remove[mt.col_key]), keep=False)

    return mt

def delete_related_samples(mt):
    """
    Calls all relatedness functions.
    Deletes:
      - KING high-kinship samples (duplicates/twins)
      - PC-relate inferred related samples
    """
    initialize_relatedness_report()

    mt_qual, n_var = quality_control(mt)
    mt_filtered = subset_100000(mt_qual, n_var[0])

    twins_table = find_high_kin_score(mt_filtered)
    append_king_results_to_file(twins_table)

    mt = delete_related_indv(mt, mt_filtered)

    if twins_table.count() != 0:
        # remove one side of each KING pair; currently removes the second sample in each kept pair
        samples_to_remove = twins_table.aggregate(hl.agg.collect_as_set(twins_table.s))
        set_to_remove = hl.literal(samples_to_remove)
        mt = mt.filter_cols(~set_to_remove.contains(mt.s))

    after_rel_trim = mt.count()

    logging.info(
        f"Dataset size after relatedness trimming. "
        f"Variants: {after_rel_trim[0]}, Samples: {after_rel_trim[1]}"
    )

    if config['verbosity']:
        summary = []
        summary.append(["Dataset size after relatedness trimming", after_rel_trim[0], after_rel_trim[1]])
        csv_writer(summary)

    return mt


"""mt = hl.read_matrix_table(config['mt_from_vcf'])
delete_related_samples(mt)"""
