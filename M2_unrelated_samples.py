import hail as hl
import yaml
import os 
from utils import *

config = load_config()
load_logging()

def quality_control(mt):

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

    logging.info("Creating Subset of arround 100.000 variants passing relatedness QC")
    
    # Calculate the fraction needed to get approximately 100,000 variants
    target_variants = 100000
    if n_var <= target_variants:
        logging.info(f"Dataset already less than {target_variants} variants, keeping all")
        fraction = 1.0
    else:
        fraction = target_variants / n_var
    
    logging.info(f"Sampling fraction: {fraction}    ")
    
    # Sample rows (variants)
    mt_for_pca = mt.sample_rows(fraction, seed=12345)  # Added seed for reproducibility
    
    # Write the subset variant information
    mt_for_pca.rows().write(config['table_with_subset'], overwrite=True)
    
    # Read back the rows table and semi-join with original matrix
    subset_rows = hl.read_table(config['table_with_subset'])
    mt_filtered = mt.semi_join_rows(subset_rows)
    
    # Verify the final count
    final_n_rows, final_n_cols = mt_filtered.count()
    logging.info(f"Number of variants and samples for PCA: {final_n_rows}, {final_n_cols}   ")
    
    return mt_filtered


def find_high_kin_score(mt):
    
    king_ht = hl.king(mt.GT)  # KING operates on the genotype data (afterQC_vcfmatrix.GT)

    high_kinship_score = king_ht.filter_entries(king_ht['phi'] > 0.45) # get the pairs of higly releated samples (auto-comparison of samples and twins)

    high_kinship_score_table = high_kinship_score.entries()

    twins_table = high_kinship_score_table.filter(high_kinship_score_table['s_1'] != high_kinship_score_table['s']) # delete the auto-comparion of samples (phi score of 0.5)

    basename = os.path.basename(config['mt_from_vcf'].rstrip('/'))
    twins_table.export(f"./{basename}_related_samples_deleted.txt") 
    
    logging.info(f"High kinship samples deleted by King saved in: ./{os.path.splitext(os.path.basename(config['mt_from_vcf']))[0]}_related_samples_deleted")

    return twins_table

def delete_related_indv(mt, mt_filtered):

    """
    Less principal components (k) -- less memory
    block_size parameter in hl.pc_relate determines how many samples are processed at a time. By reducing this value, you can lower memory usage at the expense of computational time.
    """
    logging.info("hwe_normalized_pca running with subset QC relatedness matrix")
    pca_eigenvalues, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt_filtered.GT, k=5, compute_loadings=False)  # Compute PCA

    logging.info("pc_relate running with subset QC relatedness matrix")

    relatedness_ht = hl.pc_relate(mt_filtered.GT, min_individual_maf=0.01, scores_expr=pca_scores[mt_filtered.col_key].scores,
                                block_size=1000, min_kinship=0.1, statistics='all')  # Compute relatedness

    pairs = relatedness_ht.filter(relatedness_ht['kin'] > 0.1) # Filter pairs based on kinship threshold

    related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, False) # Get the maximal independent set of related samples to remove

    samples_to_remove = related_samples_to_remove.collect() 

    basename = os.path.basename(config['mt_from_vcf'].rstrip('/'))
    with open(f"./{basename}_related_samples_deleted", "a") as file: # Save removed sample ID into a file
        for sample in samples_to_remove:
            file.write(f"{sample['node']}\n")

    mt = mt.filter_cols(hl.is_defined(related_samples_to_remove[mt.col_key]), keep=False) # Remove related samples from the matrix table
    
    return mt 

def delete_related_samples(mt):
    
    mt_qual, n_var = quality_control(mt)

    mt_filtered = subset_100000(mt_qual, n_var[0])
    
    twins_table = find_high_kin_score(mt_filtered)

    mt = delete_related_indv(mt, mt_filtered)

    if twins_table.count() != 0:
        samples_to_remove = twins_table.aggregate(hl.agg.collect_as_set(twins_table.s)) # Remove the high kinship score samples
        set_to_remove = hl.literal(samples_to_remove)
        mt = mt.filter_cols(~set_to_remove.contains(mt['s']))

    after_rel_trim = mt.count() 

    logging.info(f"Dataset size after relatedness trimming. Variants: {after_rel_trim[0]}, Samples: {after_rel_trim[1]}   ")

    return mt 
    

"""mt = hl.read_matrix_table(config['mt_from_vcf'])
delete_related_samples(mt)"""