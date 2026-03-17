import hail as hl
import logging
from utils import *

config = load_config()

def stats_by_sex(mt):
    """
    Creates stats: number of female, number of male and undefined 
    :params: mt 
    :return: stats about sexes in data 
    """
   
    results_sex_agg = mt.aggregate_cols(hl.agg.counter(mt.imputed_sex))

    num_females = results_sex_agg.get("XX", 0)
    num_males = results_sex_agg.get("XY", 0)
    num_missing = results_sex_agg.get("undefined", 0)

    logging.info(f"Sex Statistics")
    logging.info(f"Number of females: {num_females}")
    logging.info(f"Number of males: {num_males}")
    logging.info(f"Number of samples without sex information: {num_missing}")
    sex_summary = []
    sex_summary.append("")
    sex_summary.append(["Sex Statistics"])
    sex_summary.append(["Number of females", num_females])
    sex_summary.append(["Number of males", num_males])
    sex_summary.append(["Number of samples without sex information:", num_missing])
    csv_writer(sex_summary)

    return mt, results_sex_agg


def stats_by_ancestry(mt):   
    """
    Creates stats: number of samples per population in the data
    :params: mt 
    :return: stats about ancestry groups in the data 
    """
    results_ancestry_agg = mt.aggregate_cols(hl.agg.counter(mt.ancestry))  # agg samples (columns) by ancestry

    if None in results_ancestry_agg:
        logging.error("Some samples without ancestry information")

    logging.info("Ancestry Statistics")
    logging.info(f"Ancestry super-populations: {results_ancestry_agg}") 
    anc_summary = []
    anc_summary.append("")
    anc_summary.append(["Ancestry statistics"])
    for ancestry, count in results_ancestry_agg.items():
        anc_summary.append([ancestry, count])
    csv_writer(anc_summary)
    return mt, results_ancestry_agg

def af_by_sex_ancestry(mt, results_sex_agg, results_ancestry_agg):
    """
    Recalculates AF fields grouped by sex and ancestry 
    :params:
    - mt: Hail MatrixTable containing genotype data
    - results_sex_agg: Dictionary with sex categories present in the dataset
    - results_ancestry_agg: Dictionary with ancestry groups present in the dataset
    
    :return:
    - mt: Updated MatrixTable with calculated statistics
    - AF_total: Overall allele frequency for alternate allele
    - AC_total: Overall allele count for alternate allele  
    - AN_total: Overall allele number (total chromosomes analyzed)
    - nhom_alt_total: Overall count of homozygous alternate individuals
    - AF_sex: Dictionary of sex-stratified frequency statistics
    - AF_ancestry: Dictionary of ancestry-stratified frequency statistics
    """
    logging.info("Re-calculating AFs by sex and ancestry if available")

    has_sex = 'imputed_sex' in mt.col.dtype.keys()

    # PRE-STEP : Ensure haploid calls are correctly written 

    if has_sex:
        mt = mt.annotate_entries(
            GT = hl.case()
            .when(
                # Condition: Male sample AND non-PAR region AND diploid homozygous variant call
                (mt.imputed_sex == "XY") & 
                (~mt.locus.in_autosome_or_par()) & 
                mt.GT.is_diploid() & 
                (mt.GT.is_hom_ref() | mt.GT.is_hom_var()),
                # Convert to haploid: take first allele only
                hl.call(mt.GT[0])
            )
            .default(mt.GT)  # Keep original GT for all other cases
        )


    # STEP 1: Calculate comprehensive genotype statistics
    row_annotations = {
        'gt_stats': hl.agg.call_stats(mt.GT, mt.alleles)
    }
    
    # STEP 2: Conditionally calculate sex-stratified statistics
    if has_sex:
        row_annotations['gt_stats_by_sex'] = hl.agg.group_by(mt.imputed_sex, hl.agg.call_stats(mt.GT, mt.alleles))
    
    # STEP 3: Conditionally calculate ancestry-stratified statistics
    # Only compute ancestry stats if the ancestry field exists in the dataset
    if 'ancestry' in mt.col:
        row_annotations['gt_stats_by_ancestry'] = hl.agg.group_by(mt.ancestry, hl.agg.call_stats(mt.GT, mt.alleles))
        # Only do nested grouping if BOTH exist
        if has_sex:
            row_annotations['gt_stats_by_ancestry_sex'] = hl.agg.group_by(
                hl.struct(ancestry=mt.ancestry, sex=mt.imputed_sex),
                hl.agg.call_stats(mt.GT, mt.alleles)
            )

    mt = mt.annotate_rows(**row_annotations)

    # STEP 4: Extract overall population statistics
    # Index [1] refers to the alternate allele (index [0] would be reference allele)
    AF_total = mt.gt_stats.AF[1] # Alternate allele frequency
    AC_total = mt.gt_stats.AC[1] # Alternate allele count
    AN_total = mt.gt_stats.AN # Total allele number (2 * sample count for diploid)
    nhom_total = mt.gt_stats.homozygote_count[1] # Homozygous alternate count

    # HEMIZYGOTES (Total) - Default to 0 if no sex info
    if has_sex:
        nhemi_total = hl.if_else(~mt.locus.in_autosome_or_par(),
            hl.or_else(mt.gt_stats_by_sex.get("XY").AC[1], 0), 0)
    else:
        nhemi_total = hl.int32(0)

    nhet_total = AC_total - 2*nhom_total - nhemi_total # heterozygous count

    # STEP 5: Extract sex-stratified statistics with hemizygote handling
    AF_sex = {}
    if has_sex:
        sexes_avail = list(results_sex_agg.keys())  # Get available sex categories from input
        logging.warning("Samples where sex couldn't be inferred will be ignored")
        
        for sex in sexes_avail: # create stats by sex
            
            if sex == "undefined":
                continue  # Skip samples with unknown sex
            
            AF_sex[f"AF_{sex}_recalc"] = mt.gt_stats_by_sex[sex].AF[1]
            AF_sex[f"AC_{sex}_recalc"] = mt.gt_stats_by_sex[sex].AC[1]
            AF_sex[f"nhom_{sex}_recalc"] = mt.gt_stats_by_sex[sex].homozygote_count[1]
            AF_sex[f"AN_{sex}_recalc"] = mt.gt_stats_by_sex[sex].AN

            # HEMIZYGOTE HANDLING FOR SEX CHROMOSOMES
            # Add hemizygote allele count if variant is on non-autosomal/non-PAR regions
            # Only XY (male) individuals can be hemizygous for X chromosome variants
            AF_sex[f"nhemi_{sex}_recalc"] = hl.if_else(
                # Check if variant is on sex chromosome outside pseudoautosomal regions
                ~mt.locus.in_autosome_or_par(),
                # If on sex chromosome and this is XY group, use the allele count as hemizygote count
                hl.if_else(
                    sex == "XY", 
                    mt.gt_stats_by_sex[sex].AC[1],  # Males: AC equals hemizygote count
                    hl.int32(0)  # Females: no hemizygotes possible
                ),
                hl.int32(0)  # Autosomal variants: no hemizygotes
            )

            AF_sex[f"nhet_{sex}_recalc"] = AF_sex[f"AC_{sex}_recalc"] - 2*(AF_sex[f"nhom_{sex}_recalc"]) - AF_sex[f"nhemi_{sex}_recalc"] # heterozygous count



    # STEP 6: Extract ancestry-stratified statistics (if available)
    if results_ancestry_agg != "":
        AF_ancestry = {}
        ancestry_avail = list(results_ancestry_agg.keys()) # gets the ancestries in the VCF
        
        for ancestry in ancestry_avail: # create stats by ancestry
            
            if ancestry == "Multi-ancestry":
                continue  # Skip samples with multi-ancestry
            
            AF_ancestry[f"AF_{ancestry}_recalc"] = mt.gt_stats_by_ancestry[ancestry].AF[1]
            AF_ancestry[f"AC_{ancestry}_recalc"] = mt.gt_stats_by_ancestry[ancestry].AC[1]
            AF_ancestry[f"nhom_{ancestry}_recalc"] = mt.gt_stats_by_ancestry[ancestry].homozygote_count[1]
            AF_ancestry[f"AN_{ancestry}_recalc"] = mt.gt_stats_by_ancestry[ancestry].AN

            # HEMIZYGOTE HANDLING FOR ANCESTRY GROUPS
            # For ancestry groups, we need to extract the XY subset within each ancestry
            # This requires a nested grouping: ancestry -> sex -> stats
            # Create the key for this ancestry + XY combination
            if has_sex:
                anc_xy_key = hl.struct(ancestry=ancestry, sex="XY")
                AF_ancestry[f"nhemi_{ancestry}_recalc"] = hl.if_else(
                    ~mt.locus.in_autosome_or_par(),
                    hl.or_else(mt.gt_stats_by_ancestry_sex.get(anc_xy_key).AC[1], 0), 0)
            else:
                AF_ancestry[f"nhemi_{ancestry}_recalc"] = hl.int32(0)

            AF_ancestry[f"nhet_{ancestry}_recalc"] = AF_ancestry[f"AC_{ancestry}_recalc"] - 2*AF_ancestry[f"nhom_{ancestry}_recalc"] - AF_ancestry[f"nhemi_{ancestry}_recalc"]
    else:
        AF_ancestry={}

        # STEP 7: Extract ancestry-by-sex stratified statistics (if both are available)
    # Produces fields like AF_EUR_XY_recalc, AF_EUR_XX_recalc, etc.
    AF_ancestry_sex = {}

    if results_ancestry_agg != "" and has_sex:
        
        for ancestry in ancestry_avail:

            if ancestry == "Multi-ancestry":
                continue

            for sex in sexes_avail:

                if sex == "undefined":
                    continue

                key = hl.struct(ancestry=ancestry, sex=sex)

                AF_ancestry_sex[f"AF_{ancestry}_{sex}_recalc"]   = hl.or_else(mt.gt_stats_by_ancestry_sex.get(key).AF[1],   hl.float64(0))
                AF_ancestry_sex[f"AC_{ancestry}_{sex}_recalc"]   = hl.or_else(mt.gt_stats_by_ancestry_sex.get(key).AC[1],   hl.int32(0))
                AF_ancestry_sex[f"AN_{ancestry}_{sex}_recalc"]   = hl.or_else(mt.gt_stats_by_ancestry_sex.get(key).AN,      hl.int32(0))
                AF_ancestry_sex[f"nhom_{ancestry}_{sex}_recalc"] = hl.or_else(mt.gt_stats_by_ancestry_sex.get(key).homozygote_count[1], hl.int32(0))

                # HEMIZYGOTE HANDLING FOR ANCESTRY x SEX
                # Only XY individuals outside autosomal/PAR regions can be hemizygous
                AF_ancestry_sex[f"nhemi_{ancestry}_{sex}_recalc"] = hl.if_else(
                    ~mt.locus.in_autosome_or_par(),
                    hl.if_else(
                        sex == "XY",
                        hl.or_else(mt.gt_stats_by_ancestry_sex.get(key).AC[1], hl.int32(0)),
                        hl.int32(0)  # XX individuals: no hemizygotes
                    ),
                    hl.int32(0)  # Autosomal variants: no hemizygotes
                )

                AF_ancestry_sex[f"nhet_{ancestry}_{sex}_recalc"] = (
                    AF_ancestry_sex[f"AC_{ancestry}_{sex}_recalc"]
                    - 2 * AF_ancestry_sex[f"nhom_{ancestry}_{sex}_recalc"]
                    - AF_ancestry_sex[f"nhemi_{ancestry}_{sex}_recalc"]
                )

    return mt, AF_total, AC_total, AN_total, nhom_total, nhet_total, nhemi_total, AF_sex, AF_ancestry, AF_ancestry_sex



def annotate_new_vcf(mt, AF_total, AC_total, AN_total, homozygote_count_total, heterozygous_count, hemizygotes_total, AF_sex, AF_ancestry, AF_ancestry_sex):
    """
    Annotates mt with AF fields by sex and ancestry
    :params: mt, stast about total AF and AF fields by sex and ancestry
    :return: mt with AF fields by sex and ancestry
    """
    mt_af = mt.annotate_rows(
        info=mt.info.annotate(  # Extend the original info field
            AF_total_recalc=AF_total,
            AC_total_recalc=AC_total,
            nhom_alt_total_recalc=homozygote_count_total,
            nhet_total_recalc = heterozygous_count,
            nhemi_total_recalc=hemizygotes_total,
            AN_total_recalc= AN_total,
            **AF_sex,
            **AF_ancestry, 
            **AF_ancestry_sex
        )
    )
    return mt_af


def export_new_vcf(mt_af):
    """
    Export mt with ancestry and sex AF fields to VCF 
    :params: mt with AF fields by sex and ancestry
    :return: VCF with AF fields by sex and ancestry
    """
    # Get all VCF file paths in the directory
    vcfs = sorted([
        os.path.join(config['vcf_dir'], f)
        for f in os.listdir(config['vcf_dir'])
        if f.endswith('.vcf') or f.endswith('.vcf.bgz') or f.endswith('.vcf.gz')
    ])
    
    header_metadata = hl.get_vcf_metadata(vcfs[0])  # get original header
    
    if config['summary_VCF'] == True:
        logging.info("Exporting final summary VCF with recalculated AF fields and no sample information")
        # Remove all sample columns
        variants_table = mt_af.rows()
        hl.export_vcf(
        variants_table,
        config['final_vcf_AF'],
        metadata=header_metadata,
        tabix=True
    )

    else:
        logging.info("Exporting final VCF with recalculated AF fields and sample information")
        hl.export_vcf(
        mt_af,
        config['final_vcf_AF'],
        metadata=header_metadata,
        tabix=True
        )
    

