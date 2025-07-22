import hail as hl
import logging
from utils import *

load_logging()
config = load_config()

def stats_by_sex(mt):
    """
    Creates stats: number of female, number of male and undefined 
    :params: mt 
    :return: stats about sexes in data 
    """
   
    results_sex_agg = mt.aggregate_cols(hl.agg.counter(mt.imputed_sex))

    num_females = results_sex_agg.get("female", 0)
    num_males = results_sex_agg.get("male", 0)
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
    :params: mt, sexes present in the data, ancestries present in the data
    :return: AF fields by sex and ancestry
    """
    logging.info("Re-calculating AFs by sex and ancestry if available")

    mt = mt.annotate_rows(
        gt_stats=hl.agg.call_stats(mt.GT, mt.alleles), # calculate total AF
        gt_stats_by_sex=hl.agg.group_by(mt.imputed_sex, hl.agg.call_stats(mt.GT, mt.alleles)), # calculate by sex AF
    )
    # check ancestry field exists before computing ancestry AF
    if 'ancestry' in mt.col:
        mt = mt.annotate_rows(
            gt_stats_by_ancestry=hl.agg.group_by(mt.ancestry, hl.agg.call_stats(mt.GT, mt.alleles))
        )

    AF_total = mt.gt_stats.AF[1] # get total stats
    AC_total = mt.gt_stats.AC[1]
    AN_total = mt.gt_stats.AN
    AC_hom_total = mt.gt_stats.homozygote_count[1]

    AF_sex = {}
    sexes_avail = list(results_sex_agg.keys()) # gets the sexes in the VCF
    for sex in sexes_avail: # create stats by sex
        if sex == "undefined":
            continue  # Skip samples with unknown sex
        AF_sex[f"AF_{sex}_recalc"] = mt.gt_stats_by_sex[sex].AF[1]
        AF_sex[f"AC_{sex}_recalc"] = mt.gt_stats_by_sex[sex].AC[1]
        AF_sex[f"AC_hom_{sex}_recalc"] = mt.gt_stats_by_sex[sex].homozygote_count[1]
        AF_sex[f"AN_{sex}_recalc"] = mt.gt_stats_by_sex[sex].AN

    logging.warning("Samples where sex couldn't be inferred will be ignored")

    if config['ancestry'] and results_ancestry_agg != "":
        AF_ancestry = {}
        ancestry_avail = list(results_ancestry_agg.keys()) # gets the ancestries in the VCF
        for ancestry in ancestry_avail: # create stats by ancestry
            if ancestry == "Multi-ancestry":
                continue  # Skip samples with multi-ancestry
            AF_ancestry[f"AF_{ancestry}_recalc"] = mt.gt_stats_by_ancestry[ancestry].AF[1]
            AF_ancestry[f"AC_{ancestry}_recalc"] = mt.gt_stats_by_ancestry[ancestry].AC[1]
            AF_ancestry[f"AC_hom_{ancestry}_recalc"] = mt.gt_stats_by_ancestry[ancestry].homozygote_count[1]
            AF_ancestry[f"AN_{ancestry}_recalc"] = mt.gt_stats_by_ancestry[ancestry].AN
    else:
        AF_ancestry={}

    return mt, AF_total, AC_total, AN_total, AC_hom_total, AF_sex, AF_ancestry



def annotate_new_vcf(mt, AF_total, AC_total, AN_total, AC_hom_total, AF_sex, AF_ancestry):
    """
    Annotates mt with AF fields by sex and ancestry
    :params: mt, stast about total AF and AF fields by sex and ancestry
    :return: mt with AF fields by sex and ancestry
    """
    mt_af = mt.annotate_rows(
        info=mt.info.annotate(  # Extend the original info field
            AF_total_recalc=AF_total,
            AC_total_recalc=AC_total,
            AC_hom_total_recalc=AC_hom_total,
            AN_total_recalc= AN_total,
            **AF_sex,
            **AF_ancestry
        )
    )
    return mt_af


def export_new_vcf(mt_af):
    """
    Export mt with ancestry and sex AF fields to VCF 
    :params: mt with AF fields by sex and ancestry
    :return: VCF with AF fields by sex and ancestry
    """
    header_metadata = hl.get_vcf_metadata(config['vcf_for_header'])  # get original header
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
    

