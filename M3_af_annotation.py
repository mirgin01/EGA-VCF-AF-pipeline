import hail as hl
import yaml
import logging
from csv_utils import *

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


def stats_by_sex(mt):

    imputed_sex = hl.impute_sex(mt.GT, aaf_threshold=0.05, female_threshold=0.5, male_threshold=0.75)  # Imputed sex with suggested thresholds
    mt = mt.annotate_cols(pheno=imputed_sex[mt.s])  # Annotate matrix table with imputed sex

    results_sex_agg = mt.aggregate_cols(hl.agg.counter(mt.pheno.is_female))

    num_females = results_sex_agg.get(True, 0)
    num_males = results_sex_agg.get(False, 0)
    num_missing = results_sex_agg.get(None, 0)

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
    csv_writer(sex_summary, "QCsummary.csv")

    return mt, results_sex_agg


def stats_by_ancestry(mt):

    results_ancestry_agg = mt.aggregate_cols(hl.agg.counter(mt.ancestry))  # agg samples (columns) by ancestry

    if None in results_ancestry_agg:
        logging.error("Some samples without ancestry information")

    logging.info("Ancestry Statistics")
    logging.info(f"Ancestry super-populations: {results_ancestry_agg}")  ## TODO add to csv
    anc_summary = []
    anc_summary.append("")
    anc_summary.append(["Ancestry statistics"])
    for ancestry, count in results_ancestry_agg.items():
        anc_summary.append([ancestry, count])
    csv_writer(anc_summary, "QCsummary.csv")
    return mt, results_ancestry_agg

def af_by_sex_ancestry(mt, results_sex_agg, results_ancestry_agg):
    logging.info("Re-calculating AFs by sex and ancestry if available")

    mt = mt.annotate_rows(
        gt_stats=hl.agg.call_stats(mt.GT, mt.alleles), # calculate total AF
        gt_stats_by_sex=hl.agg.group_by(mt.pheno.is_female, hl.agg.call_stats(mt.GT, mt.alleles)), # calculate by sex AF
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
        if sex is None:
            continue  # Skip samples with unknown sex
        AF_sex[f"AF_{sex}"] = mt.gt_stats_by_sex[sex].AF[1]
        AF_sex[f"AC_{sex}"] = mt.gt_stats_by_sex[sex].AC[1]
        AF_sex[f"AC_hom_{sex}"] = mt.gt_stats_by_sex[sex].homozygote_count[1]
        AF_sex[f"AN_{sex}"] = mt.gt_stats_by_sex[sex].AN

    logging.warning("Samples where sex couldn't be inferred will be ignored")

    if config['ancestry']:
        AF_ancestry = {}
        ancestry_avail = list(results_ancestry_agg.keys()) # gets the ancestries in the VCF
        for ancestry in ancestry_avail: # create stats by ancestry
            AF_ancestry[f"AF_{ancestry}"] = mt.gt_stats_by_ancestry[ancestry].AF[1]
            AF_ancestry[f"AC_{ancestry}"] = mt.gt_stats_by_ancestry[ancestry].AC[1]
            AF_ancestry[f"AC_hom_{ancestry}"] = mt.gt_stats_by_ancestry[ancestry].homozygote_count[1]
            AF_ancestry[f"AN_{ancestry}"] = mt.gt_stats_by_ancestry[ancestry].AN
    else:
        AF_ancestry={}

    return mt, AF_total, AC_total, AN_total, AC_hom_total, AF_sex, AF_ancestry



def annotate_new_vcf(mt, AF_total, AC_total, AN_total, AC_hom_total, AF_sex, AF_ancestry):
    mt_af = mt.annotate_rows(
        info=mt.info.annotate(  # Extend the original info field
            AF_total=AF_total,
            AC_total=AC_total,
            AC_hom_total=AC_hom_total,
            AN_total=AN_total,
            **AF_sex,
            **AF_ancestry
        )
    )
    return mt_af


def export_new_vcf(mt_af):
    header_metadata = hl.get_vcf_metadata(config['vcf_for_header'])  # get original header
    print("exporting")
    hl.export_vcf(
        mt_af,
        config['final_vcf_QC_AF'],
        metadata=header_metadata,
        tabix=True
    )

