from M1_preprocessing import *
from M2_ancestry import *
from M3_af_annotation import *
import csv
import yaml
from csv_utils import *


def load_config(path="config.yaml"):
    with open(path, "r") as file:
        return yaml.safe_load(file)


def main():
    config = load_config()

    ## PREPROCESSING

    if config['preprocessing']:

        if config['convert_vcfs']:
            mt = convert_and_merge_vcfs(
                config['vcf_dir'],
                config['mt_from_vcf'],
                config['ref_gen']
            )

        mt = hl.read_matrix_table(config['mt_from_vcf'])

        csv_creator("QCsummary.csv")

        general_counts = []
        general_counts.append([])
        general_counts.append(["", "Variants", "Samples"])
        original_size = mt.count()


        if config['split_multiallelic']:
            mt = split_multiallelic(mt)
        if config['variant_filtering']:
            mt = gnomad_variant_filtering(mt)
        if config['sample_filtering']:
            mt = gnomad_sample_filtering(mt, config['seq_type'])

        general_counts.append(["Original dataset size", original_size[0], original_size[1]])
        final_size = mt.count()
        general_counts.append(["Final dataset size", final_size[0], final_size[1]])
        csv_writer(general_counts, "QCsummary.csv")

        logging.info(f"Final number of variants: {final_size[0]}")
        logging.info(f"Final number of samples: {final_size[1]}")

        mt.write(config['mt_afterQC'], overwrite=True)  # write matrix with QC
        logging.info(f"MT with QC written in: {config['mt_afterQC']}")

    else:
        mt = hl.read_matrix_table(config['mt_from_vcf'])

        csv_creator("QCsummary.csv")

        general_counts = []
        general_counts.append([])
        general_counts.append(["", "Variants", "Samples"])
        original_size = mt.count()
        general_counts.append(["Original dataset size", original_size[0], original_size[1]])
        csv_writer(general_counts, "QCsummary.csv")

    ## ANCESTRY

    if config['ancestry']:
        ancestry_vcf = subset_matrix(config['mt_afterQC'])
        ancestry_results = call_grafanc(ancestry_vcf, config['mt_afterQC'])
        mt = annotate_ancestry(ancestry_results, config['mt_afterQC'])
    else:
        mt = hl.read_matrix_table(config['mt_afterQC'])  # read matrix for AF recalc if ancestry is not used

    ## AF RECALC

    if config['af_annotation']:
        mt, results_sex_agg = stats_by_sex(mt)
        if config['ancestry']:
            mt, results_ancestry_agg = stats_by_ancestry(mt)
        else:
            results_ancestry_agg = ""

        mt, AF_total, AC_total, AN_total, AC_hom_total, AF_sex, AF_ancestry = af_by_sex_ancestry(mt,
                                                                                             results_sex_agg,
                                                                                             results_ancestry_agg)

        mt = annotate_new_vcf(mt, AF_total, AC_total, AN_total, AC_hom_total, AF_sex, AF_ancestry)
        export_new_vcf(mt)




if __name__ == "__main__":
    main()


