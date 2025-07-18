import csv
import yaml 
import logging
import hail as hl 
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import os 

def get_name():
    ## TODO create folder and save every output file there 
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    basename = os.path.basename(config['mt_from_vcf'].rstrip('/'))  # create name for results file
    name_only = f"{os.path.splitext(basename)[0]}_{timestamp}"

    return name_only
 
def csv_creator():
    """
    Create CSV and add proper header
    :param output_csv:
    :return: csv with header and ready to be added new rows
    """
    config['csv_filename'] = f"WorkflowStats_{get_name()}.csv"
    
    with open(config['csv_filename'], mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Summary of Outputs created by EGA standard VCF workflow v1"])
        writer.writerow([])  # Empty line for readability


def csv_writer(new_row):
    
    with open(config['csv_filename'], mode="a", newline="") as file:
        writer = csv.writer(file)
        for row in new_row:
            writer.writerow(row)

def load_config(path="config.yaml"):
    with open(path, "r") as file:
        return yaml.safe_load(file)           

def load_logging():
    name_only = get_name()
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s]: %(message)s",
        handlers=[
            logging.FileHandler(f"{name_only}.log", mode='w'),
            logging.StreamHandler()
        ]
    )

def create_sample_plot(mt): 

    sample_thresholds = config['sample_filters']

    print(sample_thresholds)

    table = mt.cols().select('sample_qc', 'imputed_sex', 'charr')
    df = table.to_pandas()

    # List of metrics to plot ## TODO make it work when not all the fields are available 
    stats = ['sample_qc.dp_stats.mean', 'sample_qc.call_rate',
             'sample_qc.r_het_hom_var', 'sample_qc.n_singleton', 'charr', 'sample_qc.r_ti_tv']

    for metric in stats: 
        plt.figure()

        # Plot boxplot
        df.boxplot(column = metric, by = "imputed_sex")

        # Extract the field name as used in the config
        if metric == "sample_qc.dp_stats.mean": 
            config_key = '.'.join(metric.split('.')[-2:]).upper()
        else:
            
            config_key = metric.split('.')[-1].upper()
        
        threshold = sample_thresholds.get(f"{config_key}_{config['seq_type']}_threshold", None)
        
        if threshold is None: # get the threshold for the fields that doesnt depend on the sequencing type
            threshold = sample_thresholds.get(f"{config_key}_threshold", None)
        
        print(threshold)
        if threshold is not None:
            plt.axhline(y=threshold, color='red', linestyle='--', label='Threshold')
            plt.legend() 

        # Improve layout and labels
        plt.title(f'{metric}')
        plt.ylabel(metric)

        name_only = get_name() # get name depending on the file working with 
        
        plt.savefig(f'{name_only}_{metric}.png')
        plt.clf()



def create_stats(mt, field_path, aggregate):
    # Get the field using bracket notation to handle nested paths
    field_expr = mt[field_path] if '.' not in field_path else eval(f"mt.{field_path}")

    try:
        if aggregate == "rows": # variable stats
            # Use Hail aggregation functions to compute statistics
            stats = mt.aggregate_rows(hl.struct(
                mean=hl.agg.mean(field_expr),
                std=hl.agg.stats(field_expr).stdev,  # or use hl.sqrt(hl.agg.stats(field_expr).variance)
                min=hl.agg.min(field_expr),
                max=hl.agg.max(field_expr),
                count=hl.agg.count_where(hl.is_defined(field_expr)),
                n_missing=hl.agg.count_where(hl.is_missing(field_expr))
            ))
        
        else: # sample stats
            # Use Hail aggregation functions to compute statistics
            stats = mt.aggregate_cols(hl.struct(
                mean=hl.agg.mean(field_expr),
                std=hl.agg.stats(field_expr).stdev,  # or use hl.sqrt(hl.agg.stats(field_expr).variance)
                min=hl.agg.min(field_expr),
                max=hl.agg.max(field_expr),
                count=hl.agg.count_where(hl.is_defined(field_expr)),
                n_missing=hl.agg.count_where(hl.is_missing(field_expr))
            ))
        
        # Convert to regular Python dict
        stats_dict = {
            'mean': stats.mean,
            'std': stats.std,
            'min': stats.min,
            'max': stats.max,
            'n_missing': stats.n_missing
        }
        
        logging.info(f"Statistics for {field_path}: {stats_dict}")

        return stats_dict
        
    except Exception as e:
        logging.error(f"Error calculating statistics for {field_path}: {e}")
        return None
        
    
    




config = load_config()

"""mt = hl.read_matrix_table(config['mt_from_vcf'])
create_stats_samples(mt, mt.sample_qc.dp_stats)"""