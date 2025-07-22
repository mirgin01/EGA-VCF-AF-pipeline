import csv
import yaml 
import logging
import hail as hl 
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import os 

# Generate timestamp once when module is imported
_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
_name_only = None 

def get_name():
    """
    Checks if name_only has been created. If not, create a global variable for the log and workflow name
    :return: _name_only global variable
    """
    global _name_only
    if _name_only is None:
        basename = os.path.basename(config['mt_from_vcf'].rstrip('/'))  
        _name_only = f"{os.path.splitext(basename)[0]}_{_timestamp}"
        print(_name_only)
    
    return _name_only
 
def csv_creator():
    """
    Create CSV and add proper header
    :return: csv with header and ready to be added new rows
    """
    config['csv_filename'] = f"WorkflowStats_{get_name()}.csv"

    logging.info(f"Writing CSV with step statistics: WorkflowStats_{get_name()}.csv")
    
    with open(config['csv_filename'], mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Summary of Outputs created by EGA standard VCF workflow v1"])
        writer.writerow([])  # Empty line for readability


def csv_writer(new_row):
    """
    add info to csv created in csv_creator
    :params: new row to be added to csv
    :return: csv with new rows
    """
    with open(config['csv_filename'], mode="a", newline="") as file:
        writer = csv.writer(file)
        for row in new_row:
            writer.writerow(row)

def load_config(path="config.yaml"):
    """
    Load config with pipeline parameters
    :params: path of the conf.yaml
    :return: config parameters will be available
    """
    with open(path, "r") as file:
        return yaml.safe_load(file)           

def load_logging():
    """
    Create log file 
    :return: log file with all the outputs
    """
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
    logging.info(f"LOG SAVED IN: {name_only}.log")

def create_sample_plot(mt): 
    """
    Creates plots of the sample quality control metrics, grouped by sex and with the proper thresholds. 
    :params: mt without the sample quality control applied
    :return: boxplot for all the available sample quality control metrics. 
    """
    sample_thresholds = config['sample_filters']

    try: 
        table = mt.cols().select('sample_qc', 'imputed_sex', 'charr')
        df = table.to_pandas()

    except LookupError as e: # if charr is not computed
        table = mt.cols().select('sample_qc', 'imputed_sex')
        df = table.to_pandas()
    
        stats = ['sample_qc.dp_stats.mean', 'sample_qc.call_rate',
                'sample_qc.r_het_hom_var', 'sample_qc.n_singleton', 'charr', 'sample_qc.r_ti_tv']

        for metric in stats: 
            
            try:
                plt.figure()

                # Plot boxplot
                df.boxplot(column = metric, by = "imputed_sex")

                # Extract the field name as used in the config
                if metric == "sample_qc.dp_stats.mean": 
                    config_key = '.'.join(metric.split('.')[-2:]).upper()
                else:
                    
                    config_key = metric.split('.')[-1].upper()
                
                # get threshold from config
                threshold = sample_thresholds.get(f"{config_key}_{config['seq_type']}_threshold", None) # fields where treshold depends on seq type
                
                if threshold is None: # threshold for the fields that doesnt depend on the seq type
                    threshold = sample_thresholds.get(f"{config_key}_threshold", None)
                
                if threshold is not None:
                    if isinstance(threshold, (list, tuple)) and len(threshold) == 2: # when threshold is an interval
                        lower, upper = threshold
                        plt.axhline(y=lower, color='red', linestyle='--', label='Lower Threshold')
                        plt.axhline(y=upper, color='orange', linestyle='--', label='Upper Threshold') 
                    else:    
                        plt.axhline(y=threshold, color='red', linestyle='--', label='Threshold')
                    
                    plt.legend() 

                plt.title(f'{metric}') # set title and label for the plot
                plt.ylabel(metric)

                name_only = get_name() # get name depending on the file working with 
                
                plt.savefig(f'{name_only}_{metric}.png')
                plt.clf()
            
            except KeyError:
                print(f"Skipping metric due to KeyError: {metric}")
                continue
                



def create_stats(mt, field_path, aggregate):
    """
    Creates basic stats about the quality control fields. Informative to understand why variants/samples/genotypes are deleted
    :params: mt without the qc applied for that stat, QC field, is the field from samples or variants
    :return: basic stats about the quality control fields 
    """
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
        
    
def verbosity_counts_variants(mt, filter_name, filter_step, summary, stats):
    """
    Calculates the number of variants deleted per qc step
    :params: mt without QC, filter name, filter command, CSV summary, basic stats to add them to CSV summary
    :return: counts about removed variants and appropiate logs.
    """
    repo_stats = mt.aggregate_rows(hl.struct(
        total = hl.agg.count(),
        passing=hl.agg.count_where(filter_step)
        )
    )   
    logging.info(f"{filter_name} filtering done - Variants removed: {repo_stats.total - repo_stats.passing}  ")
    summary.append([filter_name, repo_stats.total, repo_stats.total - (repo_stats.total-repo_stats.passing), repo_stats.total-repo_stats.passing, stats])

    return summary 



config = load_config()

"""mt = hl.read_matrix_table(config['mt_from_vcf'])
create_stats_samples(mt, mt.sample_qc.dp_stats)"""