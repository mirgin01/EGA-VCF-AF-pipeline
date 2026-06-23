import csv
import yaml 
import logging
import hail as hl 
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import os 
import pandas as pd 
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy import stats
import matplotlib.patches as mpatches

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
            logging.FileHandler(f"{name_only}.log", mode='w', encoding="utf-8", delay=True),
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

    if config["infer_sex"]: 
        # Impute and export inferred sex for the samples 
        mt = impute_sex(mt)
        
    # 1. Check what columns actually exist in the MatrixTable
    available_cols = list(mt.col.dtype.keys())
    has_sex = 'imputed_sex' in available_cols
    has_charr = 'charr' in available_cols

    # 2. Build the selection list dynamically
    select_fields = ['sample_qc']
    if has_sex:
        select_fields.append('imputed_sex')
    if has_charr:
        select_fields.append('charr')

    # Convert to pandas
    table = mt.cols().select(*select_fields)
    df = table.to_pandas()

    stats = ['sample_qc.dp_stats.mean', 'sample_qc.call_rate',
            'sample_qc.r_het_hom_var', 'sample_qc.n_singleton', 'sample_qc.r_ti_tv']
    
    if has_charr:
        stats.append('charr')

    for metric in stats: 
        
        try:
            plt.figure()

            # 3. Conditional plotting based on sex availability
            if has_sex:
                df.boxplot(column=metric, by="imputed_sex")
            else:
                df.boxplot(column=metric) # Single group plot

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

# ── helpers ──────────────────────────────────────────────────────────────────

QC_METRICS = [
    'sample_qc.dp_stats.mean',
    'sample_qc.call_rate',
    'sample_qc.r_het_hom_var',
    'sample_qc.n_singleton',
    'sample_qc.r_ti_tv',
]

METRIC_LABELS = {
    'sample_qc.dp_stats.mean': 'Mean DP',
    'sample_qc.call_rate':     'Call Rate',
    'sample_qc.r_het_hom_var': 'Het/Hom',
    'sample_qc.n_singleton':   'N Singleton',
    'sample_qc.r_ti_tv':       'Ti/Tv',
    'charr':                   'CHARR',
}


def _prepare_df(mt, include_charr: bool = True):
    """Extract QC metrics + optional metadata from a Hail MatrixTable."""
    available = list(mt.col.dtype.keys())
    has_charr = 'charr' in available and include_charr

    select_fields = ['sample_qc']
    if has_charr: select_fields.append('charr')

    df = mt.cols().select(*select_fields).to_pandas()

    metrics = list(QC_METRICS)
    if has_charr:
        metrics.append('charr')

    # keep only metrics that are actually present
    metrics = [m for m in metrics if m in df.columns]
    print(metrics)

    meta_cols = ['s'] if 's' in df.columns else []

    return df, metrics, meta_cols

def _clean_metric_matrix(df: pd.DataFrame, metrics: list) -> tuple:
    """Return (X_clean, valid_idx) with no NaNs, cast to float64."""
    sub = df[metrics].copy().apply(pd.to_numeric, errors='coerce')
    valid_mask = sub.notna().all(axis=1)
    return sub[valid_mask].values.astype(float), sub[valid_mask].index


def plot_pca_qc(mt, z_thresh: float = 3.0, name_prefix: str = 'sample_qc',
                collected: dict = None, sample_ids: list = None):
    """
    PCA across all QC metrics to detect global outliers.

    Parameters
    ----------
    mt          : Hail MatrixTable (after sample_qc, optional impute_sex / charr)
    z_thresh    : PC-score Z-score threshold to flag outliers (default 3.0)
    name_prefix : prefix for saved .png files
    """
    df, metrics, meta_cols = _prepare_df(mt)
    X, valid_idx = _clean_metric_matrix(df, metrics)

    # standardise + PCA 
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    n_components = min(len(metrics), X_scaled.shape[0], 5)
    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(X_scaled)

    pc_df = pd.DataFrame(
        pcs,
        columns=[f'PC{i+1}' for i in range(n_components)],
        index=valid_idx
    )

    # flag outliers by Z-score on PC1 & PC2
    for pc in ['PC1', 'PC2']:
        pc_df[f'z_{pc}'] = np.abs(stats.zscore(pc_df[pc]))
    pc_df['outlier'] = (pc_df['z_PC1'] > z_thresh) | (pc_df['z_PC2'] > z_thresh)
        
    # start - print 5 samples furthest from the PCA centre with their values for the QC metrics, and cohort-wide stats for those metrics

    if collected is not None and sample_ids is not None:
            # Euclidean distance from origin across PC1, PC2, PC3
            pc_df['distance'] = np.sqrt(
                pc_df['PC1']**2 + pc_df['PC2']**2 + pc_df['PC3']**2
            )
            top5_idx = pc_df['distance'].nlargest(5).index
            top5_names = [sample_ids[i] if i < len(sample_ids) else 'NA' for i in top5_idx]
            logging.info(f"Top 5 sample IDs: {', '.join(top5_names)}")

            # cohort-wide stats for reference
            logging.info("=== Top 5 samples furthest from PCA centre ===")
            logging.info(f"{'Metric':<20} {'Cohort mean':>12} {'Cohort min':>12} {'Cohort max':>12} {'Sample values (top5)':>30}")
            logging.info("-" * 90)

            metric_labels = {
                'call_rate':     'Call Rate',
                'r_ti_tv':       'Ti/Tv ratio',
                'r_het_hom_var': 'Het/Hom ratio',
                'n_singleton':   'N Singletons',
                'dp_stats_mean': 'Mean DP',
                'charr':         'CHARR',
            }

            for key, label in metric_labels.items():
                vals = collected.get(key)
                if vals is None:
                    continue
                valid = vals[~np.isnan(vals)]
                if len(valid) == 0:
                    continue

                cohort_mean = float(np.mean(valid))
                cohort_min  = float(np.min(valid))
                cohort_max  = float(np.max(valid))

                # get values for the top 5 samples by their index in sample_ids
                top5_vals = []
                for idx in top5_idx:
                    # idx is a pandas integer index into pc_df
                    # we need to map it back to position in sample_ids
                    if idx < len(sample_ids):
                        pos = idx
                        top5_vals.append(f"{vals[pos]:.4f}" if not np.isnan(vals[pos]) else "NA")
                    else:
                        top5_vals.append("NA")

                logging.info(
                    f"{label:<20} {cohort_mean:>12.4f} {cohort_min:>12.4f} "
                    f"{cohort_max:>12.4f}   {', '.join(top5_vals)}"
                )

            logging.info("=" * 90)
    # fin - print 5 samples furthest from the PCA centre with their values for the QC metrics, and cohort-wide stats for those metrics
    
    # attach sample IDs
    if 's' in df.columns:
        pc_df['sample_id'] = df.loc[valid_idx, 's'].values

    n_out = pc_df['outlier'].sum()
    print(f"[PCA] {n_out} / {len(pc_df)} samples flagged as outliers "
          f"(|Z| > {z_thresh} on PC1 or PC2)")

    # figure 1 : PC1 vs PC2 
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('PCA on QC Metrics — Global Outlier Detection', fontsize=14, fontweight='bold')

    ax = axes[0]
    ax.scatter(pc_df['PC1'], pc_df['PC2'],
               c=pc_df['outlier'].map({True: 'red', False: '#4C9BE8'}),
               alpha=0.7, edgecolors='none', s=40)
    handles = [
        mpatches.Patch(color='#4C9BE8', label='Pass'),
        mpatches.Patch(color='red',     label='Outlier'),
    ]
    ax.legend(handles=handles, fontsize=8)

    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)', fontsize=10)
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)', fontsize=10)
    ax.set_title('PC1 vs PC2', fontsize=11)
    ax.axhline(0, color='grey', linewidth=0.5, linestyle='--')
    ax.axvline(0, color='grey', linewidth=0.5, linestyle='--')

    # label outlier sample IDs if available
    if 'sample_id' in pc_df.columns:
        for _, row in pc_df[pc_df['outlier']].iterrows():
            ax.annotate(str(row['sample_id']),
                        xy=(row['PC1'], row['PC2']),
                        fontsize=6, color='darkred',
                        xytext=(4, 4), textcoords='offset points')

    # right panel : scree + cumulative variance
    ax2 = axes[1]
    var_exp  = pca.explained_variance_ratio_ * 100
    cum_var  = np.cumsum(var_exp)
    pc_names = [f'PC{i+1}' for i in range(n_components)]

    bars = ax2.bar(pc_names, var_exp, color='#4C9BE8', alpha=0.8, label='Individual')
    ax2b = ax2.twinx()
    ax2b.plot(pc_names, cum_var, 'o--', color='#E87B4C', linewidth=2, label='Cumulative')
    ax2b.set_ylabel('Cumulative Variance (%)', color='#E87B4C', fontsize=10)
    ax2b.tick_params(axis='y', labelcolor='#E87B4C')
    ax2b.set_ylim(0, 110)

    ax2.set_xlabel('Principal Component', fontsize=10)
    ax2.set_ylabel('Explained Variance (%)', fontsize=10)
    ax2.set_title('Scree Plot', fontsize=11)

    plt.tight_layout()
    out_path = f'{name_prefix}_pca_overview.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.clf()
    print(f"  → saved {out_path}")

    # figure 2 : loadings heatmap 
    loadings = pd.DataFrame(
        pca.components_[:, :].T,
        index=[METRIC_LABELS.get(m, m) for m in metrics],
        columns=pc_names
    )

    fig2, ax3 = plt.subplots(figsize=(max(6, n_components * 1.2), len(metrics) * 0.7 + 1))
    sns.heatmap(loadings, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
                linewidths=0.5, ax=ax3, vmin=-1, vmax=1,
                cbar_kws={'label': 'Loading'})
    ax3.set_title('PCA Loadings — contribution of each metric to each PC',
                  fontsize=11, fontweight='bold')
    plt.tight_layout()
    lpath = f'{name_prefix}_pca_loadings.png'
    plt.savefig(lpath, dpi=150, bbox_inches='tight')
    plt.clf()
    print(f"  → saved {lpath}")

    return pc_df   # caller can inspect / further filter


def plot_zscore_heatmap(mt, z_thresh: float = 3.0,
                        max_samples_display: int = 200,
                        name_prefix: str = 'sample_qc'):
    """
    Per-sample Z-score heatmap across all QC metrics.

    Parameters
    ----------
    mt                  : Hail MatrixTable
    z_thresh            : |Z| threshold to flag a metric as an outlier
    max_samples_display : cap for legibility — only the worst samples are shown
                          when the cohort is large
    name_prefix         : prefix for saved .png files
    """
    df, metrics, meta_cols = _prepare_df(mt)
    X, valid_idx = _clean_metric_matrix(df, metrics)

    z_matrix = pd.DataFrame(
        stats.zscore(X, axis=0),
        index=valid_idx,
        columns=[METRIC_LABELS.get(m, m) for m in metrics]
    )

    # flag samples: number of metrics with |Z| > threshold
    z_matrix['n_outlier_metrics'] = (z_matrix.abs() > z_thresh).sum(axis=1)
    z_matrix_sorted = z_matrix.sort_values('n_outlier_metrics', ascending=False)

    # attach sample IDs for the y-axis
    if 's' in df.columns:
        sample_ids = df.loc[valid_idx, 's']
        z_matrix_sorted.index = sample_ids[z_matrix_sorted.index].values

    n_flagged = (z_matrix_sorted['n_outlier_metrics'] > 0).sum()
    print(f"[Z-score] {n_flagged} / {len(z_matrix_sorted)} samples have at least "
          f"one metric with |Z| > {z_thresh}")

    metric_cols = [METRIC_LABELS.get(m, m) for m in metrics]

    # figure 1 : full cohort overview (clipped display) 
    display_df = z_matrix_sorted[metric_cols].head(max_samples_display)

    row_h   = max(0.18, min(0.5, 12 / max(len(display_df), 1)))
    fig_h   = max(5, len(display_df) * row_h + 2)
    fig_w   = max(8, len(metrics) * 1.2 + 2)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        display_df,
        cmap='RdBu_r', center=0,
        vmin=-4, vmax=4,
        linewidths=0.3 if len(display_df) < 80 else 0,
        ax=ax,
        yticklabels=(len(display_df) <= 80),
        cbar_kws={'label': 'Z-score', 'shrink': 0.6}
    )

    # draw threshold contour lines
    ax.axvline(x=0, color='black', linewidth=0.5)

    title_suffix = (f' (top {max_samples_display} by #flagged metrics)'
                    if len(z_matrix_sorted) > max_samples_display else '')
    ax.set_title(
        f'Z-score Heatmap — QC Metrics per Sample{title_suffix}\n'
        f'|Z| > {z_thresh} = outlier   |   sorted by number of flagged metrics',
        fontsize=11, fontweight='bold'
    )
    ax.set_xlabel('QC Metric', fontsize=10)
    ax.set_ylabel('Sample', fontsize=10)
    plt.tight_layout()

    hpath = f'{name_prefix}_zscore_heatmap.png'
    plt.savefig(hpath, dpi=150, bbox_inches='tight')
    plt.clf()
    print(f"  → saved {hpath}")

    # figure 2 : bar chart — how many samples fail N metrics 
    counts = z_matrix_sorted['n_outlier_metrics'].value_counts().sort_index()

    fig2, ax2 = plt.subplots(figsize=(8, 4))
    bars = ax2.bar(counts.index.astype(str), counts.values,
                   color=['#4C9BE8' if i == 0 else '#E87B4C' for i in counts.index],
                   edgecolor='white', linewidth=0.8)
    for bar, val in zip(bars, counts.values):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                 str(val), ha='center', va='bottom', fontsize=9)

    ax2.set_xlabel('Number of metrics with |Z| > threshold', fontsize=10)
    ax2.set_ylabel('Number of samples',                      fontsize=10)
    ax2.set_title(f'Samples by Number of Flagged QC Metrics (|Z| > {z_thresh})',
                  fontsize=11, fontweight='bold')
    ax2.set_xticks(range(len(counts)))
    ax2.set_xticklabels(counts.index.astype(str))
    plt.tight_layout()

    bpath = f'{name_prefix}_zscore_flagcount.png'
    plt.savefig(bpath, dpi=150, bbox_inches='tight')
    plt.clf()
    print(f"  → saved {bpath}")

    # figure 3 : per-metric Z-score distributions 
    n_met = len(metrics)
    ncols = min(3, n_met)
    nrows = int(np.ceil(n_met / ncols))

    fig3, axes3 = plt.subplots(nrows, ncols,
                                figsize=(ncols * 4, nrows * 3),
                                squeeze=False)
    for i, metric in enumerate(metric_cols):
        row, col = divmod(i, ncols)
        ax3 = axes3[row][col]
        z_vals = z_matrix_sorted[metric]

        ax3.hist(z_vals, bins=40, color='#4C9BE8', edgecolor='white', alpha=0.8)
        ax3.axvline( z_thresh, color='red',  linestyle='--', linewidth=1.2, label=f'+{z_thresh}σ')
        ax3.axvline(-z_thresh, color='red',  linestyle='--', linewidth=1.2, label=f'-{z_thresh}σ')
        ax3.axvline( 0,        color='grey', linestyle='-',  linewidth=0.8)

        n_out = (z_vals.abs() > z_thresh).sum()
        ax3.set_title(f'{metric}\n(n outliers = {n_out})', fontsize=9)
        ax3.set_xlabel('Z-score', fontsize=8)
        ax3.set_ylabel('Count',   fontsize=8)
        ax3.legend(fontsize=7)

    # hide unused subplots
    for j in range(i + 1, nrows * ncols):
        row, col = divmod(j, ncols)
        axes3[row][col].set_visible(False)

    fig3.suptitle('Per-Metric Z-score Distributions', fontsize=12, fontweight='bold')
    plt.tight_layout()
    dpath = f'{name_prefix}_zscore_distributions.png'
    plt.savefig(dpath, dpi=150, bbox_inches='tight')
    plt.clf()
    print(f"  → saved {dpath}")

    return z_matrix_sorted   # caller can inspect


def create_outlier_plots(mt, z_thresh: float = 3.0, name_prefix: str = 'sample_qc', collected: dict = None, sample_ids: list = None):
    """
    Run both PCA and Z-score analyses and save all plots.

    Returns
    -------
    pc_df       : DataFrame with PC scores + outlier flag
    zscore_df   : DataFrame with Z-scores + n_outlier_metrics column
    """
    
    print("=== PCA on QC metrics ===")
    pc_df = plot_pca_qc(mt, z_thresh=z_thresh, name_prefix=name_prefix, collected=collected, sample_ids=sample_ids)

    print("\n=== Z-score heatmap ===")
    zscore_df = plot_zscore_heatmap(mt, z_thresh=z_thresh, name_prefix=name_prefix)

    # ── combined outlier summary ──────────────────────────────────────────
    pca_outliers    = set(pc_df[pc_df['outlier']].index)
    zscore_outliers = set(zscore_df[zscore_df['n_outlier_metrics'] > 0].index)
    both            = pca_outliers & zscore_outliers

    print(f"\n{'='*50}")
    print(f"  PCA outliers          : {len(pca_outliers)}")
    print(f"  Z-score outliers      : {len(zscore_outliers)}")
    print(f"  Flagged by BOTH       : {len(both)}  ← high-confidence")
    print(f"{'='*50}")

    return pc_df, zscore_df

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

"""def verbosity_counts_variants(mt, filter_name, filter_step, summary, stats):
    # filter_step is the REMOVAL condition, so survivors are where it's False or missing
    keep_expr = hl.is_missing(filter_step) | ~filter_step
    
    repo_stats = mt.aggregate_rows(hl.struct(
        total=hl.agg.count(),
        passing=hl.agg.count_where(hl.or_else(keep_expr, True))
    ))
    logging.info(f"{filter_name} filtering done - Variants removed: {repo_stats.total - repo_stats.passing}")
    summary.append([filter_name, repo_stats.total, repo_stats.passing, repo_stats.total - repo_stats.passing, stats])
    return summary"""


def verbosity_counts_variants(mt, filter_name, keep_condition, summary, stats):
    repo_stats = mt.aggregate_rows(hl.struct(
        total=hl.agg.count(),
        passing=hl.agg.count_where(hl.or_else(keep_condition, False))
    ))
    removed = repo_stats.total - repo_stats.passing
    logging.info(f"{filter_name} filtering done - Variants removed: {removed}")
    summary.append([filter_name, repo_stats.total, repo_stats.passing, removed, stats])
    return summary

def run_charr(mt, reference):
    if "AF" in mt.info:
        gnomad_sites = hl.read_table(reference)
        charr_result = hl.compute_charr(
            mt,
            ref_AF=(1 - gnomad_sites[mt.row_key].freq[0].AF)
        )
        #annotate charr results
        mt = mt.annotate_cols(charr=charr_result[mt.col_key].charr)
    
    return mt 

def impute_sex(mt):
    """
    Imputes sex for all the samples
    :params: Imputes sex with original data
    :return: Imputed sex by sample
    """
    
    imputed_sex = hl.impute_sex(mt.GT, aaf_threshold=0.05, female_threshold=0.5, male_threshold=0.75)  # Imputed sex with suggested thresholds
    sex_expr = hl.if_else(hl.is_defined(imputed_sex.is_female), hl.if_else(imputed_sex.is_female, # rename imputed sex
                                                                        "XX",
                                                                        "XY"),
                                                                        "undefined")
    sex_ht = imputed_sex.annotate(imputed_sex=sex_expr)

    # annotate input (all chroms) mt with imputed sex
    sex_colnames = ["f_stat", "is_female", "imputed_sex"]
    mt = mt.annotate_cols(**sex_ht.select(*sex_colnames)[mt.col_key])

    # Save the sex imputation results
    basename = os.path.basename(config['mt_from_vcf'].rstrip('/'))
    sex_ht.select("imputed_sex", "f_stat", "is_female").export(f"{basename}_imputed_sex_results.tsv")
    logging.info(f"Sex inference saved in: {basename}_imputed_sex_results.tsv")

    return mt

def print_config(d, prefix=""):
    """
    Add config options to log
    """
    if isinstance(d, dict):
        for k, v in d.items():
            if isinstance(v, dict):
                logging.info(f" {prefix}{k}:")
                print_config(v, prefix + "  ")
            else:
                logging.info(f" {prefix}{k} : {v}")
    else:
        logging.info(f" {prefix}{d}")

config = load_config()

def rename_chr(mt, ref_gen): 
    """
    Rename chromosomes 23 and 24 to X and Y respectively if they exist.
    
    :parms:
    mt : MatrixTable
        Input Hail MatrixTable
    reference_genome : str
        Reference genome name (default: 'GRCh38')
    
    :return: mt with renamed chromosomes
    """
    
    # Get unique contigs in the MatrixTable
    contigs = mt.aggregate_rows(hl.agg.collect_as_set(mt.locus.contig))
    
    # Check if '23' or '24' exist in the contigs
    has_23 = '23' in contigs
    has_24 = '24' in contigs
    
    # Create mapping dictionary
    contig_map = {}
    if has_23:
        contig_map['23'] = 'X'
        logging.info("Renaming chromosome 23 to X")
    if has_24:
        contig_map['24'] = 'Y'
        logging.info("Renaming chromosome 24 to Y")
    
    # Create new contig name using conditional expression
    new_contig = mt.locus.contig
    for old_name, new_name in contig_map.items():
        new_contig = hl.if_else(
            mt.locus.contig == old_name,
            new_name,
            new_contig
        )
    
    # Create new locus with renamed contigs
    new_locus = hl.locus(
        new_contig,
        mt.locus.position,
        reference_genome=ref_gen
    )
    
    # Check if alleles is part of the key
    original_key = list(mt.row_key)
    has_alleles_key = 'alleles' in original_key
    
    # Use key_rows_by to modify the locus key
    if has_alleles_key:
        mt = mt.key_rows_by(locus=new_locus, alleles=mt.alleles)
    else:
        mt = mt.key_rows_by(locus=new_locus)

    
    return mt

import plotly.graph_objects as go
from plotly.subplots import make_subplots

def create_distribution_plots(mt, sequencingType, basename):
    """
    Creates interactive histogram plots for all QC metrics with threshold lines.
    :param mt: MatrixTable after hl.sample_qc() and run_charr()
    :param sequencingType: 'WES' or 'WGS'
    :param basename: prefix for output file
    """

    # collect metrics
    has_charr  = 'charr' in mt.col
    extra      = ['charr'] if has_charr else []
    rows       = mt.cols().key_by().select('s', 'sample_qc', *extra).collect()

    metrics = [
        {
            'name':           'Mean DP',
            'extractor':      lambda r: r.sample_qc.dp_stats.mean,
            'config_key':     f'DP_STATS.MEAN_{sequencingType}_threshold',
            'threshold_type': 'lower',
        },
        {
            'name':           'Ti/Tv ratio',
            'extractor':      lambda r: r.sample_qc.r_ti_tv,
            'config_key':     f'R_TI_TV_{sequencingType}_threshold',
            'threshold_type': 'interval',
        },
        {
            'name':           'Call Rate',
            'extractor':      lambda r: r.sample_qc.call_rate,
            'config_key':     'CALL_RATE_threshold',
            'threshold_type': 'lower',
        },
        {
            'name':           'N Singletons',
            'extractor':      lambda r: r.sample_qc.n_singleton,
            'config_key':     f'N_SINGLETON_{sequencingType}_threshold',
            'threshold_type': 'upper',
        },
        {
            'name':           'Het/Hom ratio',
            'extractor':      lambda r: r.sample_qc.r_het_hom_var,
            'config_key':     f'R_HET_HOM_VAR_{sequencingType}_threshold',
            'threshold_type': 'upper',
        },
        {
            'name':           'CHARR',
            'extractor':      lambda r: r.charr if has_charr else None,
            'config_key':     'CHARR_threshold',
            'threshold_type': 'upper',
        },
    ]

    ncols = 3
    nrows = int(np.ceil(len(metrics) / ncols))

    fig = make_subplots(
        rows=nrows, cols=ncols,
        subplot_titles=[m['name'] for m in metrics],
        vertical_spacing=0.15,
        horizontal_spacing=0.08
    )

    for i, mdef in enumerate(metrics):
        row = i // ncols + 1
        col = i %  ncols + 1

        # extract values
        vals = []
        for r in rows:
            try:
                v = mdef['extractor'](r)
                vals.append(float(v) if v is not None else np.nan)
            except Exception:
                vals.append(np.nan)

        vals = np.array(vals, dtype=float)
        valid_vals = vals[~np.isnan(vals)]

        if len(valid_vals) == 0:
            continue

        # histogram
        fig.add_trace(
            go.Histogram(
                x=valid_vals,
                nbinsx=40,
                marker_color='#4C9BE8',
                opacity=0.8,
                name=mdef['name'],
                showlegend=False,
                hovertemplate='Value: %{x:.4f}<br>Count: %{y}<extra></extra>',
            ),
            row=row, col=col
        )

        # threshold lines
        threshold = config['sample_filters'].get(mdef['config_key'])
        if threshold is not None:
            if mdef['threshold_type'] == 'lower':
                fig.add_vline(
                    x=threshold,
                    line_dash='dash', line_color='red', line_width=1.5,
                    row=row, col=col
                )
            elif mdef['threshold_type'] == 'upper':
                fig.add_vline(
                    x=threshold,
                    line_dash='dash', line_color='red', line_width=1.5,
                    row=row, col=col
                )
            elif mdef['threshold_type'] == 'interval':
                lower, upper = threshold
                fig.add_vline(
                    x=lower,
                    line_dash='dash', line_color='red', line_width=1.5,
                    row=row, col=col
                )
                fig.add_vline(
                    x=upper,
                    line_dash='dash', line_color='orange', line_width=1.5,
                    row=row, col=col
                )

    fig.update_layout(
        title=dict(
            text=f'Sample QC Metric Distributions — {sequencingType} | {basename}',
            font=dict(size=16)
        ),
        height=400 * nrows,
        width=1400,
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(size=11),
        bargap=0.05,
    )

    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=True, gridcolor='#eeeeee', title_text='Count')

    outpath = f'{basename}_sample_qc_distribution.html'
    fig.write_html(outpath)
    logging.info(f"Interactive distribution plot saved to: {outpath}")
    
def create_distribution_plots_points(mt, sequencingType, basename):
    """
    Creates interactive distribution plots (one point per sample) for all QC metrics,
    with hard thresholds drawn as lines and outlier samples tagged on hover.
    :param mt: MatrixTable after hl.sample_qc() and run_charr()
    :param sequencingType: 'WES' or 'WGS'
    :param basename: prefix for output file
    """

    # ── collect metrics ───────────────────────────────────────────────────────
    has_charr  = 'charr' in mt.col
    extra      = ['charr'] if has_charr else []
    rows       = mt.cols().key_by().select('s', 'sample_qc', *extra).collect()
    sample_ids = [r.s for r in rows]

    # ── metric definitions ────────────────────────────────────────────────────
    # (display_name, extractor, config_key_or_keys, lower_or_upper_or_interval)
    # threshold_type: 'lower', 'upper', 'interval'
    metrics = [
        {
            'name':           'Mean DP',
            'extractor':      lambda r: r.sample_qc.dp_stats.mean,
            'config_key':     f'DP_STATS.MEAN_{sequencingType}_threshold',
            'threshold_type': 'lower',
        },
        {
            'name':           'Ti/Tv ratio',
            'extractor':      lambda r: r.sample_qc.r_ti_tv,
            'config_key':     f'R_TI_TV_{sequencingType}_threshold',
            'threshold_type': 'interval',
        },
        {
            'name':           'Call Rate',
            'extractor':      lambda r: r.sample_qc.call_rate,
            'config_key':     'CALL_RATE_threshold',
            'threshold_type': 'lower',
        },
        {
            'name':           'N Singletons',
            'extractor':      lambda r: r.sample_qc.n_singleton,
            'config_key':     f'N_SINGLETON_{sequencingType}_threshold',
            'threshold_type': 'upper',
        },
        {
            'name':           'Het/Hom ratio',
            'extractor':      lambda r: r.sample_qc.r_het_hom_var,
            'config_key':     f'R_HET_HOM_VAR_{sequencingType}_threshold',
            'threshold_type': 'upper',
        },
        {
            'name':           'CHARR',
            'extractor':      lambda r: r.charr if has_charr else None,
            'config_key':     'CHARR_threshold',
            'threshold_type': 'upper',
        },
    ]

    # ── build subplots ────────────────────────────────────────────────────────
    n_metrics = len(metrics)
    ncols     = 3
    nrows     = int(np.ceil(n_metrics / ncols))

    fig = make_subplots(
        rows=nrows, cols=ncols,
        subplot_titles=[m['name'] for m in metrics],
        vertical_spacing=0.12,
        horizontal_spacing=0.08
    )

    for i, mdef in enumerate(metrics):
        row = i // ncols + 1
        col = i %  ncols + 1

        # extract values
        vals = []
        for r in rows:
            try:
                v = mdef['extractor'](r)
                vals.append(float(v) if v is not None else np.nan)
            except Exception:
                vals.append(np.nan)

        vals = np.array(vals, dtype=float)

        # get threshold from config
        threshold = config['sample_filters'].get(mdef['config_key'])

        # determine pass/fail per sample
        def is_outlier(v, threshold, threshold_type):
            if np.isnan(v) or threshold is None:
                return False
            if threshold_type == 'lower':
                return v < threshold
            elif threshold_type == 'upper':
                return v > threshold
            elif threshold_type == 'interval':
                lower, upper = threshold
                return v < lower or v > upper

        outlier_mask = np.array([
            is_outlier(v, threshold, mdef['threshold_type']) for v in vals
        ])

        x_pass    = np.where(~outlier_mask)[0]
        x_fail    = np.where( outlier_mask)[0]
        ids_pass  = [sample_ids[j] for j in x_pass]
        ids_fail  = [sample_ids[j] for j in x_fail]
        vals_pass = vals[~outlier_mask]
        vals_fail = vals[ outlier_mask]

        # passing samples
        fig.add_trace(
            go.Scatter(
                x=list(range(len(x_pass))),
                y=vals_pass,
                mode='markers',
                marker=dict(color='#4C9BE8', size=4, opacity=0.7),
                name='Pass',
                text=ids_pass,
                hovertemplate='<b>%{text}</b><br>Value: %{y:.4f}<extra></extra>',
                showlegend=(i == 0),
                legendgroup='pass'
            ),
            row=row, col=col
        )

        # failing samples
        if len(x_fail) > 0:
            fig.add_trace(
                go.Scatter(
                    x=list(range(len(x_fail))),
                    y=vals_fail,
                    mode='markers',
                    marker=dict(color='red', size=6, opacity=0.9, symbol='x'),
                    name='Outlier',
                    text=ids_fail,
                    hovertemplate='<b>%{text}</b><br>Value: %{y:.4f}<extra></extra>',
                    showlegend=(i == 0),
                    legendgroup='outlier'
                ),
                row=row, col=col
            )

        # threshold lines
        if threshold is not None:
            if mdef['threshold_type'] in ('lower', 'upper'):
                fig.add_hline(
                    y=threshold,
                    line_dash='dash', line_color='red', line_width=1.5,
                    row=row, col=col
                )
            elif mdef['threshold_type'] == 'interval':
                lower, upper = threshold
                fig.add_hline(
                    y=lower,
                    line_dash='dash', line_color='red', line_width=1.5,
                    row=row, col=col
                )
                fig.add_hline(
                    y=upper,
                    line_dash='dash', line_color='orange', line_width=1.5,
                    row=row, col=col
                )

    # ── layout ────────────────────────────────────────────────────────────────
    fig.update_layout(
        title=dict(
            text=f'Sample QC Metrics — {sequencingType} | {basename}',
            font=dict(size=16)
        ),
        height=400 * nrows,
        width=1400,
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(size=11),
    )

    fig.update_xaxes(showticklabels=False, showgrid=False)
    fig.update_yaxes(showgrid=True, gridcolor='#eeeeee')

    # ── save ──────────────────────────────────────────────────────────────────
    outpath = f'{basename}_sample_qc_distribution.html'
    fig.write_html(outpath)
    logging.info(f"Interactive distribution plot saved to: {outpath}")

    return sample_ids

import matplotlib.pyplot as plt
import seaborn as sns

def plot_sample_qc_histograms(collected, thresholds, sequencingType, basename):
    """
    Plot histograms of sample QC metrics with MAD-based thresholds overlaid.
    
    :param collected: dict of sample metrics, e.g. collected['dp_stats_mean'] = list of values
    :param thresholds: dict of MAD thresholds, e.g. thresholds['dp_stats_mean'] = (lower, upper)
    :param sequencingType: 'WES' or 'WGS' (used for plot titles)
    :param basename: string used for filenames
    """
    sns.set(style="whitegrid")
    
    for metric, values in collected.items():
        plt.figure(figsize=(8, 5))
        sns.histplot(values, bins=40, kde=False, color='skyblue')
        
        # Add thresholds if available
        thresh = thresholds.get(metric, None)
        if thresh is not None:
            if isinstance(thresh, tuple):
                lower, upper = thresh
                if lower is not None:
                    plt.axvline(lower, color='red', linestyle='--', label=f'Lower threshold: {lower:.2f}')
                if upper is not None:
                    plt.axvline(upper, color='green', linestyle='--', label=f'Upper threshold: {upper:.2f}')
            else:  # single-sided threshold
                plt.axvline(thresh, color='red', linestyle='--', label=f'Threshold: {thresh:.2f}')
        
        plt.title(f'{metric} distribution — {sequencingType} samples')
        plt.xlabel(metric)
        plt.ylabel('Number of samples')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{basename}_{metric}_histogram.png', dpi=150)
        plt.close()

def apply_filter(mt, condition, metric_path, metric_name, summary):
    """
    Applies a single sample QC filter step.
    :param mt: MatrixTable to filter
    :param condition: Hail boolean expression for filtering
    :param metric_path: string path to metric for create_stats e.g. "sample_qc.call_rate"
    :param metric_name: display name for logging and summary e.g. "Call Rate"
    :param summary: summary list to append results to
    :return: filtered MatrixTable
    """
    repo_stats = mt.aggregate_cols(hl.struct(
        total=hl.agg.count(),
        passing=hl.agg.count_where(condition)
    ))

    if config['verbosity']:
        stats = create_stats(mt, metric_path, "cols")
        summary.append([metric_name, repo_stats.total, repo_stats.passing,
                        repo_stats.total - repo_stats.passing, stats])

    logging.info(f"{metric_name} filtering done - Samples removed: {repo_stats.total - repo_stats.passing}")
    mt = mt.filter_cols(condition)
    return mt

def compute_mad_threshold(values, k=4.0, one_sided=None):
    """
    Computes MAD-based thresholds for a metric.
    :param values: numpy array of metric values
    :param k: MAD multiplier (gnomAD default = 4.0)
    :param one_sided: None (two-sided), 'lower' (flag below only), 'upper' (flag above only)
    :return: (lower, upper, median, mad)
    """
    median = float(np.median(values))
    mad    = float(np.median(np.abs(values - median)))

    # floor to avoid threshold collapse when MAD = 0
    if mad == 0:
        mad = 0.01 * abs(median) if median != 0 else 0.01

    spread = k * mad
    lower  = median - spread if one_sided != 'upper' else -np.inf
    upper  = median + spread if one_sided != 'lower' else  np.inf

    logging.info(f"MAD threshold | median={median:.4f}, MAD={mad:.4f}, "
                 f"lower={lower:.4f}, upper={upper:.4f}")

    return lower, upper, median, mad

def collect_metrics(mt):
    """
    Collects all QC metrics from the MatrixTable to numpy arrays in one Hail action.
    :param mt: MatrixTable after hl.sample_qc() and run_charr()
    :return: dictionary of {metric_name: numpy array}
    """
    has_charr = 'charr' in mt.col
    extra     = ['charr'] if has_charr else []
    rows      = mt.cols().key_by().select('s', 'sample_qc', *extra).collect()

    extractors = {
        #'call_rate':     lambda r: r.sample_qc.call_rate,
        'r_ti_tv':       lambda r: r.sample_qc.r_ti_tv,
        'r_het_hom_var': lambda r: r.sample_qc.r_het_hom_var
        #'n_singleton':   lambda r: r.sample_qc.n_singleton,
        #'dp_stats_mean': lambda r: r.sample_qc.dp_stats.mean,
        #'charr':         lambda r: r.charr if has_charr else None,
    }

    collected = {}
    for key, extractor in extractors.items():
        try:
            vals = [float(extractor(r)) if extractor(r) is not None else np.nan for r in rows]
            collected[key] = np.array(vals, dtype=float)
            n_valid = int(np.sum(~np.isnan(collected[key])))
            logging.info(f"Collected '{key}': {n_valid}/{len(rows)} valid values")
        except Exception as e:
            logging.warning(f"Could not collect metric '{key}': {e}")
            collected[key] = np.full(len(rows), np.nan)

    return collected

def apply_mad_filter(mt, metric_key, hail_expr, metric_path, metric_name,
                      summary, thresholds):
    """
    Applies a pre-computed MAD-based filter to the MatrixTable.
    :param mt: MatrixTable to filter
    :param metric_key: key to look up in thresholds dict e.g. 'call_rate'
    :param hail_expr: Hail expression for the metric e.g. mt.sample_qc.call_rate
    :param metric_path: string path for create_stats e.g. "sample_qc.call_rate"
    :param metric_name: display name for logging and summary
    :param summary: summary list to append results to
    :param thresholds: dict of pre-computed thresholds from _compute_mad_threshold
    :return: filtered MatrixTable
    """
    # Avoids breaking when metric is not in matrix 
    if thresholds[metric_key] is None:
        logging.warning(f"{metric_name} MAD filter skipped — threshold could not be computed")
        return mt                  # return mt unchanged
    
    lower, upper, median, mad = thresholds[metric_key]

    condition = (hail_expr >= lower) & (hail_expr <= upper)

    logging.info(f"{metric_name} MAD filter | lower={lower:.4f}, upper={upper:.4f}")

    mt = apply_filter(mt,
                       condition=condition,
                       metric_path=metric_path,
                       metric_name=metric_name,
                       summary=summary)
    return mt


def safe_mad_threshold(collected, key, k=4.0, one_sided=None):
    """
    Wraps _compute_mad_threshold with a NaN check.
    Returns None if the metric is unavailable, threshold tuple otherwise.
    """
    values = collected.get(key)
    if values is None:
        logging.warning(f"MAD | '{key}' not found in collected metrics — skipping")
        return None
    
    valid = values[~np.isnan(values)]
    if len(valid) < 2:
        logging.warning(f"MAD | '{key}' has fewer than 2 valid values — skipping")
        return None

    return compute_mad_threshold(valid, k=k, one_sided=one_sided)



    """
    Wraps _compute_mad_threshold with a NaN check.
    Returns None if the metric is unavailable, threshold tuple otherwise.
    """
    values = collected.get(key)
    if values is None:
        logging.warning(f"MAD | '{key}' not found in collected metrics — skipping")
        return None
    
    valid = values[~np.isnan(values)]
    if len(valid) < 2:
        logging.warning(f"MAD | '{key}' has fewer than 2 valid values — skipping")
        return None

    return _compute_mad_threshold(valid, k=k, one_sided=one_sided)

        
"""mt = hl.read_matrix_table(config['mt_from_vcf'])
create_stats_samples(mt, mt.sample_qc.dp_stats)"""