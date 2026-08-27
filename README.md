# VCF Processing Pipeline

This pipeline generates a VCF annotated with allele frequencies stratified by ancestry, sex and ancestry by sex. It is designed to facilitate downstream analyses that require population-aware frequency data. The pipeline consists of the following main modules:

- Variant, Genotype, and Sample Quality Control

- Sex Inference

- Ancestry Inference

- Allele Frequency Recalculation and Annotation by Ancestry and Sex

The final output will be a VCF file with all population allele frequencies (`AF`) and supporting metrics annotated.

**Bear in mind that the original AF values (if present in the input VCF) are NOT deleted, the output VCF will have the original AFs and the recalc ones.**

**Note:**  
If sex and ancestry cannot be inferred from genomic data, or ancestry information is not submitted, sex-based and ancestry-based grouping will be skipped.

## Requirements and Installation

## Option 1: Locally

### 1. Install Hail 

On a recent Debian-like system, the following should suffice

```
apt-get install -y \
    openjdk-11-jre-headless \
    g++ \
    python3.9 python3-pip \
    libopenblas-base liblapack3
python3.9 -m pip install hail
```

If more information is, please visit [Hail's documentation page](https://hail.is/docs/0.2/getting_started.html).

### 2. Install Required Python Packages

This workflow is written in Python and requires the following additional packages: 

- PyYAML
- pandas

You can install them using `pip`: 

```
pip install PyYAML pandas
```

### 3. CHARR Contamination Filtering Reference

If you plan to run CHARR contamination filtering (recommended), you’ll need to download the reference database for the genome version of your data. Bear in mind that, if all your data belongs to the same version, you'll only need to download that reference. 

* **If your data was aligned to GRCh37**

To check the size and contents of the folder before downloading:

```
gsutil ls -l gs://gcp-public-data--gnomad/release/2.1.1/ht/genomes/gnomad.genomes.r2.1.1.sites.ht/
```

To download the reference data: 

```
gsutil cp -r gs://gcp-public-data--gnomad/release/2.1.1/ht/genomes/gnomad.genomes.r2.1.1.sites.ht/ .
```

* **If your data was aligned to GRCh38**

To check the size and contents of the folder before downloading:

```
gsutil ls -l gs://gcp-public-data--gnomad/release/4.1/ht/genomes/gnomad.genomes.v4.1.sites.ht/
```

To download the reference data: 

```
gsutil cp -r gs://gcp-public-data--gnomad/release/4.1/ht/genomes/gnomad.genomes.v4.1.sites.ht/ .
```

### 4. Clone the repository 

Once all the requirements have been installed, clone this repository to your local computer or cluster. When cloning, you will download all the code for the pipeline, as well as the binaries for GrafAnc — the tool used to infer the ancestry of the samples.

## Option 2: Docker

See [README-for-docker.md](https://github.com/EGA-archive/AF_hail_pipeline/blob/main/README--for-docker.md) to use this pipeline inside a Docker container. 

## HOW TO RUN THE PIPELINE

All parameters and module executions are controlled via `config.yaml`. Example:

```yaml
## ============================================================
## HAIL PIPELINE CONFIGURATION FILE
## ============================================================
## Docker mount points:
##   Input files  -> /data/input
##   Output files -> /data/output
##   Work files   -> /data/work
##
## Suggested reading order:
##   1) General setup
##   2) Modules to run
##   3) Preprocessing substeps
##   4) Thresholds
##   5) Diagnostics
##   6) Spark / runtime
## ============================================================


## ============================================================
## 1. GENERAL SETUP
## Core paths and dataset-wide settings
## ============================================================

vcf_dir: ""
# Directory containing all input VCF files.
# IMPORTANT: all VCFs must use the same reference genome.
# Docker path: /data/input

vcf_for_header: ""
# Representative VCF used to build the final output VCF header.
# TODO: ideally remove this and use one file automatically from vcf_dir.
# Docker path: /data/input/filename.vcf.gz

ref_gen: ""
# Reference genome build used by the input VCFs.
# Accepted values: "GRCh37" | "GRCh38"

seq_type: ""
# Sequencing strategy.
# Affects sample QC thresholds in sample_filters.
# Accepted values: "WGS" | "WES"


## ============================================================
## 2. INTERMEDIATE AND FINAL OUTPUTS
## Files produced by the pipeline
## ============================================================

mt_from_vcf: ""
# MatrixTable generated during preprocessing.
# Docker path: /data/work/filename.mt

mt_afterQC: ""
# MatrixTable after QC steps.
# Docker path: /data/work/filename_afterQC.mt

final_vcf_AF: ""
# Final output VCF annotated with recalculated allele frequencies.
# Only relevant when af_annotation = true.
# Docker path: /data/output/your_output_name.vcf.bgz

summary_VCF: true
# If true, the output VCF will contain only variant-level information
# (no sample columns or genotype data).

use_allele_counts: true
# If true, the final VCF will include AC_hom, AC_het, and AC_hemi fields instead of genotye counts.
# AC_hom = 2*nhomalt, AC_het = nhet, AC_hemi = nhemi.


## ============================================================
## 3. MAIN MODULES TO RUN
## Turn complete modules on/off
## ============================================================
## Set to true to run a module, false to skip it.
##
## Notes:
## - prfalseeprocessing = false assumes mt_from_vcf already exists
## - infer_ancestry and submit_ancestry are mutually exclusive
## ============================================================

preprocessing: true
# Runs preprocessing: VCF conversion, splitting, and initial filtering.

infer_sex: true
# Infers biological sex from genotype data and flags mismatches.

delete_related: true
# Detects and removes related samples (duplicates / close relatives).

infer_ancestry: true
# Infer ancestry from genotype data using GrafAnc.
# Must be false if submit_ancestry = true.

submit_ancestry: false
# Read ancestry labels from ancestry_information file.
# Must be false if infer_ancestry = true.

ancestry_information: ""
# CSV with known ancestry labels.
# Only used when submit_ancestry = true.
# Required columns (comma-separated header): SampleID,Population

af_annotation: true
# Recalculate allele frequencies and annotate the final VCF.


## ============================================================
## 4. PREPROCESSING SUBSTEPS
## Only relevant when preprocessing = true
## ============================================================
## Typical first run:
##   convert_vcfs: true
##   split_multiallelic: true
## Then enable filtering steps as needed.
## ============================================================

convert_vcfs: true
# Convert input VCFs into a Hail MatrixTable.
# Must be true on first run.

split_multiallelic: true
# Split multiallelic variants into biallelic rows.

genotype_filtering: true
# Apply genotype-level filters using genotype_filters below.

variant_filtering: true
# Apply variant-level filters using variant_filters below.

sample_filtering: true
# Apply sample-level QC using sample_filters below.


## ============================================================
## 5. VARIANT FILTERING THRESHOLDS
## Used only if variant_filtering = true
## ============================================================
## To disable a filter, leave it commented like:
##   QD_threshold: #2.0
##   or set it to false
## ============================================================

variant_filters:
  QD_threshold: 2.0
  # Quality by Depth. Low values suggest low-confidence variants.

  DP_threshold: 15
  # Minimum total read depth across the cohort at a site.

  QUAL_threshold: 0.4
  # Minimum Phred-scaled variant quality score.

  MQ_threshold: 40
  # Mapping quality. Filters poorly aligned sites.

  FS_threshold: 60
  # Fisher Strand bias. High values suggest strand-specific artifacts.

  READPOSRANKSUM_threshold: -8.0
  # Filters variants where alt alleles cluster near read ends.


## ============================================================
## 6. GENOTYPE FILTERING THRESHOLDS
## Used only if genotype_filtering = true
## ============================================================
## To disable a filter, leave it commented.
## ============================================================

genotype_filters:
  GQ_threshold: 20
  # Genotype Quality. Calls below this are set to missing.

  AB_threshold: 0.2
  # Allele Balance for heterozygous calls.
  # Calls below this threshold are set to missing.


## ============================================================
## 7. SAMPLE FILTERING THRESHOLDS
## Used only if sample_filtering = true
## ============================================================
## Some thresholds depend on seq_type (WGS vs WES).
## To disable a filter, leave it commented.
## ============================================================

sample_filters:
  DP_STATS.MEAN_WGS_threshold: 15
  # Minimum mean read depth for WGS samples.

  DP_STATS.MEAN_WES_threshold: 10
  # Minimum mean read depth for WES samples.

  CALL_RATE_threshold: 0.95
  # Minimum fraction of non-missing genotype calls per sample.

  N_SINGLETON_WGS_threshold: 100000
  # Maximum singleton count for WGS.

  N_SINGLETON_WES_threshold: 5000
  # Maximum singleton count for WES.

  CHARR_threshold: 0.05
  # Contamination estimate threshold.

  R_TI_TV: true
  # Expected Ti/Tv ratio. Thresholds calculated as 4 median absolute deviations from Ti/Tv ratio of the cohort.

  R_HET_HOM_VAR: true
  # Maximum heterozygous / homozygous-variant ratio for WGS. Thresholds calculated as 4 median absolute deviations from het/hom ratio of the cohort.


## ============================================================
## 8. EXTERNAL RESOURCES
## Optional files needed by some QC calculations
## ============================================================

gnomad_sites_GRCh37: ""
# Path to gnomAD sites Hail Table for GRCh37.
# Used for CHARR computation.
# Docker path: /data/input/gnomad.genomes.r2.1.1.sites.ht

gnomad_sites_GRCh38: ""
# Path to gnomAD sites Hail Table for GRCh38.
# Used for CHARR computation.
# Docker path: /data/input/gnomad.genomes.v4.1.sites.ht


## ============================================================
## 9. LOGGING AND QC DIAGNOSTICS
## Optional outputs for auditing and visual inspection
## ============================================================

verbosity: true
# If true, generate a CSV tracking how many variants are removed at each QC step.

plots: true
# If true, generate sample-level QC plots.


## ============================================================
## 10. SPARK / HAIL RUNTIME SETTINGS
## Used when running on a cluster or with custom Spark settings
## ============================================================
## If cluster = false, default Spark settings are used.
## ============================================================

cluster: False
# Set to True to apply the custom Spark configuration below.

spark_driver_memory: "50g"
# Memory for Spark driver.

spark_executor_memory: "20g"
# Memory for each Spark executor.

spark_executor_instances: "4"
# Number of executors.

spark_executor_cores: "4"
# CPU cores per executor.

spark_rpc_askTimeout: "300s"
# Timeout for Spark RPC operations.

spark_sql_shuffle_partitions: "200"
# Number of partitions used during shuffles.

spark_network_timeout: "800s"
# Maximum Spark network wait time.

spark_local_dir: "/work/tmp"
# Spark temporary directory.
# Docker path: /data/work/tmp

tmp_dir: "/work/tmp"
# Hail temporary directory.
# Docker path: /data/work/tmp

local_tmpdir: "/work/tmp"
# Local temporary directory for Hail filesystem operations.
# Docker path: /data/work/tmp
```

---

**Modular Design:**  

Each module and function can be run independently:
- If a module is set to False, it will be skipped
- If a preprocessing step is set to False it will be skipped
- If a filtering step has the threshold value commented it will be skipped
- Thresholds can be adjusted via `config.yaml`. If a threshold is not set (threshold value commented), that filtering step will be skipped. 

---


**Once the `conf.yaml` is adjusted to your needs, you only need to run:**

```
python vcf-af-pipeline.py
```
In the pipeline diagram you'll find the different paths your data can follow with this pipeline. In purple you'll find highlighted our proposed path, where all the quality control steps are performed, related samples are deleted and ancestry is inferred. 

![Pipeline Diagram](VCF_EGA_pipeline.png)

## REFERENCES

jimmy-penn/grafanc: GrafPop from dbSNP [GitHub]. [cited 2025 Jul 24]. Available from: https://github.com/jimmy-penn/grafanc/tree/master

Hail Team. Hail 0.2. [cited 2025 Jul 24]. Available from: https://github.com/hail-is/hail
