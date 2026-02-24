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
## PATHS
vcf_dir : " "           # All VCFs in this folder must be from the same reference genome
vcf_for_header : ""     # The final VCF will have parts of this header
ref_gen : " "           # Reference genome of the VCFs (OPTIONS: GRCh37 / GRCh38)
seq_type : " "          # Sequencing type (OPTIONS: WGS / WES)
mt_from_vcf : " "       # Path where the MatrixTable produced by preprocessing will be saved
mt_afterQC : " "        # Path where the post-QC MatrixTable will be saved
final_vcf_AF : " "      # Path for the output VCF annotated with recalculated AFs
summary_VCF : false     # If true, the output VCF will contain no sample or genotype information

## MODULES TO RUN
preprocessing : true   # Run preprocessing module if true
infer_sex : true       # Run sex inference module if true
delete_related: true   # Remove related samples if true

# Ancestry analysis options:
# - infer_ancestry and submit_ancestry are mutually exclusive.
# - If infer_ancestry = true, submit_ancestry MUST be false.
# - Both can be false to skip ancestry processing entirely.
infer_ancestry : true
submit_ancestry: false
ancestry_information : ""  # Path to CSV with ancestry labels. Required header: SampleID\tPopulation

af_annotation : true  # Run allele frequency annotation if true

## PREPROCESSING STEPS
convert_vcfs : true        # Convert input VCFs to a Hail MatrixTable
split_multiallelic : true  # Split multiallelic sites into biallelic rows
genotype_filtering : true  # Apply genotype-level filters
variant_filtering : true   # Apply variant-level filters
sample_filtering : true    # Apply sample-level QC filters

# Note: to disable a specific filter, comment out its threshold value
# e.g.  QD_threshold : #2.0

## VARIANT FILTERING THRESHOLDS
variant_filters:
  QD_threshold : 2.0
  DP_threshold : 15
  QUAL_threshold : 40
  MQ_threshold : 40
  FS_threshold : 60
  READPOSRANKSUM_threshold : -8.0

## GENOTYPE FILTERING THRESHOLDS
genotype_filters:
  GQ_threshold : 20
  AB_threshold : 0.2

## SAMPLE FILTERING THRESHOLDS
sample_filters:
  DP_STATS.MEAN_WGS_threshold : 15
  DP_STATS.MEAN_WES_threshold : 10
  CALL_RATE_threshold : 0.95
  R_HET_HOM_VAR_WGS_threshold : 3.3
  R_HET_HOM_VAR_WES_threshold : #0
  N_SINGLETON_WGS_threshold : 100000
  N_SINGLETON_WES_threshold : 5000
  CHARR_threshold : 0.05
  R_TI_TV_WES_threshold : #[3.0 , 3.3]
  R_TI_TV_WGS_threshold : #[2.0 , 2.1]
gnomad_sites_GRCh37 : ""  # Path to gnomAD sites Hail Table for GRCh37
gnomad_sites_GRCh38 : ""  # Path to gnomAD sites Hail Table for GRCh38

## LOGs
verbosity : true  # If true, a CSV with variants removed per step will be created (increases execution time)
plots: false      # If true, box plots showing the distribution of each QC sample parameter will be createdpo)

## SPARK CONFIGURATION
# Allocating available memory to Spark helps avoid crashes when working with large datasets.
cluster : False  # If true, the Spark configuration below will be applied
spark_driver_memory: "50g"
spark_executor_memory : "20g"
spark_executor_instances: "4"
spark_executor_cores: "4"
spark_rpc_askTimeout: "300s"
spark_sql_shuffle_partitions: "200"
spark_network_timeout: "800s"
spark_local_dir: "./tmp"
tmp_dir: "./tmp"
local_tmpdir: "./tmp"
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
