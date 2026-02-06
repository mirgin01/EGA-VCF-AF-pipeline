# VCF Processing Pipeline

This pipeline generates a VCF annotated with allele frequencies stratified by ancestry and sex. It is designed to facilitate downstream analyses that require population-aware frequency data. The pipeline consists of the following main modules:

- Variant, Genotype, and Sample Quality Control

- Sex Inference

- Ancestry Inference

- Allele Frequency Recalculation and Annotation by Ancestry and Sex

The final output will be a VCF file with all population allele frequencies (`AF`) and supporting metrics annotated:

```
AF_total_recalc
AC_total_recalc
nhomalt_total_recalc
AN_total_recalc
AF_male_recalc
AC_male_recalc
nhomalt_male_recalc
AN_male_recalc
AF_female_recalc
AC_female_recalc
nhomalt_female_recalc
AN_female_recalc
AF_pop_n_recalc
AC_pop_n_recalc
nhomalt_pop_n_recalc
AN_pop_n_recalc
```

**Bear in mind that the original AF values (if present in the input VCF) are NOT deleted, the output VCF will have the original AFs and the recalc ones.**

**Note:**  
If sex and ancestry cannot be inferred from genomic data, sex-based and ancestry-based grouping will be skipped.

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
## Option 2: Docker

See ADD LINK to use this pipeline inside a Docker container. 

### 4. Clone the repository 

Once all the requirements have been installed, clone this repository to your local computer or cluster. When cloning, you will download all the code for the pipeline, as well as the binaries for GrafAnc — the tool used to infer the ancestry of the samples.

## HOW TO RUN THE PIPELINE

All parameters and module executions are controlled via `config.yaml`. Example:

```yaml
## PATHS

vcf_dir : " " # all the VCFs in this folder will be converted to a Hail matrix. They must be from the same reference genome.
vcf_for_header : ""  # the final VCF will have parts of this header
ref_gen : " " # reference genome from the VCFs (OPTIONS: GRCh37 / GRCh38)
mt_from_vcf : " " # path where the original matrix will be saved
seq_type : " " # sequencing type (OPTIONS: WGS / WES)
mt_afterQC : " " # path where the after QC matrix will be

## LOGs
verbosity : true # if true a csv with variants deleted per step will be create. This increases the execution time.
plots: false # create box plot showing the distribution of each QC sample parameter

## MODULES TO RUN
preprocessing : true # if true the module will be run
delete_related: true
ancestry : true
af_annotation : true

## PREPROCESSING STEPS
convert_vcfs : true # if true the module will be run
split_multiallelic : true
genotype_filtering : true
variant_filtering : true
sample_filtering : true

## VARIANT FILTERING THRESHOLDS
variant_filters:
  QD_threshold : 2.0 # threshold used during QC
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
  R_HET_HOM_VAR_WES_threshold : 10 
  R_HET_HOM_VAR_WGS_threshold : 3.3
  N_SINGLETON_WGS_threshold : 100000
  N_SINGLETON_WES_threshold : 5000 
  CHARR_threshold : 0.05
  R_TI_TV_WES_threshold : [3.0 , 3.3]
  R_TI_TV_WGS_threshold : [2.0 , 2.1]
gnomad_sites_GRCh37 : "" # reference for GRCh37 
gnomad_sites_GRCh38 : "" # reference for GRCh38

## ANCESTRY
ancestrySNPs : "path/to/GrafAnc_SNPs/" # update with your local path of GrafAnc_SNPs (Downloaded when cloning this repo)

## AF RECALC
final_vcf_AF : " " # path for VCF annotated with AFs
summary_VCF : false # if true the output VCF wont contain ANY information about the samples and their genotypes

## SPARK CONFIGURATION 

# To work with big datasets allocating the available memory into spark avoids crashes. 

spark_diver_memory: "50g" # Allocate sufficient memory for the driver
spark_executor_memory : "20g" # Allocate memory for executors
spark_executor_memory: 4 # Use multiple executors
spark_executor_cores: 4 # Assign cores per executor
spark_rpc_askTimeout': "300s" # Increase timeout for slow operations
spark_sql_shuffle_partitions: 200 # Reduce shuffle partitions for large data
spark_memory_fraction: 0.8 # Use most of the JVM heap for Spark execution
spark_local_dir: "./tmp" # Specify a temp directory for disk spill
spark_network_timeout: 800s # Avoid network timeouts
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
