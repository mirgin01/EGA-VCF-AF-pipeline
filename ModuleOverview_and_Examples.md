## Module Overview

### MODULE 1: Preprocessing

#### Step 1: Convert VCF to Hail MatrixTable

Convert all VCF files in a folder into a single Hail MatrixTable for downstream processing.

#### Step 2: Split Multi-Allelic Variants

#### Step 3: Variant Quality Control

Apply the following quality filters:

| Metric                   | Threshold         | Description                                                                                 |
|---------------------------|--------------------|---------------------------------------------------------------------------------------------|
| Quality by Depth (`QD`)   | `< 2.0`            | Variant quality normalized by depth.                                                       |
| Depth of Coverage (`DP`)  | `< 15`             | Total sequencing reads supporting the position.                                            |
| Variant Quality (`QUAL`)  | `< 30`             | Confidence score for the variant.                                                          |
| Mapping Quality (`MQ`)    | `< 40`             | Read alignment confidence.                                                                 |
| Fisher Strand Bias (`FS`) | `> 60`             | Measures strand bias in sequencing reads.                                                  |
| Read Position Bias        | `< -8.0`           | Measures whether alleles occur at read ends (potential bias).                              |

*Note:* The VCF must include read-level (INFO column) information (QD, DP, QUAL, MQ, FS and ReadPosRankSum) for these filters to work.

#### Step 4: Genotype Quality Control

Apply the following quality filters to mask out low-quality genotypes:

| Metric                   | Threshold         | Description                                                                                 |
|---------------------------|--------------------|---------------------------------------------------------------------------------------------|
| Genotype Quality (`GQ`)   | `< 20`             | Confidence in the assigned genotype. 
| Allele Balance (`AB`)     | `< 0.2`            | Ratio of alternative reads to total reads.                                                 |                                                      |

*Note:* The VCF must include genotype-level information (GQ and AD) for these filters to work.

#### Step 5: Sample Quality Control

| Metric                                | Threshold                      | Description                                                                                   |
| ------------------------------------- | ------------------------------ | --------------------------------------------------------------------------------------------- |
| Minimum Coverage                      | `WGS < 15 ; WES < 10`            | Minimum acceptable read depth per sample; ensures sufficient coverage.                        |
| Transition/Transversion Ratio (Ti/Tv) | `WGS: \~2.0–2.1; WES: \~3.0–3.3` | Expected ratio of transitions to transversions, indicating variant quality.                   |
| Het/Hom Ratio                         | `WGS > 3.3; WES > 10`            | Ratio of heterozygous to homozygous variants; detects abnormal variant patterns.              |
| Call Rate                             | `< 95%`                       | Proportion of variants successfully genotyped; low values suggest poor-quality samples.       |
| Singletons                            | `WGS > 100k; WES > 5k `          | Number of variants found only once in the dataset |
| Contamination (CHARR)                 | `WGS > 5%; WES > 0.015%  `       | Estimate of contamination in the sample.                  |

*Note:* The VCF must include sample-level (FORMAT column) information (DP,)

#### Step 6: Sex Inference

Infer sample sex using its genomic information. If sex cannot be determined, sex-based grouping will be skipped.

**Note:** Currently, only XX and XY sexes are inferred. Aneuploides will land under undefined.  

### MODULE 2: Delete related samples

Including related individuals in a cohort can distort AF calculations

### MODULE 3: Ancestry

#### Ancestry Inference

**Subset Hail Matrix for Ancestry SNPs**  
Extract ~282,424 ancestry-informative SNPs (if available).

**Run GRAF-Anc**  
Use GRAF-Anc for ancestry assignment.

**Process GRAF-Anc Results**  
GRAF-Anc assigns ancestry at two levels:

- **Continental**
- **Subcontinental**

Only the **continental level** will be used to ensure reliable population tagging, aligned with gnomAD v4 practices.

| Code | Description                  | gnomAD Equivalent                |
|------|------------------------------|----------------------------------|
| AFR  | African                      | African / African American       |
| MEN  | Middle East and North Africa | Middle Eastern                   |
| EUR  | European                     | Non-Finnish European             |
| SAS  | South Asian                  | South Asian                      |
| EAS  | East Asian                   | East Asian                       |
| AMR  | Admixed American             | Admixed American                 |
| OCN  | Oceania                      | —                                |
| MIX  | Multi-ancestry               | Remaining samples                |

*Note:* GRAF-Anc does **not** tag: Ashkenazi Jewish, Amish, or Finnish populations.

If there are not enought ancestry-informative SNPs in the input VCF the ancestry-based grouping will be skipped. 

### MODULE 4: Allele Frequency Recalculation

It calculates the allele frequencies statistics for the total entries of the dataset and stratified by male, female an each ancestry group. 
It outputs: 
- Allele Count (AC) : alternate allele count in high quality genotypes
- Allele Number (AN) : total number of called high quality genotypes
- Allele Frequency (AF) : number of individuals for alternate allele
- Number of homozygotes (nhomalt) : alternate allele frequency in high quality genotypes


And finally, export a single VCFs with all the AF fields annotated. 

## Examples 

### 1. Obtain Input Files

To test the pipeline, first retrieve the necessary VCF files:

Option A: Install and use [pyega3](https://ega-archive.org/access/download/files/pyega3/)

Option B: Use the [live outbox](https://ega-archive.org/access/download/files/live-outbox/)


Example with pyega3:

```
pyega3 -d -t datasets

pyega3 -d -t fetch EGAF00007243776
```

### 2. Running the Pipeline

When running with the provided conf.yaml, the following log is produced:

```
2025-10-03 09:57:57,643 [INFO]: LOG SAVED IN: testing-pipeline_20251003_095750.log
2025-10-03 09:57:57,667 [INFO]: === Pipeline settings ===
...
2025-10-03 09:58:14,977 [INFO]: Total number of entries with GT information: 0
2025-10-03 09:58:14,977 [ERROR]: GQ information not available - terminating pipeline
Traceback (most recent call last):
  ...
ValueError: GQ information is not available. To proceed with filtering, add a hash sign to the threshold value to inactivate the metric and perform quality control without this step.
```

### 3. Understanding the Error

The error indicates that genotype quality (GQ) is missing from the input VCF. For the genotype filtering module to run, your VCF must include both:

- GQ (Genotype Quality)

- AD (Allele Depth)

If these fields are missing, the pipeline cannot perform genotype-based QC. In this case, you must disable genotype filtering in the configuration.

Additionally, since INFO fields related to variant quality are also absent, only the QUAL filter can be applied.

Regarding sample filtering, all genotypes are missing (./.), so sample-based QC must also be disabled.

The adapted conf.yaml for EGAF00007243776 twould be:

```
## PREPROCESSING STEPS
convert_vcfs : false   # disable VCF conversion
split_multiallelic : true
genotype_filtering : false
variant_filtering : true
sample_filtering : false

# Note: To deactivate a specific threshold, comment it out with a hash (#).

## VARIANT FILTERING THRESHOLDS
variant_filters:
  QD_threshold : #2.0
  DP_threshold : #15
  QUAL_threshold : 40
  MQ_threshold : #40
  FS_threshold : #60
  READPOSRANKSUM_threshold : #-8.0

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
  R_TI_TV_WES_threshold : [3.0, 3.3]
  R_TI_TV_WGS_threshold : [2.0, 2.1]

gnomad_sites_GRCh37 : ""   # reference for GRCh37
gnomad_sites_GRCh38 : ""   # reference for GRCh38
```
