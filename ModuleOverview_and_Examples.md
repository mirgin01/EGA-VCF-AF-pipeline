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

*Note:* The VCF must include read-level information (DP, GQ, MQ, FS, etc.) for these filters to work.

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

