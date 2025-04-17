# VCF Processing Workflow for the EGA — Creating a Standard Format

This workflow aims to standardize VCF files for the European Genome-phenome Archive (EGA).  
The resulting **VSE** (VCF Standard EGA) files will have undergone:

- Variant and sample quality control
- Metadata sanity checks (sex and ancestry)
- Allele frequencies recalculation and annotation

The final output will be a VCF file with all population allele frequencies (`AF`) and supporting metrics annotated:

```
AF_total
AC_total
AC_hom_total
AN_total
AF_male
AC_male
AC_hom_male
AN_male
AF_female
AC_female
AC_hom_female
AN_female
AF_pop_n
AC_pop_n
AC_hom_pop_n
AN_pop_n
```

**Note:**  
If sex cannot be inferred from genomic data, sex-based grouping will be skipped.

---

## Module Overview

### MODULE 1: Preprocessing

#### Step 1 — Convert VCF to Hail MatrixTable

Convert all VCF files in a folder into a single Hail MatrixTable for downstream processing.

---

#### Step 2 — Split Multi-Allelic Variants

When exporting to VCF, only the allele frequency (`AF`) of the **first alternative allele** is typically annotated.  
Hail by default exports all AF values per variant, e.g.:

```
Variant: 1:1234:A-T  
AF_total = [0.99, 0.01]  # [reference, alternative]
```

Since the Beacon protocol only accepts a single AF per variant, multi-allelic variants are **split**. After splitting, the AF for the alternative allele is retained, ensuring accurate Beacon compatibility.

---

#### Step 3 — Variant Quality Control

Apply the following quality filters:

| Metric                   | Threshold         | Description                                                                                 |
|---------------------------|--------------------|---------------------------------------------------------------------------------------------|
| Quality by Depth (`QD`)   | `< 2.0`            | Variant quality normalized by depth.                                                       |
| Depth of Coverage (`DP`)  | `< 15`             | Total sequencing reads supporting the position.                                            |
| Variant Quality (`QUAL`)  | `< 30`             | Confidence score for the variant.                                                          |
| Genotype Quality (`GQ`)   | `< 20`             | Confidence in the assigned genotype.                                                       |
| Mapping Quality (`MQ`)    | `< 40`             | Read alignment confidence.                                                                 |
| Fisher Strand Bias (`FS`) | `< 40`             | Measures strand bias in sequencing reads.                                                  |
| Read Position Bias        | `< -8.0`           | Measures whether alleles occur at read ends (potential bias).                              |
| Allele Balance (`AB`)     | `< 0.2`            | Ratio of alternative reads to total reads.                                                 |

*Note:* The VCF must include read-level information (DP, GQ, MQ, FS, etc.) for these filters to work.

---

#### Step 4 — Sample Quality Control

| Metric                     | Threshold                            | Notes                                                |
|-----------------------------|--------------------------------------|------------------------------------------------------|
| Minimum Coverage            | `WGS < 15` ; `WES < 10`              | Computed with `hl.sample_qc()` if `DP` or `MIN_DP` exist. |
| Transition/Transversion Ratio (`Ti/Tv`) | `WGS: ~2.0-2.1`; `WES: ~3.0-3.3` | Computed with `hl.sample_qc()`. |
| Het/Hom Ratio               | `WGS > 3.3`; `WES > 10`              | Computed with `hl.sample_qc()`. |
| Call Rate                   | `< 95%`                              | Computed with `hl.sample_qc()`. |
| Singletons                  | `> 2 SD` or `WGS > 100k`, `WES > 5k` | Computed with `hl.sample_qc()`. |
| Contamination (CHARR)       | `WGS > 5%`, `WES > 0.015%`           | Estimated using `hl.compute_charr()` *(TODO)*.        |

---

### MODULE 2: Ancestry

#### Step 5 — Ancestry Inference

**5.1 Subset Hail Matrix for Ancestry SNPs**  
Extract ~282,424 ancestry-informative SNPs (if available).

**5.2 Run GRAF-Anc**  
Use GRAF-Anc for ancestry assignment.

**5.3 Process GRAF-Anc Results**  
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

---

### MODULE 3: Allele Frequency Recalculation

#### Step 6 — Sex Inference

Infer sample sex using:

```python
hl.impute_sex()
```

If sex cannot be determined, sex-based grouping will be skipped.

---

#### Step 7 — Ancestry Annotation

Annotate the VCF with GRAF-Anc ancestry assignments and recalculate allele frequencies (`AF`) per ancestry group.

---

#### Step 8 — AF Recalculation

Calculate `AF`, `AC`, `AC_hom`, and `AN` for:

- Total population
- Males
- Females
- Each ancestry group

---

#### Step 9 — Export Hail Matrix to VCF

Convert the processed MatrixTable back into a VCF, with allele frequencies fully annotated.

---

## Configuration

All parameters are controlled via `config.yaml`. Example:

```yaml
## TODO decide title
vcf_dir : "/home/mireia/Bioinfo/Beacon/beacon2-ri-tools-v2/files/vcf/files_to_read" # all the VCFs in this folder will be converted into a Hail matrix
vcf_for_header : "/home/mireia/Bioinfo/Beacon/beacon2-ri-tools-v2/files/vcf/files_to_read/GCAT-EGAD00001007774-af_annotated.vcf.bgz"  # the final VCF will have parts of this header
ref_gen : "GRCh37" # reference genome from the VCFs
mt_from_vcf : "/home/mireia/Bioinfo/Hail_folder/GCAT.mt" # path where the original matrix will be saved
seq_type : "WGS" # sequencing typw
mt_afterQC : "/home/mireia/Bioinfo/Hail_folder/GCAT_afterQC.mt" # path where the after QC matrix will be saved

## MODULES TO RUN
preprocessing : false # if true the module will be run
ancestry : true
af_annotation : true

## FUNCTIONS TO RUN
convert_vcfs : true # if false the module will be run
split_multiallelic : true
variant_filtering : true
sample_filtering : true

## VARIANT FILTERING THRESHOLDS
QD_threshold : 2.0 # threshold used during the QC
DP_threshold : 15 
QUAL_threshold : 30
MQ_threshold : 40
FS_threshold : 40
ReadPosRankSum_threshold : -8.0
GQ_threshold : 20
AB_threshold : 0.2

## SAMPLE FILTERING THRESHOLDS
DP_WGS_threshold : 15
DP_WES_threshold : 10
TITV_WGS_threshold : [2.0 , 2.1]
TITV_WES_threshold : [3.0 , 3.3]
callRate_threshold : 0.95
singletons_WGS_threshold : 5000
singletons_WES_threshold : 100000
hethom_WES_threshold : 10
hethom_WGS_threshold : 3.3

## ANCESTRY
ancestrySNPs : "/home/mireia/GitHub/Hail/GrafAnc_SNPs/" # update with your local path of GrafAnc_SNPs

## AF RECALC
final_vcf_QC_AF : "/home/mireia/Bioinfo/Hail_folder/GCAT_AF.vcf.bgz" # path for VCF annotated with AFs
```

---

✅ **Modular Design:**  
Each module and function can be run independently. Thresholds can be adjusted via `config.yaml`.

---



