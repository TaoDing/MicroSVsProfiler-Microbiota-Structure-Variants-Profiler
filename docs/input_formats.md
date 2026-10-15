# Input Data Formats
This document describes the input files of the MicroSVsProfiler workflow. Replace them with the names of your own files as appropriate.
## 1. Input file summary
| Input | Expected content | Format / orientation |
| --- | --- | --- |
| Raw metagenomic reads | Sequencing reads for individual samples | FASTQ |
| Clinical metadata | Clinical phenotypes associated with each sample | Tab-delimited; samples as rows, variables as columns |
| Genome annotations | Annotations for more than 3,000 representative microbial species | GTF |
| Metabolite abundance | Abundance measurements for metabolites | Tab-delimited |
| dSV profile | SGV-Finder deletion-like SV feature values | CSV; SV features as rows, samples as columns |
| vSV profile | SGV-Finder variable SV feature values | CSV; SV features as rows, samples as columns |
| Taxonomy table | Microbial taxonomy reference | Tab-delimited |
## 2. Input details

### 2.1 Raw metagenomic sequencing data
Provide FASTQ files containing the raw shotgun metagenomic reads for each sample.
### 2.2 Clinical metadata (`data/metadata.txt`)
The metadata table contains one row per individual/sample and one column per clinical variable. Sample identifiers are described in the document as SRR/DRR accessions. The R example uses a column called `Disease_Stage` with labels `HC` and `CRC`.
### 2.3 Genome annotation data
A GTF annotation resource is used for more than 3,000 representative microbial species.
### 2.4 Metabolite abundance (`data/metabolomic.tab`)
The document describes an abundance matrix containing **325 metabolites per individual**. No additional normalization is performed at the import stage, so values are assumed to be comparable across samples.
### 2.5 SV profiles (`data/ALL.dsgv.csv`, `data/ALL.vsgv.csv`)
These files contain normalized coverage-based feature profiles output by SGV-Finder. In the documentation, **rows are SV features and columns are samples**. They are read as comma-separated matrices with the first column used as row names.

The documentation uses `DSGV` and `VSGV`; the revised manuscript uses **dSV** (deletion-like SV region) and **vSV** (variable SV region). Feature identifiers should not be interpreted as precise genomic coordinates or breakpoint positions without supporting annotation.
### 2.6 Taxonomy reference (`data/MGS_taxonomy.csv`)

## 3. R import example
The following code preserves the input paths and delimiter choices, with explicit metadata row names and basic sample-ID checks. Adjust the metadata grouping column to match your input file.
```r
options(stringsAsFactors = FALSE)

phenotypes <- read.table(
  "data/metadata.txt", header = TRUE, sep = "\t",
  row.names = 1, check.names = FALSE
)

# Grouping convention:
if ("Disease_Stage" %in% names(phenotypes)) {
  ctrl <- rownames(phenotypes)[phenotypes$Disease_Stage == "HC"]
  crc  <- rownames(phenotypes)[phenotypes$Disease_Stage == "CRC"]
}

metabolomic <- read.table(
  "data/metabolomic.tab", header = TRUE, sep = "\t",
  row.names = 1, check.names = FALSE
)

mgs_taxonomy <- read.table(
  "data/MGS_taxonomy.csv", header = TRUE, sep = ",",
  row.names = 1, check.names = FALSE
)

df_dSV <- read.table(
  "data/ALL.dsgv.csv", header = TRUE, sep = ",",
  row.names = 1, check.names = FALSE
)

df_vSV <- read.table(
  "data/ALL.vsgv.csv", header = TRUE, sep = ",",
  row.names = 1, check.names = FALSE
)

# Example checks: SV matrices should have samples in columns.
stopifnot(!anyDuplicated(rownames(phenotypes)))
stopifnot(!anyDuplicated(colnames(df_dSV)))
stopifnot(!anyDuplicated(colnames(df_vSV)))

cat("Metadata samples:", nrow(phenotypes), "\n")
cat("dSV features / samples:", nrow(df_dSV), ncol(df_dSV), "\n")
cat("vSV features / samples:", nrow(df_vSV), ncol(df_vSV), "\n")
cat("dSV samples matched to metadata:",
    sum(colnames(df_dSV) %in% rownames(phenotypes)), "\n")
cat("vSV samples matched to metadata:",
    sum(colnames(df_vSV) %in% rownames(phenotypes)), "\n")
```

## 4. Items to verify before running
- Check each file's actual delimiter, header, first identifier column, and numeric-data orientation.
- Verify that sample IDs are consistent between metadata, raw sequencing files, and feature matrices.
- Record missing metadata and any cohort-specific exclusions before analyses.
- Confirm normalization and sample orientation of the metabolomics matrix before integrative analyses.

For software installation, see [Installation](installation.md).
