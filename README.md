# MicroSVsProfiler-Microbiota-Structure-Variants-Profiler
Based on previous software and advanced analysis process of microbial SV, we synthesized and implemented a pipeline called MicroSVsProfiler (Microbiota Structure Variants Profiler), which can perform preprocessing of metagenomic data, detect structural variations.
## 1. Overview
MicroSVsProfiler is a computational workflow for characterizing microbial genomic structural variation patterns from shotgun metagenomic sequencing data and investigating their associations with colorectal cancer (CRC).

The workflow integrates shell, R, and Python scripts for metagenomic preprocessing, microbial taxonomic and functional profiling, structural variation detection, statistical analysis, and diagnostic modeling.

This repository accompanies the study *Characterizing Gut Microbial Structural Variation Reveals Strain Heterogeneity and Metabolic Dysfunction in Colorectal Cancer* and provides scripts for reproducing the major computational analyses and adapting the workflow to other metagenomic datasets.

**Scope of SV detection**: MicroSVsProfiler uses SGV-Finder to identify coverage-defined deletion-like and variable genomic regions relative to reference genomes. These features do not necessarily represent breakpoint-resolved structural variants.
## 2. Workflow
The analytical workflow includes the following major stages:
1. Metagenomic quality control and host-read removal
2. Microbial taxonomic and functional profiling
3. Microbial genomic SV detection using SGV-Finder
4. SV-based distance calculation, PCoA and clustering
5. Differential SV and species abundance analysis
6. Hierarchical clustering and cohort heterogeneity analysis
7. Random forest diagnostic modeling and leakage-free validation
8. Microbial SV functional annotation and metabolite associations

Detailed scripts and step-specific instructions are provided in this repository.
## 3. Computational Environment and System Requirements

### 3.1 Hardware Environment
The computational analyses were performed on a high-performance computing (HPC) cluster comprising three compute nodes. One node was equipped with an Intel Xeon Gold 6130 processor (32 physical cores) and 512 GB RAM, while the other two nodes were each equipped with Intel Xeon Platinum 9242 processors (96 physical cores) and 384 GB RAM. The cluster used a 100 Gb/s Omni-Path Architecture (OPA) interconnect for communication between nodes.

The hardware specifications above describe the available computing environment rather than the minimum requirements for running MicroSVsProfiler.
### 3.2 Operating System
- Operating system: CentOS Linux release 7.9.2009 (Core)
- Architecture: x86_64
## 4. Installation and Dependencies
For detailed software requirements, dependencies, and installation instructions, see the [Installation Guide](docs/installation.md).
### 4.1 Command-line software
The workflow uses the following software:
| Software | Version |
| --- | --- |
| KneadData | 0.7.5 |
| Trimmomatic | 0.39 |
| Bowtie2 | 2.2.3 |
| seqtk | 1.3 |
| Kraken2 | 2.1.2 |
| Bracken | 2.6.2 |
| MetaPhlAn | 3.0 |
| HUMAnN | 3.0 |
| SGV-Finder | 1.0 |
| Bracken | 2.6.2 |
### 4.2 R dependencies
The R-based analyses use packages including:
- Data processing: data.table, dplyr, tidyr, plyr
- Statistical analysis: ppcor, vegan, rfPermute
- Clustering: cluster, fpc, flashClust
- Machine learning: randomForest, caret, pROC
- Visualization: ggplot2, pheatmap, cowplot, gplots
- Functional enrichment: MicrobiomeProfiler, clusterProfiler
## 5. Input Data
The workflow accepts shotgun metagenomic sequencing data and associated sample metadata.

Primary input files include:
- Paired-end metagenomic sequencing reads in FASTQ format
- Sample metadata containing sample identifiers and disease status
- Cohort information for multi-cohort analyses

Depending on the analytical step, intermediate inputs may include species abundance tables, microbial SV matrices, and functional annotation results.

For detailed input file requirements, formats, and R data import examples, see the [Input Data Formats](docs/input_formats.md).
## 6. Running MicroSVsProfiler
The repository contains stepwise scripts covering the major analytical stages, including preprocessing, microbial profiling, SV detection, differential analysis, clustering, diagnostic modeling, and functional annotation.

Users should follow the steps in their documented execution order and configure the required input paths, output directories, software locations, and database paths before execution.
## 7. Output Files
Depending on the analysis, MicroSVsProfiler generates:

Microbial species abundance matrices

Coverage-defined microbial SV feature matrices

Differential species and SV analysis results

Clustering and heterogeneity analysis outputs

Random forest diagnostic model results

Functional annotation and metabolite association results
## 8. Reproducibility, Limitations and References
**Reproducibility**

The pipeline was developed for the analyses described in the associated manuscript. Software versions, reference databases, filtering parameters, and statistical settings should be documented to support reproducibility.

Results may vary with software versions, reference database releases, sample composition, and sequencing depth.

**Limitations**

SGV-Finder identifies microbial genomic regions based on sequencing coverage patterns rather than directly resolving structural variant breakpoints. Sequence divergence, reference mismatch, species abundance, and sequencing depth may influence the resulting features.

Therefore, SV-based genomic clusters should not automatically be interpreted as phylogenetically distinct strains without independent validation.

**References**

Zeevi, D. et al., Structural variation in the gut microbiome associates with host health. *Nature*, 2019. **568**, 43-48. <https://doi.org/10.1038/s41586-019-1065-y>

Wang, D. et al., Characterization of gut microbial structural variations as determinants of human bile acid metabolism. *Cell Host Microbe*, 2021. **29**, 1802-1814.e1805. <https://doi.org/10.1016/j.chom.2021.11.003>

**Version Information**

Documentation and pipeline release: [2026-10-14].
