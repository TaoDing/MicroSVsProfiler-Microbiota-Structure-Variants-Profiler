# Step 6 - Hierarchical clustering and cohort heterogeneity analysis

## Overview

This step combines two analyses:

1. **Hierarchical clustering of differential SV features**, to compare healthy controls (HC) with early-stage CRC and with advanced-stage CRC **within each cohort**.
2. **Rank-based ANOVA-type heterogeneity analysis**, to quantify how much variation in each microbial SV feature is individually associated with disease status and cohort membership.

All executable R code is included below; no separate file under `scripts/` is required.

## Requirements

- R
- CRAN packages: `pheatmap`, `ggplot2`
- Completed SV profiles from Step 3 and differential SV feature IDs from Step 5
- Metadata with matching sample IDs and the fields `Type2` and `Study`

Install once in R:

```r
install.packages(c("pheatmap", "ggplot2"))
```

## Input files

The following are self-contained **synthetic examples** (not patient data):

| File | Purpose | Required orientation |
| --- | --- | --- |
| `examples/step6/sv_matrix.example.tsv` | Combined dSV/vSV coverage or presence profile | Rows = sample IDs; columns = SV features |
| `examples/step6/metadata.example.tsv` | Matched clinical grouping and cohort | Columns = `ID`, `Type2`, `Study` |
| `examples/step6/differential_sv_features.example.tsv` | Preselected differential SV feature list from Step 5 | Column = `Feature` |

For real analyses, replace these paths in the configuration block below. The first column of the SV matrix must be `ID`; SV identifiers must exactly match those in the differential feature list. Do not confuse an SGV-Finder SV feature identifier with a precise genomic coordinate.

## Output directory

All outputs are written to **`results/step6/`**. The code creates that folder automatically. Important files include:

| Output | Description |
| --- | --- |
| `*_differential_SV_heatmap.pdf` | Within-cohort heatmap for HC vs Early or HC vs Advanced |
| `*_feature_clusters.tsv` | Hierarchical row-cluster assignment for differential SV features |
| `heatmap_summary.tsv` | Samples and retained SV features in each heatmap |
| `rank_based_heterogeneity.tsv` | Marginal disease- and cohort-associated rank R² for each eligible SV feature |
| `heterogeneity_disease_vs_cohort.pdf` | Disease vs cohort heterogeneity scatterplot |
| `matched_samples.tsv` | Final matched sample IDs and annotations |

## Hierarchical clustering of differential SV features

Each cohort is analysed separately, using HC vs Early and HC vs Advanced comparisons. The feature list must originate from **Step 5**, not from an arbitrary selection of all SVs. The heatmap depicts per-feature row-wise standardized values for visualization only. Hierarchical clustering uses the default `pheatmap` distance/linkage settings (Euclidean distance and complete linkage). A maximum of seven row clusters is requested; for small demonstration matrices the number is capped at the available SV features.

For heatmap clustering, features with any missing observation in the selected samples, or without variability, are omitted. This avoids silently converting missing coverage to absence. With real data, document and justify alternative missing-data filtering if used.

## Rank-based ANOVA-type heterogeneity analysis

The method ranks each feature across samples and reports `1 - (within-group sum of squares / total sum of squares)` for **disease status** and **cohort** separately. This corresponds to a *marginal, one-factor explained-variance fraction* on rank-transformed data, rather than a joint two-factor ANOVA partition or covariate-adjusted model. The two columns do **not** add up to 100%, and observed associations can reflect cohort–disease imbalance.

## Complete R workflow

Run the following R block from the **repository root** using RStudio or by pasting it into an R session.

```r
# Run from the MicroSVsProfiler repository root.
# Install once: install.packages(c("pheatmap", "ggplot2"))
suppressPackageStartupMessages({
  library(pheatmap)
  library(ggplot2)
})

# ---------- Configuration ----------
input_matrix <- "examples/step6/sv_matrix.example.tsv"
metadata_file <- "examples/step6/metadata.example.tsv"
differential_features_file <- "examples/step6/differential_sv_features.example.tsv"
out_dir <- "results/step6"
group_col <- "Type2"                  # HC / Early / Advanced
cohort_col <- "Study"                  # independent cohorts
comparisons <- c("Early", "Advanced") # each compared with HC
min_samples_per_group <- 2
heatmap_row_clusters <- 7
min_nonzero_for_heterogeneity <- 2 # adjust for real data
set.seed(2024)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE,
             stringsAsFactors = FALSE, quote = "", comment.char = "")
}

# Matrix orientation REQUIRED: samples as rows; SV features as columns.
sv_raw <- read_tsv(input_matrix)
meta <- read_tsv(metadata_file)
selected <- read_tsv(differential_features_file)
stopifnot("ID" %in% names(sv_raw), "ID" %in% names(meta),
          group_col %in% names(meta), cohort_col %in% names(meta),
          "Feature" %in% names(selected))
if (anyDuplicated(sv_raw$ID) || anyDuplicated(meta$ID)) {
  stop("Sample IDs must be unique in both input files")
}
feature_cols <- setdiff(names(sv_raw), "ID")
if (anyDuplicated(feature_cols)) stop("Duplicate SV feature IDs")
sv_mat <- as.matrix(sv_raw[, feature_cols, drop = FALSE])
suppressWarnings(storage.mode(sv_mat) <- "numeric")
if (any(!is.finite(sv_mat[!is.na(sv_mat)]))) stop("Non-finite SV values")
rownames(sv_mat) <- sv_raw$ID
common <- meta$ID[meta$ID %in% rownames(sv_mat)]
if (!length(common)) stop("No matching sample IDs")
meta <- meta[match(common, meta$ID), , drop = FALSE]
sv_mat <- sv_mat[common, , drop = FALSE]
stopifnot(identical(meta$ID, rownames(sv_mat)))
if (anyNA(meta[[group_col]]) || anyNA(meta[[cohort_col]])) {
  stop("Missing disease/cohort values among matched samples")
}
selected_ids <- unique(as.character(selected$Feature))
selected_ids <- intersect(selected_ids, colnames(sv_mat))
if (length(selected_ids) < 2) stop("Fewer than two selected SV features matched")
cat("Matched samples:", nrow(sv_mat), "Selected SVs:",
    length(selected_ids), "\n")
write.table(data.frame(ID=common, Group=meta[[group_col]], Study=meta[[cohort_col]]),
            file.path(out_dir, "matched_samples.tsv"), sep="\t",
            quote=FALSE, row.names=FALSE)

# ---------- Part A: hierarchical clustering of differential SV features ----------
# Run separately by cohort and case-control comparison.
# All feature selection must have been completed beforehand (Step 5).
make_heatmap <- function(cohort, case_group) {
  ix <- which(meta[[cohort_col]] == cohort &
              meta[[group_col]] %in% c("HC", case_group))
  if (length(ix) == 0) return(NULL)
  m <- meta[ix, , drop=FALSE]
  if (any(table(factor(m[[group_col]], levels=c("HC",case_group))) <
          min_samples_per_group)) return(NULL)
  x <- sv_mat[m$ID, selected_ids, drop=FALSE]
  # Missing coverage must not be silently coded as 0.
  # For clustering, retain features fully observed in these samples and varying.
  good <- colSums(is.na(x)) == 0 &
          apply(x, 2, function(z) length(unique(z[!is.na(z)])) >= 2)
  x <- x[, good, drop=FALSE]
  if (ncol(x) < 2) return(NULL)
  h <- t(x)                     # pheatmap: SV features x samples
  # Row-wise z-score is for visualization, NOT the input for statistical tests.
  h <- t(scale(t(h)))
  ann <- data.frame(Group=factor(m[[group_col]], levels=c("HC",case_group)))
  rownames(ann) <- m$ID
  tag <- paste(cohort, paste0("HC_vs_", case_group), sep="_")
  k <- min(heatmap_row_clusters, nrow(h))
  obj <- pheatmap(h, annotation_col=ann, cluster_rows=TRUE,
                  cluster_cols=TRUE, show_colnames=FALSE,
                  show_rownames=nrow(h) <= 80, border_color=NA,
                  cutree_rows=if (k >= 2) k else NA,
                  silent=TRUE)
  pdf(file.path(out_dir, paste0(tag, "_differential_SV_heatmap.pdf")),
      width=10, height=max(6, min(22, nrow(h)*0.20+3)))
  grid::grid.newpage(); grid::grid.draw(obj$gtable); dev.off()
  rowcluster <- if (k >= 2) cutree(obj$tree_row, k=k) else
    setNames(1L, rownames(h))
  write.table(data.frame(Feature=names(rowcluster), Cluster=as.integer(rowcluster)),
              file.path(out_dir, paste0(tag,"_feature_clusters.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
  data.frame(Cohort=cohort, Comparison=paste("HC",case_group,sep="_vs_"),
             Samples=nrow(x), Features=ncol(x), Row_clusters=k)
}
heatmap_qc <- list()
for (cohort in unique(meta[[cohort_col]])) {
  for (case_group in comparisons) {
    result <- make_heatmap(cohort, case_group)
    if (!is.null(result)) heatmap_qc[[length(heatmap_qc)+1]] <- result
  }
}
if (length(heatmap_qc)) {
  write.table(do.call(rbind,heatmap_qc),
              file.path(out_dir,"heatmap_summary.tsv"), sep="\t",
              quote=FALSE,row.names=FALSE)
} else warning("No cohort/comparison had sufficient samples for heatmaps")

# ---------- Part B: rank-based, marginal ANOVA-type heterogeneity ----------
# Calculates marginal explained variance for Study and Type2 SEPARATELY.
# These values are not adjusted for each other and do not add to 100%.
rank_r2 <- function(values, group) {
  keep <- !is.na(values) & !is.na(group)
  x <- values[keep]; g <- factor(group[keep])
  if (length(x) < 3 || nlevels(g) < 2) return(NA_real_)
  r <- rank(x, ties.method="average") / length(x)
  total_ss <- sum((r-mean(r))^2)
  if (total_ss == 0) return(NA_real_)
  within_ss <- sum(vapply(split(r,g), function(z) sum((z-mean(z))^2),
                          numeric(1)))
  1-within_ss/total_ss
}
nonzero_count <- colSums(!is.na(sv_mat) & sv_mat != 0)
keep_features <- colnames(sv_mat)[nonzero_count > min_nonzero_for_heterogeneity]
if (!length(keep_features)) stop("No features meet the heterogeneity prevalence cutoff")
hetero <- do.call(rbind,lapply(keep_features, function(feature) {
  x <- sv_mat[,feature]
  data.frame(Feature=feature,
             Study_R2=rank_r2(x,meta[[cohort_col]]),
             Disease_R2=rank_r2(x,meta[[group_col]]),
             Nonzero_samples=sum(!is.na(x) & x != 0),
             Nonmissing_samples=sum(!is.na(x)))
}))
write.table(hetero, file.path(out_dir,"rank_based_heterogeneity.tsv"),
            sep="\t", quote=FALSE,row.names=FALSE,na="NA")
plot_df <- hetero[is.finite(hetero$Study_R2) &
                  is.finite(hetero$Disease_R2),,drop=FALSE]
if (nrow(plot_df)) {
  p <- ggplot(plot_df,aes(x=Disease_R2,y=Study_R2,
                         size=Nonzero_samples)) +
    geom_point(alpha=0.65) +
    geom_abline(slope=1,intercept=0,linetype="dashed") +
    coord_equal() + theme_bw(base_size=12) +
    labs(x="Disease status explained variance (marginal rank R2)",
         y="Cohort explained variance (marginal rank R2)",
         size="Nonzero samples")
  ggsave(file.path(out_dir,"heterogeneity_disease_vs_cohort.pdf"),
         p,width=7,height=6)
}
cat("Step 6 completed; outputs in",out_dir,"\n")
```

## Interpretation and reproducibility notes

- **SV-based feature clusters are not validated microbial strains.** They describe groups of features with similar patterns, not phylogenetically defined isolates.
- The clustering uses differential SVs; therefore, the observed grouping is **descriptive**, and the same feature-selection results should not be treated as independent validation.
- The disease and cohort R² values represent **separate, unadjusted groupings**, not mutually independent causal contributions. To estimate cohort-adjusted disease effects, fit an appropriately specified multivariable model and report it separately.
- `Study` indicates the cohort and `Type2` indicates HC, Early, or Advanced. If your final data use another field, update `group_col` and `cohort_col` explicitly.