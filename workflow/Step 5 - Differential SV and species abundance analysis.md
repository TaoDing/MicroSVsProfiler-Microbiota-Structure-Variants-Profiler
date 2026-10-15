# Step 5 - Differential SV and species abundance analysis

## 1. Purpose and scope

This workflow compares healthy controls (HC) against early-stage colorectal cancer (Early) and advanced-stage colorectal cancer (Advanced) **within each cohort** and produces:

- vSV (variable genomic region) differential abundance results using a two-sided Wilcoxon rank-sum test;
- dSV (deletion-like genomic region) presence/absence results using chi-square test;
- cohort-specific, species-level relative abundance tables and classic LEfSe inputs;
- nominal P values and Benjamini–Hochberg (BH) adjusted P values across the **complete tested feature set** for each analysis;
- optional cross-cohort UpSet plots of differential dSV and vSV features.

**Terminology:** dSVs and vSVs are *coverage-defined microbial genomic variation features* inferred by SGV-Finder, not breakpoint-resolved structural rearrangements. Tabular processing is performed in **R**, and LEfSe is called from **Shell**.

## 2. Prerequisites

Run from the repository root. Expected dependencies:

| Component | Purpose |
| --- | --- |
| R ≥4.0 | Tabular processing, statistical tests, LEfSe input generation |
| Base R (`stats`, `utils`) | Required statistical analyses; no extra R packages for the core workflow |
| `UpSetR` (optional) | Cross-cohort intersections of significant SV features |
| LEfSe (`lefse_format_input.py`, `lefse_run.py`, `lefse_plot_res.py`) | LDA-based species-level biomarker analysis |

LEfSe requires a compatible Python environment. Verify the installed CLI help and resolve NumPy/h5py ABI conflicts in that environment rather than altering the system's base environment.

## 3. Input files and formats

Place input files under `data/` (or edit the paths in the R code below):

| File | Format and orientation | Key fields |
| --- | --- | --- |
| `data/metadata.tsv` | Tab-delimited; **one sample per row** | `ID`, `Study`, `Type2` |
| `data/ALL.vsgv.csv` | Comma-delimited; **one sample per row**, vSVs in columns | `ID`, feature identifiers, numeric coverage-derived values; `NA` allowed |
| `data/ALL.dsgv.csv` | Comma-delimited; **one sample per row**, dSVs in columns | `ID`, binary 0/1 calls; `NA` means missing/uncallable |
| `data/species_reads.tsv` | Tab-delimited; **one sample per row**, species in columns | `ID`, non-negative read counts |


The example files in [`examples/step5/`](../examples/step5/) provide the **format only**; they contain synthetic values and are not biological data.

> **Important:** Run the script only after verifying that the column `Type2` uses the categories `HC`, `Early`, and `Advanced`, and that the `Study` identifiers match the intended cohorts. Do not silently map unknown categories to CRC. If your data use a different orientation, transpose deliberately before running, then confirm that row IDs are sample IDs.

## 4. Output directory

**All output files are written under `results/step5/`.**

```text
results/step5/
├── sample_match_summary.tsv
├── missing_samples_*.tsv
├── differential_sv/
│   ├── <cohort>_HC_vs_Early_vSV.tsv
│   ├── <cohort>_HC_vs_Early_dSV.tsv
│   ├── <cohort>_HC_vs_Advanced_vSV.tsv
│   └── <cohort>_HC_vs_Advanced_dSV.tsv
├── species/
│   ├── <cohort>_HC_vs_Early_species_statistics.tsv
│   ├── <cohort>_HC_vs_Early_LEfSe.tsv
│   └── ...
├── lefse/
│   ├── <cohort>_HC_vs_Early.in
│   ├── <cohort>_HC_vs_Early.res
│   └── <cohort>_HC_vs_Early_LDA.pdf
└── intersections/
    ├── dSV_membership.tsv
    ├── vSV_membership.tsv
    └── *_upset.pdf (optional)
```

Files for comparisons with no matched subjects are **not generated**. A missing group is reported rather than silently filled with fabricated samples.

## 5. Core R workflow (run the entire block)

Copy this block into an R session **from the repository root**. It is self-contained and does not require saving a separate `.R` file. It matches all data by sample `ID` and performs cohort-wise tests.

```r
# ---- Step 5 configuration -------------------------------------------------
output_root <- "results/step5"
metadata_file <- "data/metadata.tsv"
vsv_file <- "data/ALL.vsgv.csv"
dsv_file <- "data/ALL.dsgv.csv"
species_file <- "data/species_reads.tsv"

# Set to NULL for all cohorts in metadata, or e.g. c("01HK", "02AUS").
selected_cohorts <- NULL
comparisons <- list(c("HC", "Early"), c("HC", "Advanced"))

for (subdir in c("", "differential_sv", "species", "lefse", "intersections")) {
  dir.create(file.path(output_root, subdir), recursive = TRUE,
             showWarnings = FALSE)
}

# ---- Import and validation ------------------------------------------------
read_id_table <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  delimiter <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  x <- read.table(path, header = TRUE, sep = delimiter, quote = '"',
                  comment.char = "", check.names = FALSE, fill = FALSE,
                  stringsAsFactors = FALSE)
  if (!"ID" %in% names(x)) stop(path, " must include a column named ID")
  x$ID <- trimws(as.character(x$ID))
  if (anyNA(x$ID) || any(!nzchar(x$ID)) || anyDuplicated(x$ID))
    stop("Missing, empty, or duplicated ID in: ", path)
  x
}

read_numeric_table <- function(path) {
  x <- read_id_table(path)
  if (ncol(x) < 2) stop("No features in ", path)
  ids <- x$ID
  raw <- x[, -1, drop = FALSE]
  if (anyDuplicated(colnames(raw))) stop("Duplicate features in ", path)
  numeric_x <- lapply(raw, function(z) suppressWarnings(as.numeric(z)))
  for (j in seq_along(raw)) {
    before <- !is.na(raw[[j]]) & nzchar(trimws(as.character(raw[[j]])))
    if (any(before & is.na(numeric_x[[j]])))
      stop("Non-numeric values in feature: ", names(raw)[j])
  }
  ans <- as.matrix(as.data.frame(numeric_x, check.names = FALSE))
  rownames(ans) <- ids
  colnames(ans) <- names(raw)
  ans
}

meta <- read_id_table(metadata_file)
if (!all(c("Study", "Type2") %in% names(meta)))
  stop("Metadata requires Study and Type2 columns")
meta$Study <- as.character(meta$Study)
meta$Type2 <- as.character(meta$Type2)
if (anyNA(meta$Study) || anyNA(meta$Type2))
  stop("Missing Study or Type2 metadata; resolve before analysis")

vsv <- read_numeric_table(vsv_file)
dsv <- read_numeric_table(dsv_file)
species_counts <- read_numeric_table(species_file)
if (any(!is.na(dsv) & !dsv %in% c(0, 1)))
  stop("dSV data must use 0/1/NA; do not replace NA with 0")
if (any(species_counts < 0, na.rm = TRUE) || anyNA(species_counts))
  stop("Species reads must be nonnegative and complete")
if (any(rowSums(species_counts) == 0))
  stop("At least one sample has zero total species-assigned reads")

# Abundance (%) is based on ALL original species, before any feature filtering.
species_pct <- sweep(species_counts, 1, rowSums(species_counts), "/") * 100

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}

# Explicit sample matching; unmatched samples are reported and excluded.
cohorts <- if (is.null(selected_cohorts)) unique(meta$Study) else selected_cohorts
cohorts <- cohorts[!is.na(cohorts) & nzchar(cohorts)]
summary_rows <- list()
result_membership <- list(dSV = list(), vSV = list())

# For each sample and its partners, non-callable features remain NA.
wilcox_feature <- function(x, g, control, case) {
  x1 <- x[g == control]; x2 <- x[g == case]
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  p <- if (length(x1) && length(x2)) {
    if (length(unique(c(x1, x2))) < 2) 1 else
      tryCatch(wilcox.test(x1, x2, exact = FALSE)$p.value,
               error = function(e) NA_real_)
  } else NA_real_
  c(P_value = p,
    HC_mean = if (length(x1)) mean(x1) else NA_real_,
    Case_mean = if (length(x2)) mean(x2) else NA_real_,
    HC_n_observed = length(x1), Case_n_observed = length(x2))
}

dsv_feature <- function(x, g, control, case) {
  x1 <- x[g == control]; x2 <- x[g == case]
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  n1 <- length(x1); n2 <- length(x2)
  if (!n1 || !n2) {
    return(data.frame(P_value = NA_real_, test = "not_testable",
      HC_present = sum(x1 == 1), HC_absent = sum(x1 == 0),
      Case_present = sum(x2 == 1), Case_absent = sum(x2 == 0),
      HC_n_observed = n1, Case_n_observed = n2))
  }
  tab <- matrix(c(sum(x1 == 1), sum(x1 == 0),
                  sum(x2 == 1), sum(x2 == 0)), nrow = 2, byrow = TRUE)
  if (length(unique(c(x1, x2))) < 2) {
    p <- 1; method <- "no_variation"
  } else {
    expected <- suppressWarnings(chisq.test(tab, correct = FALSE)$expected)
    if (any(expected < 5)) {
      p <- fisher.test(tab)$p.value; method <- "Fisher_exact"
    } else {
      p <- chisq.test(tab, correct = FALSE)$p.value
      method <- "Pearson_chisq"
    }
  }
  data.frame(P_value = p, test = method,
    HC_present = tab[1, 1], HC_absent = tab[1, 2],
    Case_present = tab[2, 1], Case_absent = tab[2, 2],
    HC_n_observed = n1, Case_n_observed = n2)
}

for (cohort in cohorts) {
  for (contrast in comparisons) {
    control <- contrast[1]; case <- contrast[2]
    comparison <- paste0("HC_vs_", case)
    tag <- paste(cohort, comparison, sep = "_")
    subset_meta <- meta[meta$Study == cohort & meta$Type2 %in% contrast,
                        , drop = FALSE]
    # Only samples available in all required matrices are analyzed.
    common <- Reduce(intersect, list(subset_meta$ID,
      rownames(vsv), rownames(dsv), rownames(species_pct)))
    matched_meta <- subset_meta[match(common, subset_meta$ID), , drop = FALSE]
    counts <- table(factor(matched_meta$Type2, levels = contrast))
    summary_rows[[tag]] <- data.frame(cohort = cohort, comparison = comparison,
      n_metadata = nrow(subset_meta), n_matched = length(common),
      n_HC = as.integer(counts[1]), n_case = as.integer(counts[2]))

    missing <- subset_meta[!subset_meta$ID %in% common,
                            c("ID", "Study", "Type2"), drop = FALSE]
    if (nrow(missing)) write_tsv(missing,
      file.path(output_root, paste0("missing_samples_", tag, ".tsv")))

    if (any(counts == 0)) {
      message("Skipping ", tag, ": missing HC or ", case, " samples")
      next
    }
    groups <- matched_meta$Type2

    # ---- Continuous vSV features -----------------------------------------
    xv <- vsv[common, , drop = FALSE]
    vstat <- t(vapply(seq_len(ncol(xv)), function(j)
      wilcox_feature(xv[, j], groups, control, case), numeric(5)))
    vres <- data.frame(Feature = colnames(xv), vstat, check.names = FALSE)
    names(vres)[names(vres) == "HC_mean"] <- "HC_mean_coverage"
    names(vres)[names(vres) == "Case_mean"] <- "Case_mean_coverage"
    vres$BH_adjusted_P <- p.adjust(vres$P_value, "BH")
    vres$Comparison <- comparison; vres$Study <- cohort
    write_tsv(vres, file.path(output_root, "differential_sv",
                               paste0(tag, "_vSV.tsv")))
    result_membership$vSV[[tag]] <- vres$Feature[
      !is.na(vres$BH_adjusted_P) & vres$BH_adjusted_P < 0.05]

    # ---- Binary dSV features ---------------------------------------------
    xd <- dsv[common, , drop = FALSE]
    dres <- do.call(rbind, lapply(seq_len(ncol(xd)), function(j)
      cbind(Feature = colnames(xd)[j],
            dsv_feature(xd[, j], groups, control, case))))
    dres$P_value <- as.numeric(dres$P_value)
    dres$BH_adjusted_P <- p.adjust(dres$P_value, "BH")
    dres$HC_prevalence <- with(dres,
       HC_present / (HC_present + HC_absent))
    dres$Case_prevalence <- with(dres,
       Case_present / (Case_present + Case_absent))
    dres$Comparison <- comparison; dres$Study <- cohort
    write_tsv(dres, file.path(output_root, "differential_sv",
                               paste0(tag, "_dSV.tsv")))
    result_membership$dSV[[tag]] <- dres$Feature[
      !is.na(dres$BH_adjusted_P) & dres$BH_adjusted_P < 0.05]

    # ---- Species abundance statistics ------------------------------------
    xs <- species_pct[common, , drop = FALSE]
    s_p <- apply(xs, 2, function(x) {
      if (length(unique(x)) < 2) return(1)
      kruskal.test(x ~ factor(groups, levels = contrast))$p.value
    })
    sres <- data.frame(Species = colnames(xs),
      HC_mean_pct = colMeans(xs[groups == control, , drop = FALSE]),
      Case_mean_pct = colMeans(xs[groups == case, , drop = FALSE]),
      P_value = as.numeric(s_p),
      BH_adjusted_P = p.adjust(as.numeric(s_p), "BH"),
      Comparison = comparison, Study = cohort, row.names = NULL,
      check.names = FALSE)
    write_tsv(sres, file.path(output_root, "species",
                 paste0(tag, "_species_statistics.tsv")))

    # Classic LEfSe input: first row sample IDs, second row is Class.
    # In the formatter, -c 2 identifies Class; no subclass/subject is set.
    lefse_path <- file.path(output_root, "species", paste0(tag, "_LEfSe.tsv"))
    con <- file(lefse_path, open = "wt")
    writeLines(paste(c("SampleID", common), collapse = "\t"), con)
    writeLines(paste(c("Class", groups), collapse = "\t"), con)
    # Drop features that are all zero in this comparison, for LEfSe input only.
    eligible <- colnames(xs)[colSums(xs > 0) > 0]
    for (feature in eligible) {
      writeLines(paste(c(feature,
        format(xs[, feature], scientific = FALSE, trim = TRUE)),
        collapse = "\t"), con)
    }
    close(con)
    message("Wrote results for ", tag)
  }
}
write_tsv(do.call(rbind, summary_rows),
          file.path(output_root, "sample_match_summary.tsv"))

# ---- Optional cross-cohort feature intersection matrices ------------------
for (kind in c("dSV", "vSV")) {
  sets <- result_membership[[kind]]
  sets <- sets[lengths(sets) > 0]
  if (!length(sets)) next
  all_features <- sort(unique(unlist(sets, use.names = FALSE)))
  if (!length(all_features)) next
  membership <- data.frame(Feature = all_features, check.names = FALSE)
  for (nm in names(sets)) membership[[nm]] <- as.integer(all_features %in% sets[[nm]])
  write_tsv(membership, file.path(output_root, "intersections",
                                  paste0(kind, "_membership.tsv")))
  if (requireNamespace("UpSetR", quietly = TRUE) && ncol(membership) >= 3) {
    pdf(file.path(output_root, "intersections", paste0(kind, "_upset.pdf")),
        width = 11, height = 6)
    UpSetR::upset(membership[, -1, drop = FALSE],
                 nsets = ncol(membership) - 1, order.by = "freq")
    dev.off()
  }
}
message("Finished. See ", output_root)
```

**Statistical policy.** Two-sided Wilcoxon tests are used for vSV values; dSV calls are treated as binary presence/absence. For dSVs, Fisher's exact test is used when any chi-squared expected cell count is <5; otherwise, Pearson's chi-squared test is used. BH correction is applied separately within each cohort, contrast, and feature family. All tested features are included, including constant features assigned P=1; features with no observations in either group are not testable and receive `NA`. Do **not** interpret uncallable dSVs (`NA`) as absent (`0`). These rules provide a documented, consistent implementation.

## 6. Species-level LEfSe analysis

Run the preceding R block first. The generated LEfSe input is a transposed taxon-by-sample table, with **row 1 containing sample IDs, row 2 containing `Class`, and subsequent rows containing relative abundance (%)**.

For example, from the repository root:

```bash
# Example: one cohort, HC versus Early
prefix="01HK_HC_vs_Early"
lefse_format_input.py "results/step5/species/${prefix}_LEfSe.tsv" \
  "results/step5/lefse/${prefix}.in" -c 2 -o 1000000
lefse_run.py "results/step5/lefse/${prefix}.in" \
  "results/step5/lefse/${prefix}.res" -l 2.0
lefse_plot_res.py "results/step5/lefse/${prefix}.res" \
  "results/step5/lefse/${prefix}_LDA.pdf" --format pdf
```

The traditional LEfSe result (`.res`) often contains taxon name, mean abundance, enriched group, LDA score, and a nominal P value; check the actual output columns before merging. LDA scores are distinct from effect-size estimates such as the difference in group means and generally have no default confidence interval.

### Optional: merge BH-adjusted species P values with LEfSe output

The R analysis already writes `*_species_statistics.tsv` calculated using the **full count-derived species matrix**, not a preselected LEfSe hit list. To attach LEfSe results by species ID:

```r
prefix <- "01HK_HC_vs_Early"
base <- "results/step5"
stats <- read.delim(file.path(base, "species", paste0(prefix, "_species_statistics.tsv")),
                    check.names = FALSE)
lefse <- read.delim(file.path(base, "lefse", paste0(prefix, ".res")),
                    header = FALSE, fill = TRUE, check.names = FALSE,
                    na.strings = c("", "NA"))
if (ncol(lefse) != 5L) stop("Check the actual LEfSe .res layout first")
names(lefse) <- c("Species", "LEfSe_mean", "Enriched_group", "LDA_score", "LEfSe_P")
if (anyDuplicated(lefse$Species)) stop("Duplicated species names in LEfSe output")
idx <- match(stats$Species, lefse$Species)
stats$LDA_score <- as.numeric(lefse$LDA_score[idx])
stats$Enriched_group <- lefse$Enriched_group[idx]
stats$LEfSe_P <- as.numeric(lefse$LEfSe_P[idx])
write.table(stats, file.path(base, "species", paste0(prefix, "_LEfSe_with_BH.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
```

**Reporting caution:** The independently calculated `P_value` / `BH_adjusted_P` values are from a full-matrix Kruskal–Wallis analysis; `LEfSe_P` and `LDA_score` come from the LEfSe run. They are different statistical outputs and should retain separate column labels in supplementary tables. Nominal LEfSe significance does not imply BH-adjusted significance. Analyses spanning cohorts should account for study effects rather than interpreting pooled, unadjusted LEfSe as cohort-adjusted inference.

## 7. Interpretation and checks

1. Review `sample_match_summary.tsv` and missing-sample reports before interpreting any test.
2. vSV values must represent comparable SGV-Finder coverage-derived measurements; do not replace missing values with zero before testing.
3. dSV inputs must be genuinely binary: `1` = detected/present, `0` = absent with adequate coverage, `NA` = not callable. If your exported dSV matrix uses a different encoding, recode it explicitly after verifying SGV-Finder semantics.
4. Confirm that Kraken2/Bracken output counts are species-level read estimates and that species names are unique and stable across matrices. Report whether classification-based relative abundance is normalized over all species-assigned reads or a broader read universe.
5. Review independence and possible confounding by study, age, sex, and sequencing depth. Unadjusted univariate tests do not correct for clinical covariates.
6. For UpSet plots, the sample sets are defined by **BH-adjusted P < 0.05** within each cohort/contrast.