# Step 8 - Microbial SV functional annotation and metabolite associations

This workflow integrates **microbial SV functional annotation and KEGG Orthology (KO) enrichment** with **metabolite differential testing, SV–metabolite correlation, and visualization**. All executable R code is embedded below. Run the R code from the repository root. All outputs go to **`results/step8/`**.

> **Scope:** SGV-Finder reports coverage-defined variable SVs (vSVs) and deletion-like SVs (dSVs). Functional annotations are putative gene-function associations, not experimental proof that a genomic breakpoint changes a pathway or metabolite. Correlations are observational, not causal.

## 1. Requirements and provenance

- R (version 4.2), packages: `ggplot2`, `pheatmap`, `circlize`, `ComplexHeatmap`, and `MicrobiomeProfiler` (only required for KO enrichment). Install using their appropriate CRAN/Bioconductor/GitHub channels and record `sessionInfo()`.
- The study obtained proGenomes reference genomes and used **PATRIC 3.6.6** annotations. **Do not infer KO IDs or genes from SGV-Finder feature names alone.**
- Example files below are **synthetic and not biological results**.

## 2. Inputs

| File | Format and required fields | Meaning |
|---|---|---|
| `data/metadata.tsv` | TSV: `ID`, `Type2` (`HC`, `Early`, `Advanced`), `Study` | Sample-level phenotype and cohort |
| `data/metabolites.tsv` | TSV: `ID` + metabolite numeric columns | Sample × metabolite matrix, already normalized |
| `data/sv_features.tsv` | TSV: `ID` + numeric SV columns | Sample × coverage-defined SV features; dSV 0/1 where observed |
| `data/sv_ko_annotation.tsv` | TSV: `SV_ID`, `Type`, `Knum`, `Taxon`, `Pathway` | Existing *externally derived* annotation mapping, possibly multiple KOs per SV |
| `data/significant_sv.tsv` | TSV: `SV_ID`, optional `Comparison` | SVs previously selected in Step 5; used for downstream exploratory correlations |

The provided examples live at `examples/step8/`. To run the example, set `input_dir <- "examples/step8"` in the code. **Do not use synthetic example data for biological interpretation or publication.**

Check IDs, SV names, and metabolite identifiers for uniqueness. Missing SV coverage must not be automatically treated as a dSV absence. No sample should be used unless its identifier matches the metadata.

## 3. Run: differential metabolites, SV–metabolite correlations and functional annotation

Copy the **entire block** into an interactive R session or execute it using the R console while the working directory is the repository root. Change only the configuration in the first few lines. The program stops on invalid input instead of silently producing misleading results.

```r
# ---- 3.1 Configuration / dependencies ----
input_dir <- "data"                      # or "examples/step8"
out_dir <- "results/step8"
alpha <- 0.05
min_pairs <- 4L                         # increase (e.g. 10) for real analyses
correlation_method <- "spearman"
use_step5_sv_list <- TRUE
run_enrichment <- TRUE                  # requires MicrobiomeProfiler
set.seed(2026)

required <- c("ggplot2", "pheatmap")
missing_pkg <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkg)) stop("Install R packages: ", paste(missing_pkg, collapse = ", "))
for (dir in c(out_dir, file.path(out_dir, "functional_annotation"),
              file.path(out_dir, "metabolites"), file.path(out_dir, "correlations"),
              file.path(out_dir, "figures"))) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
read_tsv <- function(name) read.delim(file.path(input_dir, name), check.names = FALSE,
                                      stringsAsFactors = FALSE)
meta <- read_tsv("metadata.tsv")
met <- read_tsv("metabolites.tsv")
sv <- read_tsv("sv_features.tsv")
ann <- read_tsv("sv_ko_annotation.tsv")
sig <- read_tsv("significant_sv.tsv")
stopifnot(all(c("ID", "Type2", "Study") %in% names(meta)),
          all(c("SV_ID", "Type", "Knum", "Taxon", "Pathway") %in% names(ann)),
          "SV_ID" %in% names(sig), "ID" %in% names(met), "ID" %in% names(sv))
for (d in list(meta, met, sv)) if (anyNA(d$ID) || anyDuplicated(d$ID))
  stop("Missing or duplicated sample ID")
if (anyDuplicated(names(sv)) || anyDuplicated(names(met)))
  stop("Duplicate feature names")
common <- Reduce(intersect, list(meta$ID, met$ID, sv$ID))
if (length(common) < 4L) stop("Fewer than four matched samples")
meta <- meta[match(common, meta$ID), , drop = FALSE]
met <- met[match(common, met$ID), , drop = FALSE]
sv <- sv[match(common, sv$ID), , drop = FALSE]
stopifnot(identical(meta$ID, met$ID), identical(meta$ID, sv$ID))
if (!all(meta$Type2 %in% c("HC", "Early", "Advanced")))
  stop("Unexpected Type2 value; expected HC, Early, Advanced")
met_ids <- setdiff(names(met), "ID")
sv_ids <- setdiff(names(sv), "ID")
for (nm in met_ids) met[[nm]] <- suppressWarnings(as.numeric(met[[nm]]))
for (nm in sv_ids) sv[[nm]] <- suppressWarnings(as.numeric(sv[[nm]]))
if (anyDuplicated(ann[c("SV_ID", "Knum")]))
  ann <- unique(ann)
write.table(data.frame(ID = common, Type2 = meta$Type2, Study = meta$Study),
            file.path(out_dir, "matched_samples.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# ---- 3.2 Differential metabolites: HC vs Early, HC vs Advanced, Early vs Advanced ----
comparisons <- list(c("HC", "Early"), c("HC", "Advanced"), c("Early", "Advanced"))
met_out <- list()
for (cmp in comparisons) {
  group_a <- cmp[1]; group_b <- cmp[2]
  a <- which(meta$Type2 == group_a); b <- which(meta$Type2 == group_b)
  if (length(a) < 2L || length(b) < 2L) {
    warning("Skipping comparison with fewer than two samples per group: ", paste(cmp, collapse = " vs "))
    next
  }
  tab <- lapply(met_ids, function(id) {
    x <- met[[id]][a]; y <- met[[id]][b]
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    p <- if (length(x) >= 2L && length(y) >= 2L &&
             (length(unique(c(x, y))) > 1L))
           suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value) else NA_real_
    data.frame(Metabolite = id, Group1 = group_a, Group2 = group_b,
               N1 = length(x), N2 = length(y), Mean1 = if(length(x)) mean(x) else NA_real_,
               Mean2 = if(length(y)) mean(y) else NA_real_, P = p)
  })
  tab <- do.call(rbind, tab)
  tab$FDR_BH <- p.adjust(tab$P, method = "BH")
  key <- paste(group_a, "vs", group_b, sep = "_")
  met_out[[key]] <- tab
  write.table(tab, file.path(out_dir, "metabolites", paste0(key, "_wilcoxon.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
if (length(met_out)) {
  write.table(do.call(rbind, met_out), file.path(out_dir, "metabolites", "all_comparisons.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

# ---- 3.3 SV–metabolite Spearman correlations ----
selected <- if (use_step5_sv_list) intersect(unique(sig$SV_ID), sv_ids) else sv_ids
if (!length(selected)) stop("No selected SVs match sv_features.tsv")
if (length(setdiff(unique(sig$SV_ID), sv_ids)))
  warning("Some significant SV IDs are missing from the SV matrix")
cor_tab <- list(); k <- 0L
for (s in selected) for (m in met_ids) {
  x <- sv[[s]]; y <- met[[m]]
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  rho <- NA_real_; p <- NA_real_
  if (n >= min_pairs && length(unique(x[ok])) > 1L && length(unique(y[ok])) > 1L) {
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = correlation_method,
                                    exact = FALSE))
    rho <- unname(ct$estimate); p <- ct$p.value
  }
  k <- k + 1L
  cor_tab[[k]] <- data.frame(SV_ID = s, Metabolite = m, N_pairs = n,
                              Spearman_rho = rho, P = p)
}
cor_tab <- do.call(rbind, cor_tab)
# BH correction within the complete tested SV × metabolite pair family.
cor_tab$FDR_BH <- p.adjust(cor_tab$P, method = "BH")
write.table(cor_tab, file.path(out_dir, "correlations", "sv_metabolite_spearman.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- 3.4 Heatmap (rho) with optional significance symbols ----
rho <- matrix(NA_real_, length(selected), length(met_ids),
              dimnames = list(selected, met_ids))
qmat <- rho
for (i in seq_len(nrow(cor_tab))) {
  row <- cor_tab[i, ]
  rho[row$SV_ID, row$Metabolite] <- row$Spearman_rho
  qmat[row$SV_ID, row$Metabolite] <- row$FDR_BH
}
if (any(is.finite(rho))) {
  # Remove entirely unavailable rows/columns; pheatmap needs a finite plotting range.
  keep_r <- rowSums(is.finite(rho)) > 0
  keep_c <- colSums(is.finite(rho)) > 0
  plot_rho <- rho[keep_r, keep_c, drop = FALSE]
  plot_q <- qmat[keep_r, keep_c, drop = FALSE]
  pdf(file.path(out_dir, "figures", "sv_metabolite_spearman_heatmap.pdf"),
      width = max(7, ncol(plot_rho) * 0.42 + 3), height = max(6, nrow(plot_rho) * 0.32 + 3))
  pheatmap::pheatmap(plot_rho, cluster_rows = FALSE, cluster_cols = FALSE,
                     display_numbers = ifelse(is.finite(plot_q) & plot_q < alpha, "*", ""),
                     color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(101),
                     breaks = seq(-1, 1, length.out = 102), na_col = "grey85",
                     main = "SV–metabolite Spearman correlation (* BH q < 0.05)")
  dev.off()
}

# ---- 3.5 Functional annotation tabulation ----
ann <- ann[!is.na(ann$Knum) & nzchar(ann$Knum) & ann$Type %in% c("VSGV", "DSGV"), , drop = FALSE]
ann$selected_in_step5 <- ann$SV_ID %in% selected
write.table(ann, file.path(out_dir, "functional_annotation", "sv_ko_mapping.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (nrow(ann)) {
  pathway <- aggregate(SV_ID ~ Type + Taxon + Pathway, data = unique(ann[c("SV_ID", "Type", "Taxon", "Pathway")]),
                       FUN = function(x) length(unique(x)))
  names(pathway)[4] <- "N_SV_features"
  write.table(pathway, file.path(out_dir, "functional_annotation", "pathway_sv_counts.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  p <- ggplot2::ggplot(pathway, ggplot2::aes(x = Taxon, y = N_SV_features, fill = Pathway)) +
    ggplot2::geom_col(position = "fill", color = "black", linewidth = 0.15) +
    ggplot2::facet_wrap(~Type) + ggplot2::coord_flip() +
    ggplot2::labs(y = "Fraction of annotated SV features", x = "Taxon") +
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(out_dir, "figures", "sv_pathway_composition.pdf"),
                  p, width = 10, height = 6)
}

# ---- 3.6 KO enrichment (MicrobiomeProfiler approach) ----
# a defensible enrichment universe and KO annotation source must be specified for the study.
if (run_enrichment) {
  if (!requireNamespace("MicrobiomeProfiler", quietly = TRUE)) {
    warning("MicrobiomeProfiler not installed; skipping KO enrichment")
  } else {
    for (type in c("VSGV", "DSGV")) {
      # Only KOs linked to the selected significant SVs; collapse duplicates.
      kos <- unique(ann$Knum[ann$Type == type & ann$selected_in_step5])
      if (!length(kos)) next
      er <- tryCatch(
        MicrobiomeProfiler::enrichKO(kos, pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH", minGSSize = 10,
                                    maxGSSize = 500, qvalueCutoff = 0.2),
        error = function(e) {warning("enrichKO failed for ", type, ": ", conditionMessage(e)); NULL})
      if (is.null(er)) next
      er_df <- tryCatch(as.data.frame(er), error = function(e) NULL)
      if (!is.null(er_df)) write.table(er_df,
          file.path(out_dir, "functional_annotation", paste0(type, "_KO_enrichment.tsv")),
          sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
```

## 4. Interpretation and statistical cautions

1. **Functional enrichment universe:** An enrichment result depends on the KO universe, reference database and the way genes are mapped to SV regions.
2. **Comparison families:** Differential-metabolite BH adjustment is applied separately within each of the three prespecified clinical contrasts. Correlation BH adjustment is across **all valid SV–metabolite pairs jointly**, not just after selecting small P values. A different declared multiple-testing family will produce different q values.
3. **Cohort effects and confounding:** For multi-cohort inference, perform cohort-stratified analyses and/or adjust confounders in a suitable model. Association does not establish a mechanistic effect.
4. **Selected features:** Using Step 5 differential SVs for exploratory correlation is a data-dependent analysis. Avoid presenting those associations as independently validated. Use a held-out cohort or prespecified external validation if confirmatory claims are needed.
5. **Missing values:** Pairwise complete cases are used for each correlation. A dSV `0` is a biologically defined absent state, not a replacement for missing sequencing coverage.
6. **Example size:** The small synthetic data can demonstrate file handling, but enrichment may yield no results under `minGSSize = 10`; correlation significance with 6 participants is not informative.

## 5. Optional visualizations

The code below gives a portable chord diagram based directly on the computed results:

```r
# Optional after running Section 3; install 'circlize' first.
if (requireNamespace("circlize", quietly = TRUE)) {
  links <- subset(cor_tab, is.finite(FDR_BH) & FDR_BH < 0.05)
  if (nrow(links)) {
    links$weight <- abs(links$Spearman_rho)
    pdf(file.path(out_dir, "figures", "significant_sv_metabolite_chord.pdf"),
        width = 12, height = 12)
    circlize::circos.clear()
    circlize::chordDiagram(links[, c("SV_ID", "Metabolite", "weight")],
                            transparency = 0.4)
    circlize::circos.clear()
    dev.off()
  } else message("No BH-significant SV–metabolite links: skipping chord diagram")
}
```

## 6. Expected output layout

```text
results/step8/
├── matched_samples.tsv
├── session_info.txt
├── metabolites/
│   ├── HC_vs_Early_wilcoxon.tsv
│   ├── HC_vs_Advanced_wilcoxon.tsv
│   ├── Early_vs_Advanced_wilcoxon.tsv
│   └── all_comparisons.tsv
├── correlations/
│   └── sv_metabolite_spearman.tsv
├── functional_annotation/
│   ├── sv_ko_mapping.tsv
│   ├── pathway_sv_counts.tsv
│   ├── VSGV_KO_enrichment.tsv        # only if available
│   └── DSGV_KO_enrichment.tsv        # only if available
└── figures/
    ├── sv_pathway_composition.pdf
    ├── sv_metabolite_spearman_heatmap.pdf
    └── significant_sv_metabolite_chord.pdf  # only when significant links exist
```
