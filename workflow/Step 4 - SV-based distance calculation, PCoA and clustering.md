# Step 4 - SV-based distance, ordination and clustering

## Overview

This step analyzes **coverage-defined microbial genomic SV features** from SGV-Finder. It combines deletion-like (dSV) and variable (vSV) profiles, calculates pairwise Canberra dissimilarities, visualizes samples with principal coordinates analysis (PCoA), tests phenotype/cohort associations using PERMANOVA, and applies **partitioning around medoids (PAM)** with silhouette-based selection of the cluster count. All outputs are written to **`results/step4/`**.

## Requirements

- R
- R packages: `vegan`, `cluster`, `ggplot2`
- Upstream SV matrices from [Step 3](Step%203%20-%20Microbial%20genomic%20SV%20detection%20using%20SGV-Finder.md)

Install missing packages in R:

```r
install.packages(c("vegan", "cluster", "ggplot2"))
```

## Input files and orientation

| Path used below | Required format | Meaning |
| --- | --- | --- |
| `data/ALL.dsgv.csv` | CSV, first column = SV feature ID, remaining columns = sample IDs | dSV feature values (features × samples) |
| `data/ALL.vsgv.csv` | CSV, same orientation | vSV feature values (features × samples) |
| `data/metadata.tsv` | Tab-delimited, columns `ID`, `Type2`, `Study` | sample phenotype and cohort labels |

The code **prefixes feature IDs** with `dSV:` or `vSV:` to avoid collisions. Samples must occur in both SV files and metadata.

`Type2` contain `HC`, `Early`, and `Advanced`; `Study` identifies cohorts. Change the variables in the first lines of the code if your actual filenames or metadata columns differ. If your SGV-Finder exports are **samples × features**, transpose them before running. If you are analyzing just a specific species, first restrict the feature matrices to that species rather than combining all species.

## Full R analysis (copy into an R session or RStudio)

The code below is self-contained. Copy it into an R/RStudio session and edit the three input paths. A small synthetic example is provided under `examples/step4/`; to run it, point `input_dsv`, `input_vsv`, and `metadata_file` to those files.

```r
# Configuration: edit the input filenames and metadata column names here.
input_dsv <- "data/ALL.dsgv.csv"  # features x samples; first column = feature ID
input_vsv <- "data/ALL.vsgv.csv"  # features x samples; first column = feature ID
metadata_file <- "data/metadata.tsv" # ID, Type2, Study (at minimum)
out_dir <- "results/step4"
min_shared_features <- 1L
max_clusters <- 6L
set.seed(100)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (pkg in c("vegan", "cluster", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install R package: ", pkg)
}

read_sv <- function(file, prefix) {
  x <- read.table(file, header = TRUE, sep = if (grepl("\\.csv$",file,ignore.case=TRUE)) "," else "\t",
                  row.names = 1, check.names = FALSE, quote = "", comment.char = "",
                  na.strings = c("NA", ""))
  if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) stop("Duplicated IDs in ", file)
  if (!nrow(x) || !ncol(x)) stop("Empty SV table: ",file)
  m <- as.matrix(x)
  suppressWarnings(storage.mode(m) <- "double")
  if (any(!is.finite(m[!is.na(m)]))) stop("Non-numeric SV entry in ", file)
  if (any(m < 0, na.rm = TRUE)) stop("Canberra distances require nonnegative values")
  rownames(m) <- paste0(prefix, ":", rownames(m))
  m
}

dsv <- read_sv(input_dsv, "dSV")
vsv <- read_sv(input_vsv, "vSV")
common <- intersect(colnames(dsv), colnames(vsv))
if (length(common) < 4) stop("At least four overlapping samples required")
sv <- t(rbind(dsv[,common,drop=FALSE], vsv[,common,drop=FALSE]))
meta <- read.delim(metadata_file, check.names=FALSE, stringsAsFactors=FALSE)
if (!all(c("ID", "Type2", "Study") %in% names(meta))) stop("Required metadata: ID, Type2, Study")
if (anyDuplicated(meta$ID)) stop("Duplicated metadata sample IDs")
matched <- rownames(sv)[rownames(sv) %in% meta$ID]
meta <- meta[match(matched, meta$ID),,drop=FALSE]
sv <- sv[matched,,drop=FALSE]
if (nrow(sv) < 4) stop("Too few samples matched to metadata")
write.table(data.frame(ID=matched, meta[,setdiff(names(meta),"ID"),drop=FALSE]),
            file.path(out_dir,"matched_metadata.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

# Pairwise-complete Canberra: NA means missing measurement, NOT SV absence.
# On a shared, nonnegative feature set, contributions are abs(a-b)/(a+b),
# with 0/0 contributing zero. Unlike zero imputation, missing values are ignored.
shared_canberra <- function(m, min_shared=1L) {
  n <- nrow(m); d <- matrix(NA_real_,n,n,dimnames=list(rownames(m),rownames(m)))
  overlap <- matrix(0L,n,n,dimnames=dimnames(d)); diag(d) <- 0
  for (i in seq_len(n)) for (j in seq_len(i)) {
    ok <- is.finite(m[i,]) & is.finite(m[j,])
    overlap[i,j] <- overlap[j,i] <- sum(ok)
    if (i == j) next
    if (sum(ok) < min_shared) next
    a <- m[i,ok]; b <- m[j,ok]; denom <- a+b
    d[i,j] <- d[j,i] <- sum(ifelse(denom > 0, abs(a-b)/denom, 0))
  }
  list(distance=d,shared=overlap)
}
ans <- shared_canberra(sv, min_shared_features)
write.table(ans$distance,file.path(out_dir,"canberra_pairwise_distance.tsv"),sep="\t",quote=FALSE)
write.table(ans$shared,file.path(out_dir,"shared_feature_counts.tsv"),sep="\t",quote=FALSE)

# A complete finite distance matrix is necessary for cmdscale/PAM/adonis2.
# Exclude samples with the most non-finite distances, documenting exclusions.
d <- ans$distance
removed <- character()
while (nrow(d) > 2 && any(!is.finite(d))) {
  bad <- rowSums(!is.finite(d))
  victim <- rownames(d)[which.max(bad)]
  removed <- c(removed,victim)
  d <- d[rownames(d)!=victim,colnames(d)!=victim,drop=FALSE]
}
writeLines(removed,file.path(out_dir,"excluded_samples_missing_distances.txt"))
if (nrow(d)<4 || any(!is.finite(d))) stop("Insufficient complete distance matrix")
meta2 <- meta[match(rownames(d),meta$ID),,drop=FALSE]
write.table(d,file.path(out_dir,"canberra_complete_distance.tsv"),sep="\t",quote=FALSE)

# PCoA and plots. Negative eigenvalues may occur with Canberra; use the positive eigenvalue total to label the displayed axes as an approximation.
pcoa <- cmdscale(as.dist(d), k=2, eig=TRUE, add=TRUE)
coords <- data.frame(ID=rownames(d), PCo1=pcoa$points[,1],PCo2=pcoa$points[,2])
coords <- merge(coords,meta2,by="ID",sort=FALSE)
coords <- coords[match(rownames(d),coords$ID),,drop=FALSE]
write.table(coords,file.path(out_dir,"pcoa_coordinates.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
pos <- pcoa$eig[pcoa$eig>0]
p1 <- if (length(pos)) 100*pcoa$eig[1]/sum(pos) else NA_real_
p2 <- if (length(pos)) 100*pcoa$eig[2]/sum(pos) else NA_real_
fig <- ggplot2::ggplot(coords,ggplot2::aes(x=PCo1,y=PCo2,color=Type2,shape=Study)) +
  ggplot2::geom_point(size=2.5,alpha=0.8) + ggplot2::theme_bw() +
  ggplot2::labs(x=sprintf("PCo1 (%.1f%% of positive eigenvalues)",p1),
                y=sprintf("PCo2 (%.1f%% of positive eigenvalues)",p2))
ggplot2::ggsave(file.path(out_dir,"pcoa_by_group_cohort.pdf"),fig,width=8,height=6)

# PERMANOVA. Include Study as a covariate where disease groups coexist within cohorts. This is association, not causal adjustment.
if (length(unique(meta2$Type2)) > 1 && length(unique(meta2$Study)) > 1) {
  perm <- vegan::adonis2(as.dist(d) ~ Study + Type2, data=meta2,
                         permutations=999, by="margin")
} else if (length(unique(meta2$Type2)) > 1) {
  perm <- vegan::adonis2(as.dist(d) ~ Type2,data=meta2,permutations=999)
} else { perm <- NULL }
if (!is.null(perm)) capture.output(perm,file=file.path(out_dir,"permanova.txt"))

# Partitioning around medoids uses a DISTANCE object, not raw rows of d.
# Select k by maximum average silhouette width.
k_values <- seq.int(2L,min(max_clusters,nrow(d)-1L))
fits <- lapply(k_values,function(k) cluster::pam(as.dist(d),k=k,diss=TRUE))
widths <- vapply(fits,function(z) z$silinfo$avg.width,numeric(1))
best <- which.max(widths); pam_fit <- fits[[best]]
write.table(data.frame(k=k_values,mean_silhouette=widths),
            file.path(out_dir,"pam_silhouette_by_k.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
clusters <- data.frame(ID=names(pam_fit$clustering),cluster=pam_fit$clustering)
clusters <- merge(clusters,meta2,by="ID",sort=FALSE)
write.table(clusters,file.path(out_dir,"pam_clusters.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
plot_data <- merge(coords,clusters[,c("ID","cluster")],by="ID")
plot_data$cluster <- factor(plot_data$cluster)
p <- ggplot2::ggplot(plot_data,ggplot2::aes(PCo1,PCo2,color=cluster,shape=Type2)) +
  ggplot2::geom_point(size=2.5) + ggplot2::theme_bw() +
  ggplot2::ggtitle(sprintf("PAM clusters, k=%d",k_values[best]))
ggplot2::ggsave(file.path(out_dir,"pam_pcoa.pdf"),p,width=8,height=6)

# Cluster-phenotype association using nominal tests. Sparse tables are assessed by Monte Carlo chi-square approximation rather than asymptotic P.
assoc <- function(field) {
  x <- clusters[!is.na(clusters[[field]]) & nzchar(as.character(clusters[[field]])),]
  tab <- table(x$cluster,x[[field]])
  if (nrow(tab)<2 || ncol(tab)<2) return(NULL)
  approx <- suppressWarnings(chisq.test(tab))
  fit <- if (any(approx$expected<5)) chisq.test(tab,simulate.p.value=TRUE,B=9999) else approx
  write.table(as.data.frame.matrix(tab),file.path(out_dir,paste0("cluster_",field,"_table.tsv")),sep="\t",quote=FALSE)
  data.frame(variable=field,p_value=fit$p.value,statistic=unname(fit$statistic),method=fit$method)
}
tests <- Filter(Negate(is.null),lapply(c("Type2","Study"),assoc))
if (length(tests)) {
  tt <- do.call(rbind,tests)
  tt$BH_adjusted_P <- p.adjust(tt$p_value,method="BH")
  write.table(tt,file.path(out_dir,"cluster_associations.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
}
cat("Step 4 complete:",nrow(d),"samples; PAM k=",k_values[best],"\n")
```

## Optional: species-specific clustering

For a species-specific analysis, apply the *same analytical procedure* to feature subsets. If feature identifiers are of the form `Species:region`, keep the relevant species before the distance calculation. For example, after reading `dsv` and `vsv` (before constructing `sv`):

```r
species <- "Blautia_wexlerae"  # example, replace with an actual feature prefix
select_species <- function(m, species) {
  ids <- sub("^[^:]+:", "", rownames(m)) # strip dSV:/vSV: prefix
  m[startsWith(ids, paste0(species, ":")), , drop = FALSE]
}
dsv <- select_species(dsv, species)
vsv <- select_species(vsv, species)
if (nrow(dsv) + nrow(vsv) == 0) stop("No matching features")
```

This selection convention must be verified against the actual feature IDs; it is not guaranteed by SGV-Finder. Run the full block separately per species, setting a distinct output folder such as `results/step4/species/Blautia_wexlerae/` to prevent overwriting.

## Outputs

All files are under `results/step4/`:

| File | Description |
| --- | --- |
| `matched_metadata.tsv` | Metadata of matched samples |
| `canberra_pairwise_distance.tsv` | Pairwise-complete Canberra distances, including missing pairs |
| `shared_feature_counts.tsv` | Number of jointly measured features per sample pair |
| `excluded_samples_missing_distances.txt` | Samples omitted to obtain a complete distance matrix |
| `canberra_complete_distance.tsv` | Finite distance matrix used downstream |
| `pcoa_coordinates.tsv` / `pcoa_by_group_cohort.pdf` | Ordination and visualization |
| `permanova.txt` | PERMANOVA results when eligible |
| `pam_silhouette_by_k.tsv` | Average silhouette width for each candidate k |
| `pam_clusters.tsv` / `pam_pcoa.pdf` | PAM group assignments and visualization |
| `cluster_associations.tsv` | Chi-square associations and BH-adjusted P values |

## Methods and interpretation notes

1. **Distances:** For two samples, only jointly observed SV values contribute to their Canberra distance. All-zero shared features contribute zero. Pairwise deletion creates distances computed from different feature sets; this can affect geometry and comparability. Inspect `shared_feature_counts.tsv` and consider a shared-feature threshold before interpreting clusters.
2. **PCoA:** Canberra dissimilarities can be non-Euclidean; the code requests additive correction (`cmdscale(add=TRUE)`). Axis labels indicate the share of positive eigenvalues and are descriptive.
3. **PAM:** The `cluster::pam()` call explicitly uses `diss=TRUE`.
4. **Association tests:** Chi-square/Monte Carlo P values are nominal unless corrected. The small association summary applies BH correction across tests in that file. For species-wise screens, correct across the complete species testing family instead.
5. **SV-based clusters are not confirmed strains.** Distinct phylogenetic strains require independent genomic evidence (e.g., SNVs, assemblies, or phylogenies).