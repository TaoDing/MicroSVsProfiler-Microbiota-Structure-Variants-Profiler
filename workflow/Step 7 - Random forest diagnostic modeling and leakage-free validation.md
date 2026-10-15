# Step 7 - Random forest diagnostic modeling and leakage-free validation

## 1. Purpose and provenance

Two explicitly separated evaluation routes for colorectal cancer (CRC) early detection (HC versus Early): **(A) preselected-biomarker models** and **(B) leakage-free models**. It includes species abundance, dSV, vSV, and combined feature sets; stratified five-fold cross-validation and leave-one-cohort-out (LOCO) validation; random forest classifiers; held-out ROC/AUC, accuracy, sensitivity and specificity. The protocol mentions feature ranking by AUCRF/random forest and a top-20 marker set, as well as cross-cohort validation and introduces a strictly training-only modeling route to address information leakage.

> **Important distinction.** Route A use a biomarker list that was derived from all discovery samples. In that case, cross-validation performance can be optimistic and must **not** be presented as an unbiased estimate. Route B repeats marker selection and all fitted preprocessing inside each training partition. LOCO is leakage-free only if the held-out cohort was not involved in any earlier marker discovery, preprocessing fit or tuning.

## 2. Input layout

Run from the repository root. Required files are specified at the start of the embedded R code. Example files are supplied under `examples/step7/`. For actual research, replace these paths with your own real files (e.g. in `data/`).

| File | Required structure | Notes |
| --- | --- | --- |
| `metadata.tsv` | Columns `ID`, `Type2`, `Study` | `Type2` has `HC`, `Early` and optionally `Advanced`; this workflow models `HC` vs `Early` only. `Study` identifies a cohort, **not country**. |
| `species_abundance.tsv` | Rows species, columns sample IDs | Numeric abundances, preferably consistently normalized before analysis. |
| `dsv.tsv` | Rows dSV features, columns sample IDs | Binary 0/1 or documented state values; do not confuse unknown (NA) with absence (0). |
| `vsv.tsv` | Rows vSV features, columns sample IDs | Continuous normalized SV coverage values. |
| `preselected_markers.tsv` | `feature_set`, `feature` | `feature_set` is `Species`, `SV`, or `Combined`; feature IDs have prefixed names (e.g. `dSV__feature1`). |

The example code intersects sample IDs across all three feature tables, then orders all data by metadata ID. **This is a complete-case sample intersection by modality**; modify deliberately if modality-specific missing samples must be retained. Missing *feature values* are imputed using **training-only medians**, whereas completely unobserved features in a training partition are discarded. The choice to use median imputation for binary dSVs is an implementation default and should be verified against your dSV coding protocol. The code never uses 0 as a blanket NA replacement.

## 3. Evaluation designs

### 3.1 Route A: preselected-biomarker model

1. Supply a **fixed** marker list in `preselected_markers.tsv`. This corresponds to the original analysis concept of differential microbial species and SV biomarkers, followed by random forest classification.
2. Fit imputation medians and random forest **within each training fold**, but do not repeat the historical global feature screening in CV. This design is intentionally labeled *preselected* and **not leakage-free** if marker discovery used fold-test labels.
3. Evaluate three feature views: Species, SV (dSV+vSV), and Combined. If a subset has preselected markers, include the appropriate rows in the marker file.
4. An independent external cohort remains a valid assessment only when it was not consulted at any stage of marker discovery and model development.

### 3.2 Route B: leakage-free model

1. Split samples into **outer five-fold CV** or **outer LOCO** partitions.
2. For each outer training partition, run **inner stratified CV** to choose the number of markers (`5`, `10`, `20`) and random-forest `mtry` (`2`, `5`). Within each inner split, rank markers and fit imputation **only on the inner-training data**.
3. Re-rank markers using all **outer training** samples, select the tuned number, fit a new forest with 1,000 trees, and predict the held-out partition once.
4. Aggregate strictly held-out predictions for ROC/AUC, accuracy, sensitivity and specificity. Never inspect test-fold labels to select markers, model parameters, thresholds, or decision rules.
5. For LOCO, the cohort used as the test set is excluded **before any label-dependent screening**. The provided example uses `Study` as the unit of exclusion.

**Known limitations:** Nested CV can be computationally expensive; 5-fold CV requires at least five samples per class, and each inner training group must contain both classes. Multiple samples from the same participant or center require group-aware splitting. Missingness may itself be predictive; imputation and feature filtering should be prespecified and sensitivity-tested. Pooled LOCO AUC aggregates scores from models trained on different cohort combinations, so per-cohort AUC and calibration are also recommended. These examples do not claim a fully locked prospective external validation design.

## 4. Dependencies and execution

- R
- `randomForest` (version 4.7-1.2)
- `pROC` (version 1.18.5)

Install with `install.packages(c("randomForest", "pROC"))` if needed. All executable R code is embedded directly below; **no separate `scripts/` files are created**. Copy the block into your interactive R session, or paste into a temporary local R file if a batch execution is preferred. Edit paths in the `User configuration` section.

## 5. Complete R implementation

```r
# Embedded R workflow: copy the entire block to an R session from repository root.
# Packages: install.packages(c("randomForest", "pROC"))
suppressPackageStartupMessages({library(randomForest); library(pROC)})

# ---- User configuration ----
metadata_path <- "examples/step7/metadata.tsv"
species_path  <- "examples/step7/species_abundance.tsv"
dsv_path     <- "examples/step7/dsv.tsv"
vsv_path     <- "examples/step7/vsv.tsv"
preselected_path <- "examples/step7/preselected_markers.tsv"
output_root <- "results/step7"
set.seed(2026)

# Example uses a small simulated dataset; use real discovery-cohort data for research.
# Metadata columns: ID, Type2 (HC/Early/Advanced), Study. Tables: rows=features, cols=samples.
read_features <- function(path, prefix) {
  d <- read.delim(path, check.names=FALSE, row.names=1)
  if (anyDuplicated(rownames(d)) || anyDuplicated(colnames(d))) stop("Duplicate IDs: ", path)
  a <- t(as.matrix(d)); storage.mode(a) <- "double"
  colnames(a) <- paste0(prefix, "__", colnames(a))
  if (any(!is.finite(a[!is.na(a)]))) stop("Non-finite values: ", path)
  a
}
meta <- read.delim(metadata_path, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(all(c("ID", "Type2", "Study") %in% names(meta)))
if (anyDuplicated(meta$ID)) stop("Duplicate sample IDs in metadata")
blocks <- list(Species=read_features(species_path,"Species"),
               dSV=read_features(dsv_path,"dSV"),
               vSV=read_features(vsv_path,"vSV"))
ids <- Reduce(intersect, c(list(meta$ID), lapply(blocks, rownames)))
meta <- meta[match(ids,meta$ID), ,drop=FALSE]
meta <- meta[meta$Type2 %in% c("HC","Early"), ,drop=FALSE] # Early detection only
if (nrow(meta)<10) stop("At least 10 matching HC/Early samples required")
if (length(unique(meta$Study)) < 2) warning("LOCO needs >=2 cohorts")
Xall <- do.call(cbind, lapply(blocks, function(a) a[meta$ID,,drop=FALSE]))
y <- factor(meta$Type2, levels=c("HC","Early"))
cohorts <- as.character(meta$Study)
if (any(table(y)<5)) stop("Need >=5 samples in each class for 5-fold CV")
if (anyDuplicated(colnames(Xall))) stop("Feature ID collision")

feature_sets <- list(Species=grep("^Species__",colnames(Xall), value=TRUE),
                     SV=grep("^(dSV|vSV)__",colnames(Xall), value=TRUE),
                     Combined=colnames(Xall))
dir.create(output_root,recursive=TRUE,showWarnings=FALSE)

make_folds <- function(y,k=5,seed=2026) {
  if (any(table(y)<k)) stop("Insufficient sample count per class for stratified CV")
  set.seed(seed); folds <- vector("list",k)
  for (lev in levels(y)) {
    idx <- sample(which(y==lev))
    parts <- split(idx, rep(seq_len(k),length.out=length(idx)))
    for (j in seq_len(k)) folds[[j]] <- c(folds[[j]],parts[[j]])
  }
  folds
}

# Missing values: training-column median only; all-NA features are removed.
fit_preprocess <- function(x) {
  m <- apply(x,2,function(z) if(all(is.na(z))) NA_real_ else median(z,na.rm=TRUE))
  keep <- is.finite(m)
  if(!any(keep)) stop("No nonmissing feature")
  list(cols=colnames(x)[keep], med=m[keep])
}
apply_preprocess <- function(x,pp) {
  a <- x[,pp$cols,drop=FALSE]
  for(j in seq_along(pp$cols)) a[is.na(a[,j]),j] <- pp$med[j]
  a
}

# Univariate association: Wilcoxon continuous SV/abundance, Fisher for binary dSV.
# The selection is ALWAYS fitted from training samples only in the leakage-free route.
rank_markers <- function(x,y) {
  p <- vapply(seq_len(ncol(x)),function(j) {
    a <- x[y=="HC",j]; b <- x[y=="Early",j]
    if(length(unique(c(a,b)))<=1) return(1)
    if (startsWith(colnames(x)[j],"dSV__") && all(c(a,b) %in% c(0,1))) {
      tab <- table(factor(c(rep("HC",length(a)),rep("Early",length(b))),levels=c("HC","Early")),
                   factor(c(a,b),levels=c(0,1)))
      return(fisher.test(tab)$p.value)
    }
    suppressWarnings(wilcox.test(a,b,exact=FALSE)$p.value)
  },numeric(1))
  names(p) <- colnames(x)
  names(sort(p,na.last=TRUE))
}

train_rf <- function(train_x,train_y,selected,mtry,ntree=500) {
  pp <- fit_preprocess(train_x[,selected,drop=FALSE]); x <- apply_preprocess(train_x,pp)
  mtry <- max(1,min(as.integer(mtry),ncol(x)))
  fit <- randomForest(x=x,y=train_y,ntree=ntree,mtry=mtry,importance=TRUE)
  list(model=fit,pp=pp,selected=colnames(x))
}
predict_rf <- function(fit,test_x) {
  x <- apply_preprocess(test_x,fit$pp)
  as.numeric(predict(fit$model,newdata=x,type="prob")[,"Early"])
}

# Inner CV selects marker count and mtry for the OUTER TRAINING partition.
# No feature ranking, imputation or hyperparameter selection accesses outer-test rows.
tune_inner <- function(x,y,seed=17) {
  k <- min(3L,min(table(y)))
  if(k<2L) stop("Cannot tune: fewer than 2 per class in training partition")
  folds <- make_folds(y,k,seed)
  grid <- expand.grid(nmarkers=c(5L,10L,20L),mtry=c(2L,5L),KEEP.OUT.ATTRS=FALSE)
  grid$auc <- NA_real_
  for(g in seq_len(nrow(grid))) {
    pred <- rep(NA_real_,length(y))
    for(j in seq_along(folds)) {
      tr <- setdiff(seq_along(y),folds[[j]]); te <- folds[[j]]
      pp <- fit_preprocess(x[tr,,drop=FALSE]); xx <- apply_preprocess(x[tr,,drop=FALSE],pp)
      ranked <- rank_markers(xx,y[tr]); sel <- head(ranked,grid$nmarkers[g])
      fit <- train_rf(x[tr,,drop=FALSE],y[tr],sel,grid$mtry[g],ntree=300)
      pred[te] <- predict_rf(fit,x[te,,drop=FALSE])
    }
    grid$auc[g] <- as.numeric(pROC::auc(pROC::roc(y,pred,levels=c("HC","Early"),direction="<",quiet=TRUE)))
  }
  grid[order(-grid$auc,grid$nmarkers,grid$mtry),,drop=FALSE][1,,drop=FALSE]
}

# Preselected biomarkers: fixed list chosen prior to CV.
# Report as exploratory/apparent CV, NOT leakage-free unbiased validation.
markers <- read.delim(preselected_path,stringsAsFactors=FALSE)
stopifnot(all(c("feature_set","feature") %in% names(markers)))

run_analysis <- function(route,setname,design=c("CV","LOCO")) {
  design <- match.arg(design)
  cols <- feature_sets[[setname]]
  x <- Xall[,cols,drop=FALSE]
  output <- file.path(output_root,route,setname,design)
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  if (design=="CV") {
    folds <- make_folds(y,5,seed=2026)
    splits <- lapply(seq_along(folds),function(i) list(train=setdiff(seq_along(y),folds[[i]]),test=folds[[i]],label=paste0("fold",i)))
  } else {
    splits <- lapply(unique(cohorts),function(s) list(train=which(cohorts!=s),test=which(cohorts==s),label=s))
  }
  predictions <- list(); choices <- list(); importances <- list()
  for (i in seq_along(splits)) {
    z <- splits[[i]]; tr <- z$train; te <- z$test
    if(length(unique(y[tr]))!=2) {warning("Skip split with single training class: ",z$label); next}
    if(route=="preselected") {
      sel <- unique(markers$feature[markers$feature_set==setname]); sel <- intersect(sel,colnames(x))
      if(!length(sel)) stop("No preselected features for ",setname)
      config <- list(nmarkers=length(sel),mtry=max(1,floor(sqrt(length(sel)))))
    } else {
      config <- tune_inner(x[tr,,drop=FALSE],y[tr],seed=20+i)
      pp <- fit_preprocess(x[tr,,drop=FALSE]); xtr <- apply_preprocess(x[tr,,drop=FALSE],pp)
      sel <- head(rank_markers(xtr,y[tr]),config$nmarkers)
    }
    fit <- train_rf(x[tr,,drop=FALSE],y[tr],sel,config$mtry,ntree=1000)
    prob <- predict_rf(fit,x[te,,drop=FALSE]); predicted <- ifelse(prob>=0.5,"Early","HC")
    predictions[[length(predictions)+1L]] <- data.frame(ID=meta$ID[te],Study=cohorts[te],
                 TrueLabel=as.character(y[te]),Probability_Early=prob,Predicted=predicted,Split=z$label)
    choices[[length(choices)+1L]] <- data.frame(Split=z$label,n_features=length(fit$selected),
                                  mtry=fit$model$mtry,features=paste(fit$selected,collapse=";"))
    imp <- importance(fit$model,type=2)
    importances[[length(importances)+1L]] <- data.frame(Split=z$label,Feature=rownames(imp),MeanDecreaseGini=as.numeric(imp[,1]))
  }
  if(!length(predictions)) stop("No successful folds for ",route," ",setname," ",design)
  p <- do.call(rbind,predictions)
  write.table(p,file.path(output,"held_out_predictions.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  write.table(do.call(rbind,choices),file.path(output,"selected_features.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  write.table(do.call(rbind,importances),file.path(output,"feature_importance.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  # Overall pooled out-of-fold AUC; LOCO cross-cohort pooled values may need calibration.
  confusion <- table(factor(p$TrueLabel,levels=c("HC","Early")),factor(p$Predicted,levels=c("HC","Early")))
  tn <- confusion[1,1]; fp <- confusion[1,2]; fn <- confusion[2,1]; tp <- confusion[2,2]
  metric <- data.frame(route=route,features=setname,validation=design,n=nrow(p),
    AUC=if(length(unique(p$TrueLabel))==2) as.numeric(auc(roc(p$TrueLabel,p$Probability_Early,levels=c("HC","Early"),direction="<",quiet=TRUE))) else NA_real_,
    accuracy=(tp+tn)/sum(confusion),sensitivity=if(tp+fn>0) tp/(tp+fn) else NA_real_,
    specificity=if(tn+fp>0) tn/(tn+fp) else NA_real_)
  write.table(metric,file.path(output,"metrics.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  if(length(unique(p$TrueLabel))==2) {
    r <- roc(p$TrueLabel,p$Probability_Early,levels=c("HC","Early"),direction="<",quiet=TRUE)
    pdf(file.path(output,"roc.pdf"),width=6,height=6); plot(r,legacy.axes=TRUE,main=paste(route,setname,design)); dev.off()
  }
  metric
}

summary <- list()
for(route in c("preselected","leakage_free")) for(setname in names(feature_sets)) {
  for(design in c("CV","LOCO")) {
    summary[[length(summary)+1L]] <- run_analysis(route,setname,design)
  }
}
summary <- do.call(rbind,summary)
write.table(summary,file.path(output_root,"model_summary.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
print(summary)
```

## 6. Outputs

All generated results are under **`results/step7/`**, using these route-specific subfolders:

```text
results/step7/
├── model_summary.tsv
├── preselected/
│   ├── Species/{CV,LOCO}/
│   ├── SV/{CV,LOCO}/
│   └── Combined/{CV,LOCO}/
└── leakage_free/
    ├── Species/{CV,LOCO}/
    ├── SV/{CV,LOCO}/
    └── Combined/{CV,LOCO}/
```

Each successful evaluation folder contains `held_out_predictions.tsv`, `selected_features.tsv`, `feature_importance.tsv`, `metrics.tsv` and, when both labels are represented in the aggregated test data, `roc.pdf`. Feature importance uses randomForest's MeanDecreaseGini and is descriptive, not an inferential P value. AUC is computed with `HC` as the negative class and `Early` as the positive class; a fixed threshold of 0.5 is used for accuracy, sensitivity and specificity.

## 7. Reporting checklist

- Explicitly label *preselected-biomarker* CV as potentially optimistic if the marker list came from full-dataset significance tests or rankings.
- Prioritize nested CV and fully isolated LOCO results for internal or cross-cohort claims; distinguish these from independently collected external validation.
- State total samples by cohort and disease group; report marker selection protocol, missing data processing, exact training/test partition, hyperparameter grid and random seeds.
- Report confidence intervals (e.g. participant-level bootstrap on held-out predictions, respecting study structure where applicable), plus per-cohort metrics for LOCO; the reference code does not compute CIs.
- Avoid claiming that an AUC from an exploratory globally selected signature is free of information leakage.