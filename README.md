
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dCAM

<!-- badges: start -->

<!-- badges: end -->

We present reference-free differential Convex Analysis of Mixtures
(dCAM), a joint deconvolution method of analyzing grouped bulk data of
varying composition that can improve detection of cell type specific
differential expressions in many biological contexts. With an L1-norm
differential expression penalty-based regularization, dCAM jointly and
iteratively estimates both cell type proportions and cell type specific
gene expressions in grouped bulk data, and performs cell type-specific
differential expression analysis between groups. Through realistic
grouped bulk simulations, we demonstrate that dCAM can better detect
differentially expressed genes than existing methods. We also report
real-data case studies to validate the applicability of dCAM in
biomedical research.

## Installation

You can install dCAM using `devtools`:

``` r
install.packages("devtools")
pak::pkg_install("cbil-vt/dCAM")
```

## Example

This is a basic example which shows the usage of dCAM. The synthetic
data constains 3 cell types, 7799 genes, and 20 sample (10 in group 1,
and 10 in group 2). We simulated a small portion of down- or
up-regulated genes based on a sorted cell bulk RNA-seq data (GSE73721).
This examples should take about 3 minutes to finish.

``` r
library(dCAM)

# A synthetic data
data(dataset)

# Run dCAM
fit <- dCAM::run_dcam(
  dataset,
  K = 3,
  cam_params = list(radius.thres = 0.96),
  reweight = "once",
  freeze_A = FALSE,
  max_iter = 5
)

# Performance evalution using pAUC
up <- compute_deg_metrics(fit$delta_S, dataset$truth$is_up,   0.05, "up")
dn <- compute_deg_metrics(fit$delta_S, dataset$truth$is_down, 0.05, "down")
cat(sprintf("  detection (pAUC):          up=%.4f  down=%.4f\n",
            up$pooled_pauc_norm, dn$pooled_pauc_norm))
```

In this example, for up-regulated genes, we have pAUC of 0.847, while
for down-regulated genes, the pAUC is 0.347. This performance difference
is expected, as down regulation leads to lower signal to noise ratio and
is more challenging to detect.
