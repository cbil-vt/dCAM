# ---------------------------------------------------------------------------
# >>> EDIT THIS ONE LINE <<<  the folder that contains this file
# ---------------------------------------------------------------------------
PROJECT_DIR <- 'D:/work/output/dcam_exp/'

# ---------------------------------------------------------------------------
# CAM3 SETTINGS
# ---------------------------------------------------------------------------
# Passed straight through to CAM3::CAM3Run().
CAM3_PARAMS <- list(
  dim.rdc        = 6,
  thres.low      = 0.05,
  thres.high     = 1.00,
  cluster.method = "Fixed-Radius",
  radius.thres   = 0.96,             # cluster radius
  sim.thres      = 0.98,             # merge clusters more similar than this
  MG.num.thres   = 20,               # a marker cluster needs at least this many genes
  fast.mode      = FALSE             # FALSE = full search
)
CAM_SEED <- 2024
CAM_DUP_COR_THRES <- 0.95


# ---------------------------------------------------------------------------
# dCAM SETTINGS  (used by 02_run_dcam.R)
# ---------------------------------------------------------------------------
# How many cell types to look for. Set to NA to let CAM choose automatically
# using MDL (Minimum Description Length) over the range below.
K_USE     <- 3
K_RANGE   <- 2:8        # only used when K_USE is NA

# Which penalty to use :
#   "AFUSED" = adaptive fused lasso (default)
#   "MCP", "TLP", "SCAD"
PENALTY   <- "AFUSED"

# Penalty strength. Give a VECTOR to tune over it, or a single number to fix it.
LAMBDA_GRID <- c(1, 2, 3)
# How to pick lambda from that grid:
#   "bic"    = information criterion -- works WITHOUT truth, so usable on real data
#   "truth"  = pick the lambda with the best pAUC -- SIMULATION ONLY (needs truth)
#   "fixed"  = No tuning; use LAMBDA_FIXED below
LAMBDA_RULE  <- "truth"
LAMBDA_FIXED <- 2

# How dcam_bic() decides an entry of S2 - S1 counts as a claimed change (added).
# The quadratic penalty shrinks differences toward zero without setting them
# exactly to zero, so a fixed 1e-8 cutoff counts almost everything and the BIC
# penalty term stops discriminating between lambdas.
#   "relative" = threshold is BIC_ZERO_FRAC * median(|S2 - S1|) over non-zeros
#   "absolute" = threshold is BIC_ZERO_ABS
BIC_ZERO_RULE <- "relative"
BIC_ZERO_FRAC <- 1e-3
BIC_ZERO_ABS  <- 1e-8

# Extra settings for the penalties .
PENALTY_PARAMS <- list(
  gamma = 3,      # MCP shape
  tau   = 1,      # TLP shape
  a     = 3.7,    # SCAD shape (3.7 is the standard choice)
  eps   = NA,     # AFUSED: NA = choose automatically from the data
  majorize = FALSE
)

# Keep the proportions A fixed at CAM's estimate and refine the SIGNATURES only.
#   TRUE: CAM estimates A very accurately. Re-estimating A from the
#     penalized S pulls A away from the truth, and the next
#     S update inherits that damage.
#   FALSE: original behaviour (alternate S and A updates).
FREEZE_A  <- FALSE

# REWEIGHT decides whether the penalty weights are recomputed each round.
#   "once"
#   "iterative"
REWEIGHT  <- "once"  # "iterative" = recompute penalty weights each round
# "once"      = weights fixed at start
MAX_ITER  <- 5          # most rounds of refinement, 100
TOL       <- 1e-4         # stop criteria


# ---------------------------------------------------------------------------
# EVALUATION SETTINGS  (used by 04_evaluate.R)
# ---------------------------------------------------------------------------
MAX_FPR <- 0.05          # ROC is summarised over the low false-positive region

# Which ways of ranking genes to evaluate:
#   "effsize"   = rank by the size of the estimated change (S2 - S1)
#   "pvalue"    = rank by a resampling t-test
#   "stability" = rank by how OFTEN a gene is picked across resamples
SCORE_MODES <- c("effsize")

# --- resampling (shared by "pvalue" and "stability") -------------------------
RESAMPLE_SCHEME <- "loo"   # "loo" = leave-one-out (recommended)
# "bootstrap" = resample with replacement
N_RESAMPLE      <- 4     # number of resamples when scheme = "bootstrap"
N_RESAMPLE_DCAM <- 4      # fewer for dCAM, because each round refits the penalty


# --- p-value mode only -----------------------------------------------------
PVAL_ENGINE  <- "paired_ttest"  # "paired_ttest" = ordinary t.test() on the resampled
#   values (simple; p-values are for RANKING only)
# "wald"       = t = delta / SE (p-values do not
#   depend on the number of resamples)
S0_QUANTILE  <- 0.5             # small stabiliser added to the SE

# --- stability mode only ---------------------------------------------------
STAB_TOP_FRAC   <- 0.05              # each run keeps this fraction of genes per cell type
STAB_THRESHOLDS <- c(0.6, 0.75, 0.9) # report the stable-gene table at each of these

