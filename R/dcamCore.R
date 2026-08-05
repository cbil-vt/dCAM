
#' Title
#'
#' @param dataset
#' @param penalty
#' @param K
#' @param K_range
#' @param lambda_grid
#' @param lambda_rule
#' @param lambda_fixed
#' @param params
#' @param reweight
#' @param freeze_A
#' @param max_iter
#' @param tol
#' @param per_group
#' @param max_fpr
#' @param cam_init
#' @param verbose
#' @param cam_params
#' @param cam_seed
#' @param cam_dup_cor
#' @param bic_zero_rule
#' @param bic_zero_frac
#' @param bic_zero_abs
#'
#' @returns
#' @export
#'
#' @examples
run_dcam <- function(dataset,
                     penalty = "AFUSED",
                     K = NA,
                     K_range = 2:8,
                     lambda_grid = c(0.5, 1, 2, 4, 8),
                     lambda_rule = "bic",
                     lambda_fixed = 2,
                     params = list(gamma = 3, tau = 1, a = 3.7, eps = NA),
                     reweight = "iterative",
                     freeze_A = TRUE,
                     max_iter = 100,
                     tol = 1e-4,
                     per_group = TRUE,
                     max_fpr = 0.05,
                     cam_init = NULL,
                     verbose = TRUE,
                     cam_params = list(),
                     cam_seed = 2024,
                     cam_dup_cor = 0.95,
                     bic_zero_rule = "relative",
                     bic_zero_frac = 1e-3,
                     bic_zero_abs = 1e-8) {

  s  <- split_groups(dataset$X, dataset$group)
  X1 <- s$X1; X2 <- s$X2

  # --- starting point: CAM ------------------------------------------------
  if (is.null(cam_init)) {
    if (verbose) cat("  [dCAM] running CAM for the starting point ...\n")
    cam_init <- run_cam_init(X1, X2, K = K, K_range = K_range,
                             per_group = per_group,
                             cam_params = cam_params, seed = cam_seed,
                             dup_cor = cam_dup_cor)
  } else if (verbose) cat("  [dCAM] reusing the CAM starting point provided\n")
  if (verbose) cat(sprintf("  [dCAM] using K = %d cell types\n", cam_init$K))

  # --- AFUSED's eps: pick it from the data if not given -------------------

  if (identical(penalty, "AFUSED") && (is.null(params$eps) || is.na(params$eps))) {
    d0 <- abs(cam_init$S2 - cam_init$S1)
    params$eps <- max(stats::quantile(d0[d0 > 0], 0.10, na.rm = TRUE), 1e-6)
    if (verbose) cat(sprintf("  [dCAM] AFUSED eps chosen from data: %.4g\n", params$eps))
  }

  # --- choose lambda and fit ----------------------------------------------
  if (verbose)
    cat(
      sprintf(
        "  [dCAM] penalty=%s  tuning rule=%s  freeze_A=%s\n",
        penalty,
        lambda_rule,
        freeze_A
      )
    )
  tuned <- tune_lambda(
    X1,
    X2,
    cam_init$S1,
    cam_init$S2,
    cam_init$A1,
    cam_init$A2,
    penalty = penalty,
    lambda_grid = lambda_grid,
    params = params,
    rule = lambda_rule,
    lambda_fixed = lambda_fixed,
    truth = dataset$truth,
    max_iter = max_iter,
    tol = tol,
    reweight = reweight,
    freeze_A = freeze_A,
    max_fpr = max_fpr,
    verbose = verbose,
    bic_zero_rule = bic_zero_rule,
    bic_zero_frac = bic_zero_frac,
    bic_zero_abs = bic_zero_abs
  )
  f <- tuned$fit

  if (isTRUE(verbose) && !is.null(f$n_increase) && f$n_increase > 0)
    cat(sprintf("  [dCAM] note: the reported objective rose in %d of %d rounds.\n",
                f$n_increase, f$n_iter),
        "         Inspect fit$objective. With the default settings the trace should\n",
        "         fall monotonically; repeated rises point at a numerical problem\n",
        "         (check for the 'non-finite value' notes above).\n")

  fit <- make_fit("dCAM", f$S1, f$S2, f$A1, f$A2,
                  params = list(penalty = penalty, lambda = tuned$lambda,
                                lambda_rule = tuned$rule, freeze_A = freeze_A,
                                reweight = reweight, K = cam_init$K,
                                eps = params$eps,
                                majorize = params$majorize %||% FALSE),
                  extra = list(lambda_table = tuned$table,
                               objective = f$objective, n_iter = f$n_iter,
                               n_increase = f$n_increase,
                               objective_form = f$objective_form,
                               cam_init = cam_init,
                               cam_params = cam_params))

  # Line the cell types up with the truth, when we have it, so that "cell type 1"
  # means the same thing for every method during evaluation.
  if (!is.null(dataset$truth)) fit <- align_fit_to_reference(fit, dataset$truth$S1)
  fit
}
