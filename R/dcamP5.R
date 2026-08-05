# =============================================================================
# PART 5 -- choosing lambda (the penalty strength)
# =============================================================================

# -----------------------------------------------------------------------------
# dcam_bic() -- an information criterion for one fitted model.

dcam_bic <- function(X1, X2, fit,
                     zero_rule = c("relative", "absolute"),
                     zero_frac = 1e-3, zero_abs = 1e-8) {
  zero_rule <- match.arg(zero_rule)
  n   <- length(X1) + length(X2)
  rss <- sum((X1 - fit$S1 %*% fit$A1)^2) + sum((X2 - fit$S2 %*% fit$A2)^2)
  ad  <- abs(fit$delta_S)
  thr <- if (zero_rule == "relative") {
    nz <- ad[ad > 0 & is.finite(ad)]
    if (length(nz)) zero_frac * stats::median(nz) else zero_abs
  } else zero_abs
  df  <- sum(ad > thr)                         # how many changes are claimed
  structure(n * log(rss / n) + log(n) * df,
            df = df, rss = rss, threshold = thr, zero_rule = zero_rule)
}


# -----------------------------------------------------------------------------
# tune_lambda() -- fit at every lambda in the grid and pick one.
#
#   rule = "bic"   : lowest BIC. Works without truth -> use this on real data.
#   rule = "truth" : best pAUC against the planted truth.
#   rule = "fixed" : no tuning, just use lambda_fixed.
#

# -----------------------------------------------------------------------------

#' Title
#'
#' @param X1
#' @param X2
#' @param S1i
#' @param S2i
#' @param A1i
#' @param A2i
#' @param penalty
#' @param lambda_grid
#' @param params
#' @param rule
#' @param lambda_fixed
#' @param truth
#' @param max_iter
#' @param tol
#' @param reweight
#' @param freeze_A
#' @param max_fpr
#' @param verbose
#' @param bic_zero_rule
#' @param bic_zero_frac
#' @param bic_zero_abs
#'
#' @returns
#' @export
#'
#' @examples
tune_lambda <- function(X1, X2, S1i, S2i, A1i, A2i,
                        penalty = "AFUSED", lambda_grid = c(0.5, 1, 2, 4, 8),
                        params = list(), rule = c("bic", "truth", "fixed"),
                        lambda_fixed = 2, truth = NULL,
                        max_iter = 100, tol = 1e-4,
                        reweight = "iterative", freeze_A = TRUE,
                        max_fpr = 0.05, verbose = TRUE,
                        bic_zero_rule = "relative",
                        bic_zero_frac = 1e-3, bic_zero_abs = 1e-8) {
  rule <- match.arg(rule)

  if (rule == "fixed") {
    fit <- dcam_refine(X1, X2, S1i, S2i, A1i, A2i, penalty, lambda_fixed,
                       params, max_iter, tol, reweight, freeze_A)
    return(list(fit = fit, lambda = lambda_fixed, rule = rule, table = NULL))
  }
  if (rule == "truth" && is.null(truth))
    stop("rule='truth' needs the planted truth, which real data does not have.\n",
         "  Use rule='bic' instead.")

  rows <- list(); fits <- list()
  for (lam in lambda_grid) {
    if (verbose) cat(sprintf("    lambda = %-6g ", lam))
    f <- dcam_refine(X1, X2, S1i, S2i, A1i, A2i, penalty, lam,
                     params, max_iter, tol, reweight, freeze_A)
    bic <- dcam_bic(X1, X2, f, bic_zero_rule, bic_zero_frac, bic_zero_abs)
    sc  <- NA_real_
    if (!is.null(truth)) {
      # browser()
      f <- align_fit_to_reference(f, truth$S1)  # align to truth
      up <- compute_deg_metrics(f$delta_S, truth$is_up,   max_fpr, "up")$pooled_pauc_norm
      dn <- compute_deg_metrics(f$delta_S, truth$is_down, max_fpr, "down")$pooled_pauc_norm
      sc <- mean(c(up, dn))
    }
    rows[[length(rows) + 1]] <- data.frame(lambda = lam, bic = as.numeric(bic),
                                           df = attr(bic, "df"),
                                           nonzero_frac = attr(bic, "df") / length(f$delta_S),
                                           pauc_mean = sc)
    fits[[as.character(lam)]] <- f
    if (verbose) cat(sprintf("BIC=%12.1f  df=%7d  pAUC=%s\n", as.numeric(bic),
                             attr(bic, "df"),
                             ifelse(is.na(sc), "  --  ", sprintf("%.4f", sc))))
  }
  tab  <- do.call(rbind, rows)
  # (added) if df never changes, the BIC penalty term cannot discriminate and
  # "bic" is really just choosing the smallest lambda.
  if (rule == "bic" && length(unique(tab$df)) == 1L)
    warning("BIC df is identical at every lambda (", tab$df[1],
            ") -- the penalty term is not discriminating, so BIC is effectively\n",
            "  selecting the smallest lambda. Lower BIC_ZERO_FRAC in 00_config.R.",
            call. = FALSE)
  best <- if (rule == "bic") tab$lambda[which.min(tab$bic)]
  else               tab$lambda[which.max(tab$pauc_mean)]
  if (verbose) cat(sprintf("    -> chose lambda = %g  (rule: %s)\n", best, rule))

  list(fit = fits[[as.character(best)]], lambda = best, rule = rule, table = tab)
}
