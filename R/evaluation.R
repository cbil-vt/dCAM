#
auc_proc <- function(score, label, max_fpr = 0.05) {
  ord <- order(score, decreasing = TRUE)
  y   <- label[ord]
  P   <- sum(y == 1); N <- sum(y == 0)
  if (P == 0 || N == 0)
    return(list(fpr = NA, tpr = NA, max_fpr = max_fpr,
                pauc_raw = NA_real_, pauc_norm = NA_real_))

  tpr <- c(0, cumsum(y == 1) / P)   # share of real genes found so far
  fpr <- c(0, cumsum(y == 0) / N)   # share of non-genes wrongly called

  # Keep only the part of the curve below max_fpr.
  keep  <- which(fpr <= max_fpr)
  fpr_t <- fpr[keep]; tpr_t <- tpr[keep]

  # If the kept part stops early, add one interpolated point exactly at max_fpr
  # so every method's area is measured up to the same boundary.
  if (tail(fpr_t, 1) < max_fpr) {
    nx <- max(keep) + 1
    if (nx <= length(fpr) && (fpr[nx] - fpr[nx - 1]) > 0) {
      slope <- (tpr[nx] - tpr[nx - 1]) / (fpr[nx] - fpr[nx - 1])
      fpr_t <- c(fpr_t, max_fpr)
      tpr_t <- c(tpr_t, tpr[nx - 1] + slope * (max_fpr - fpr[nx - 1]))
    }
  }

  # Area by the trapezoid rule, then divided by the box area so 0..1.
  pauc_raw <- sum(diff(fpr_t) * (head(tpr_t, -1) + tail(tpr_t, -1)) / 2)
  list(fpr = fpr_t, tpr = tpr_t, max_fpr = max_fpr,
       pauc_raw = pauc_raw, pauc_norm = pauc_raw / max_fpr)
}


# compute_deg_metrics() -- score a whole genes x cellTypes matrix at once.
#' Title
#'
#' @param score_mat
#' @param truth_mat
#' @param max_fpr
#' @param direction
#'
#' @returns
#' @export
#'
#' @examples
compute_deg_metrics <- function(score_mat, truth_mat, max_fpr = 0.05,
                                direction = c("up", "down", "both")) {
  direction <- match.arg(direction)
  s <- switch(direction, up = score_mat, down = -score_mat, both = abs(score_mat))
  curve <- auc_proc(as.numeric(s), as.numeric(truth_mat) == 1, max_fpr)

  # Also do it cell type by cell type, in case one cell type behaves oddly.
  per_ct <- vapply(seq_len(ncol(score_mat)), function(k) {
    cc <- auc_proc(s[, k], truth_mat[, k] == 1, max_fpr)
    cc$pauc_norm %||% NA_real_
  }, numeric(1))
  names(per_ct) <- colnames(score_mat) %||% paste0("cell", seq_len(ncol(score_mat)))

  list(pooled_pauc_norm = curve$pauc_norm,
       pooled_curve     = curve,
       per_celltype     = per_ct,
       direction        = direction,
       n_positive       = sum(truth_mat))
}


# estimation_quality() -- how close is the fit to the truth, ignoring ROC?
estimation_quality <- function(fit, truth) {
  if (is.null(truth)) return(NULL)
  sc <- mean(c(diag(suppressWarnings(cor(fit$S1, truth$S1))),
               diag(suppressWarnings(cor(fit$S2, truth$S2)))), na.rm = TRUE)
  ac <- mean(c(diag(suppressWarnings(cor(t(fit$A1), t(truth$A1)))),
               diag(suppressWarnings(cor(t(fit$A2), t(truth$A2))))), na.rm = TRUE)
  dc <- suppressWarnings(cor(as.numeric(fit$delta_S),
                             as.numeric(truth$S2 - truth$S1)))
  data.frame(S_cor = sc, A_cor = ac, delta_cor = dc)
}
