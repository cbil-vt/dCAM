# =============================================================================
# PART 3 -- the refinement loop
# =============================================================================

# -----------------------------------------------------------------------------
# sanitize_matrix() -- replace any NaN / Inf / NA with 0 and say so.

# -----------------------------------------------------------------------------
sanitize_matrix <- function(M, name = "matrix", verbose = TRUE) {
  bad <- !is.finite(M)
  if (any(bad)) {
    if (verbose)
      message(sprintf("    note: %s had %d non-finite value(s) (%.3f%%) -- set to 0",
                      name, sum(bad), 100 * mean(bad)))
    M[bad] <- 0
  }
  M
}


# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
dcam_refine <- function(X1, X2, S1, S2, A1, A2,
                        penalty = "AFUSED", lambda = 1,
                        params = list(gamma = 3, tau = 1, a = 3.7, eps = 1e-6),
                        max_iter = 100, tol = 1e-4,
                        reweight = c("iterative", "once"),
                        freeze_A = TRUE, verbose = FALSE) {
  reweight   <- match.arg(reweight)
  weight_fun <- get_weight_fun(penalty, lambda, params)


  obj_form <- if (reweight == "once") "weighted_l2"
  else if (isTRUE(params$majorize)) "true" else "effective"


  S1 <- sanitize_matrix(S1, "starting S1", verbose)
  S2 <- sanitize_matrix(S2, "starting S2", verbose)
  A1 <- sanitize_matrix(A1, "starting A1", verbose)
  A2 <- sanitize_matrix(A2, "starting A2", verbose)

  # If the weights are frozen, work them out now from the starting difference.
  W_fixed <- if (reweight == "once") pmax(weight_fun(S1 - S2), 0) else NULL

  obj_prev  <- dcam_objective(X1, X2, S1, S2, A1, A2, penalty, lambda,
                              params$gamma %||% 3, params$tau %||% 1,
                              params$a %||% 3.7,   params$eps %||% 1e-6,
                              obj_form, W_fixed)
  obj_trace <- obj_prev
  n_up <- 0L                 # (added) how many rounds the objective went UP
  t <- 0L

  for (t in seq_len(max_iter)) {
    # (1) turn the current difference into per-gene push strengths
    W_pen <- if (reweight == "once") W_fixed else pmax(weight_fun(S1 - S2), 0)
    W_pen[!is.finite(W_pen)] <- 0            # a zero weight = "do not push"

    # (2) re-estimate both signatures together
    su <- update_S_penalized(X1, X2, A1, A2, W_pen)
    S1 <- sanitize_matrix(su$S1, "S1", verbose = FALSE)
    S2 <- sanitize_matrix(su$S2, "S2", verbose = FALSE)

    # (3) re-estimate the proportions -- skipped when freeze_A = TRUE
    if (!freeze_A) {
      A1 <- sanitize_matrix(fit_A_nnls(S1, X1), "A1", verbose = FALSE)
      A2 <- sanitize_matrix(fit_A_nnls(S2, X2), "A2", verbose = FALSE)
    }

    # (4) has it stopped improving?
    obj <- dcam_objective(X1, X2, S1, S2, A1, A2, penalty, lambda,
                          params$gamma %||% 3, params$tau %||% 1,
                          params$a %||% 3.7,   params$eps %||% 1e-6,
                          obj_form, W_pen)
    obj_trace <- c(obj_trace, obj)

    # Guard the stopping test: if either objective is not a usable number we
    # stop here rather than crashing, and keep the best answer we have.
    if (!is.finite(obj) || !is.finite(obj_prev)) {
      if (verbose)
        message("    note: objective became non-finite -- stopping at iteration ", t)
      break
    }
    if (obj > obj_prev * (1 + 1e-12)) n_up <- n_up + 1L
    rel <- abs(obj_prev - obj) / max(abs(obj_prev), 1e-12)
    if (verbose) cat(sprintf("    iter %2d  obj=%.4e  change=%.2e\n", t, obj, rel))
    if (!is.finite(rel) || rel < tol) break
    obj_prev <- obj
  }

  dimnames(S1) <- dimnames(S2) <- list(rownames(X1), rownames(A1))
  list(S1 = S1, S2 = S2, A1 = A1, A2 = A2, delta_S = S2 - S1,
       objective = obj_trace, n_iter = t, n_increase = n_up,
       objective_form = obj_form)
}
