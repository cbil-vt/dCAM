# =============================================================================
# PART 2 -- the two update steps
# =============================================================================

# -----------------------------------------------------------------------------
# update_S_penalized()
# Re-estimate S1 and S2 TOGETHER, one gene at a time.

# -----------------------------------------------------------------------------
update_S_penalized <- function(X1, X2, A1, A2, W_pen) {
  G  <- nrow(X1); K <- nrow(A1)
  M1 <- ncol(X1); M2 <- ncol(X2)
  S_out <- matrix(0, nrow = G, ncol = 2 * K)   # cols 1..K = S1, K+1..2K = S2

  # These blocks do not change from gene to gene -- build them once.
  top_left  <- t(A1)                            # "S1 explains X1"
  top_right <- matrix(0, nrow = M1, ncol = K)
  mid_left  <- matrix(0, nrow = M2, ncol = K)
  mid_right <- t(A2)                            # "S2 explains X2"
  zero_pen  <- numeric(K)                       # the penalty rows aim at 0

  use_fcnnls <- requireNamespace("NMF", quietly = TRUE)

  for (i in seq_len(G)) {
    sw <- sqrt(pmax(W_pen[i, ], 0))             # sqrt because the maths needs it
    bot_l <- diag(sw, nrow = K); bot_r <- -bot_l

    Xi <- rbind(cbind(top_left, top_right),
                cbind(mid_left, mid_right),
                cbind(bot_l,    bot_r))
    yi <- c(X1[i, ], X2[i, ], zero_pen)

    sol <- if (use_fcnnls) {
      tryCatch(as.numeric(NMF::fcnnls(x = Xi, y = matrix(yi, ncol = 1))$x),
               error = function(e) pmax(tryCatch(qr.solve(Xi, yi),
                                                 error = function(e2) rep(0, 2 * K)), 0))
    } else if (requireNamespace("nnls", quietly = TRUE)) {
      tryCatch(nnls::nnls(Xi, yi)$x,
               error = function(e) pmax(tryCatch(qr.solve(Xi, yi),
                                                 error = function(e2) rep(0, 2 * K)), 0))
    } else {
      pmax(tryCatch(qr.solve(Xi, yi), error = function(e) rep(0, 2 * K)), 0)
    }
    S_out[i, ] <- sol
  }

  list(S1 = S_out[, 1:K,          drop = FALSE],
       S2 = S_out[, (K + 1):(2 * K), drop = FALSE])
}


# -----------------------------------------------------------------------------
# dcam_objective() -- how good is the current answer?

dcam_objective <- function(X1, X2, S1, S2, A1, A2, penalty, lambda,
                           gamma = 3, tau = 1, a = 3.7, eps = 1e-6,
                           form = "effective", W = NULL) {
  fit1 <- sum((X1 - S1 %*% A1)^2)
  fit2 <- sum((X2 - S2 %*% A2)^2)
  0.5 * (fit1 + fit2) +
    (if (identical(form, "weighted_l2") && !is.null(W))
      0.5 * sum(W * (S1 - S2)^2, na.rm = TRUE)
     else penalty_value(S1 - S2, penalty, lambda, gamma, tau, a, eps, form))
}
