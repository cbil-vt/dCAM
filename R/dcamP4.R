# =============================================================================
# PART 4 -- getting a starting point (CAM), including choosing K
# =============================================================================


run_cam_init <- function(X1, X2, K = NA, K_range = 2:8, per_group = TRUE,
                         cam_params = list(), seed = 2024, dup_cor = 0.95) {

  cam_defaults <- list(dim.rdc = 10,
                       thres.low = 0.05, thres.high = 1.00,
                       cluster.method = "Fixed-Radius",
                       radius.thres = 0.95, sim.thres = 0.95,
                       MG.num.thres = 20, fast.mode = FALSE)
  cp <- utils::modifyList(cam_defaults, as.list(cam_params))

  # Ask CAM for a single K, or for a whole range (then MDL picks).
  K_arg <- if (is.na(K)) K_range else K

  run_one <- function(X) {
    set.seed(seed)
    rc <- do.call(CAM3Run, c(list(data = X, K = K_arg), cp))
    # Which K did it settle on? With a single K there is only one answer;
    # with a range, CAM3 has already scored them by MDL.
    k_used <- if (length(K_arg) == 1) K_arg else pick_K_mdl(rc, K_arg)
    res <- rc@ASestResult[[as.character(k_used)]]

    Ar  <- res@Aest
    Akx <- if (nrow(Ar) == ncol(X)) t(Ar) else Ar     # want cellTypes x samples

    # CAM gives the signature directly; only fall back if that slot is missing.
    Sk <- tryCatch({
      s <- res@Sest
      if (is.null(s)) stop("no Sest")
      if (nrow(s) != nrow(X)) s <- t(s)
      s
    }, error = function(e) fit_S_ols(X, Akx))

    dimnames(Sk)  <- list(rownames(X), paste0("cell", seq_len(ncol(Sk))))
    dimnames(Akx) <- list(colnames(Sk), colnames(X))

    # Clean any non-finite values CAM may have produced, so they cannot poison
    # everything downstream.
    Sk  <- sanitize_matrix(Sk,  "CAM signature")
    Akx <- sanitize_matrix(Akx, "CAM proportions")


    if (ncol(Sk) > 1) {
      cc <- suppressWarnings(cor(Sk)); diag(cc) <- 0
      if (max(abs(cc), na.rm = TRUE) > dup_cor)
        warning("CAM's cell types look nearly identical (max correlation ",
                sprintf("%.2f", max(abs(cc), na.rm = TRUE)),
                "). It probably could not separate ", ncol(Sk),
                " cell types from ", ncol(X), " samples.\n",
                "  Try more samples per group, or a smaller K.",
                call. = FALSE)
    }
    list(S = Sk, A = Akx, K = k_used, cam_obj = rc)
  }

  if (per_group) {
    g1 <- run_one(X1)
    g2 <- run_one(X2)
    # Line the two groups' cell types up with each other.
    perm <- best_permutation(g2$S, g1$S); ord <- order(perm)
    list(S1 = g1$S, S2 = g2$S[, ord, drop = FALSE],
         A1 = g1$A, A2 = g2$A[ord, , drop = FALSE],
         K  = g1$K, cam = list(g1 = g1$cam_obj, g2 = g2$cam_obj))
  } else {
    both <- run_one(cbind(X1, X2))
    n1 <- ncol(X1)
    list(S1 = both$S, S2 = both$S,
         A1 = both$A[, seq_len(n1),               drop = FALSE],
         A2 = both$A[, seq_len(ncol(X2)) + n1,    drop = FALSE],
         K  = both$K, cam = list(both = both$cam_obj))
  }
}


# -----------------------------------------------------------------------------
# pick_K_mdl() -- read the MDL scores out of a CAM result and return the K with
# the lowest (best) score.
# -----------------------------------------------------------------------------
pick_K_mdl <- function(cam_obj, K_range) {
  mdl <- tryCatch({
    m <- cam_obj@MDLResult
    # Different CAM versions expose this slightly differently -- try the usual
    # names and take whichever yields one number per K.
    v <- NULL
    for (nm in c("MDLs", "mdls", "MDL", "mdl")) {
      if (.hasSlot(m, nm)) { cand <- methods::slot(m, nm)
      if (length(cand) == length(K_range)) { v <- as.numeric(cand); break } }
    }
    v
  }, error = function(e) NULL)

  if (is.null(mdl) || !length(mdl)) {
    k <- K_range[ceiling(length(K_range) / 2)]
    message("    (MDL scores unavailable -- defaulting to K = ", k,
            "; set K_USE explicitly if you know it)")
    return(k)
  }
  k <- K_range[which.min(mdl)]
  cat(sprintf("    MDL chose K = %d  (scores: %s)\n", k,
              paste(sprintf("%d:%.1f", K_range, mdl), collapse = " ")))
  k
}
