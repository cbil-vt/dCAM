# =============================================================================
# R/utils.R  --  small helpers used by several scripts
# =============================================================================
# fitting a signature given proportions, lining up cell types between two
# estimates, and saving/loading results.
# =============================================================================


# `a %||% b` -> returns a, unless a is NULL, in which case it returns b.
`%||%` <- function(a, b) if (is.null(a)) b else a


# -----------------------------------------------------------------------------
# fit_S_ols() -- given the bulk X and the proportions A, estimate the signature S
#                by ordinary least squares.
#   X is genes x samples, A is cellTypes x samples, result S is genes x cellTypes.
#   The tiny value added to the diagonal keeps the maths stable when two cell
#   types have very similar proportions.
# -----------------------------------------------------------------------------
fit_S_ols <- function(X, A) {
  AAt <- A %*% t(A)
  S <- X %*% t(A) %*% solve(AAt + 1e-8 * diag(nrow(AAt)))
  dimnames(S) <- list(rownames(X), rownames(A))
  S
}


# -----------------------------------------------------------------------------
# fit_S_nnls() -- same idea, but the signature is forced to be non-negative
#                 (expression cannot be below zero).
# -----------------------------------------------------------------------------
fit_S_nnls <- function(X, A) {
  if (!requireNamespace("nnls", quietly = TRUE))
    stop("Package 'nnls' is required. Install it with install.packages('nnls').")
  S <- t(apply(X, 1, function(y)
    tryCatch(nnls::nnls(t(A), y)$x,
             error = function(e) pmax(qr.solve(t(A), y), 0))))
  dimnames(S) <- list(rownames(X), rownames(A))
  S
}


# -----------------------------------------------------------------------------
# fit_A_nnls() -- given the bulk X and a signature S, estimate the proportions A.
#                 Non-negative, and every sample's proportions are made to sum to 1.
# -----------------------------------------------------------------------------
fit_A_nnls <- function(S, X) {
  if (!requireNamespace("nnls", quietly = TRUE))
    stop("Package 'nnls' is required. Install it with install.packages('nnls').")
  A <- apply(X, 2, function(y)
    tryCatch(nnls::nnls(S, y)$x,
             error = function(e) pmax(qr.solve(S, y), 0)))
  A <- matrix(A, nrow = ncol(S))
  cs <- colSums(A); cs[cs == 0] <- 1        # avoid dividing by zero
  A <- sweep(A, 2, cs, "/")
  dimnames(A) <- list(colnames(S), colnames(X))
  A
}


# -----------------------------------------------------------------------------
# best_permutation() -- cell types come out of a method in an arbitrary order.

# -----------------------------------------------------------------------------
#' Title
#'
#' @param M_est
#' @param M_ref
#'
#' @returns
#' @export
#'
#' @examples
best_permutation <- function(M_est, M_ref) {
  K <- ncol(M_ref)
  cost <- suppressWarnings(cor(M_est, M_ref))
  cost[!is.finite(cost)] <- 0
  # Try the proper assignment solver first; fall back to a simple greedy pick.
  if (requireNamespace("clue", quietly = TRUE)) {
    as.integer(clue::solve_LSAP(1 - cost, maximum = FALSE))
  } else {
    perm <- integer(K); used <- rep(FALSE, K)
    for (j in seq_len(K)) {                 # for each reference column...
      cand <- which(!used)                  # ...pick the best unused estimate
      pick <- cand[which.max(cost[cand, j])]
      perm[pick] <- j; used[pick] <- TRUE
    }
    perm
  }
}


# -----------------------------------------------------------------------------
# align_fit_to_reference() -- reorder a fit's cell types to match a reference,

# -----------------------------------------------------------------------------
align_fit_to_reference <- function(fit, S_ref) {
  perm <- best_permutation(fit$S1, S_ref)
  ord  <- order(perm)                       # the column order that lines things up
  fit$S1 <- fit$S1[, ord, drop = FALSE]
  fit$S2 <- fit$S2[, ord, drop = FALSE]
  fit$A1 <- fit$A1[ord, , drop = FALSE]
  fit$A2 <- fit$A2[ord, , drop = FALSE]
  colnames(fit$S1) <- colnames(fit$S2) <- colnames(S_ref)
  rownames(fit$A1) <- rownames(fit$A2) <- colnames(S_ref)
  fit$delta_S <- fit$S2 - fit$S1
  fit
}


# -----------------------------------------------------------------------------
# make_fit() --
# -----------------------------------------------------------------------------
make_fit <- function(method, S1, S2, A1, A2, params = list(), extra = list()) {
  structure(c(list(method  = method,
                   S1      = S1, S2 = S2,
                   A1      = A1, A2 = A2,
                   delta_S = S2 - S1,
                   params  = params),
              extra),
            class = c("dcam_fit", "list"))
}


# Friendly one-line summary when you print a fit at the R prompt.
print.dcam_fit <- function(x, ...) {
  cat(sprintf("<fit> %s | %d genes x %d cell types | %d + %d samples\n",
              x$method, nrow(x$S1), ncol(x$S1), ncol(x$A1), ncol(x$A2)))
  if (length(x$params)) {
    p <- paste(names(x$params), unlist(lapply(x$params, function(v)
      paste(format(v, digits = 4), collapse = ","))), sep = "=", collapse = ", ")
    cat("  settings:", p, "\n")
  }
  invisible(x)
}


# -----------------------------------------------------------------------------
# split_groups() -- cut the combined bulk matrix into the two groups.
# -----------------------------------------------------------------------------
split_groups <- function(X, group) {
  g <- sort(unique(group))
  if (length(g) != 2) stop("expected exactly 2 groups, found ", length(g))
  list(X1 = X[, group == g[1], drop = FALSE],
       X2 = X[, group == g[2], drop = FALSE])
}



# -----------------------------------------------------------------------------
# save_rds() / load_rds() -- save and load with a short message, so you always
# know where a file went.
# -----------------------------------------------------------------------------
save_rds <- function(obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, path)
  cat("  saved:", path, "\n")
  invisible(path)
}

load_rds <- function(path) {
  if (!file.exists(path)) stop("file not found: ", path)
  cat("  loaded:", path, "\n")
  readRDS(path)
}
