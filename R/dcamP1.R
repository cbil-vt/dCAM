# =============================================================================
# PART 1 -- the four penalty weight rules
# Each takes the current difference and returns a weight per gene/cell type.
# =============================================================================

.keep_shape <- function(out, like) {
  if (is.null(dim(out)) && !is.null(dim(like))) dim(out) <- dim(like)
  out
}

# MCP: strong push for small differences, fading to nothing at gamma*lambda.
deriv_mcp <- function(delta, lambda, gamma = 3) {
  thr <- gamma * lambda
  .keep_shape(lambda * pmax(1 - abs(delta) / thr, 0), delta)
}

# TLP: full push while the difference is under tau, none once it is over.
deriv_tlp <- function(delta, lambda, tau = 1) {
  .keep_shape(lambda * (abs(delta) <= tau), delta)
}

# SCAD: flat push up to lambda, then fades, then stops at a*lambda.
deriv_scad <- function(delta, lambda, a = 3.7) {
  ad <- abs(delta)
  .keep_shape(ifelse(ad <= lambda, lambda,
                     pmax((a * lambda - ad) / (a - 1), 0)), delta)
}


deriv_afused <- function(delta, lambda, eps = 1e-6, majorize = FALSE) {
  ad <- abs(delta)
  w  <- if (isTRUE(majorize)) lambda / (2 * pmax(ad, eps) * (ad + eps))
  else                  lambda / (ad + eps)
  .keep_shape(w, delta)
}

# Pick the right weight rule for a penalty name.
get_weight_fun <- function(penalty, lambda, params) {
  switch(penalty,
         "MCP"    = function(d) deriv_mcp(d,    lambda, params$gamma %||% 3),
         "TLP"    = function(d) deriv_tlp(d,    lambda, params$tau   %||% 1),
         "SCAD"   = function(d) deriv_scad(d,   lambda, params$a     %||% 3.7),
         "AFUSED" = function(d) deriv_afused(d, lambda, params$eps %||% 1e-6,
                                             params$majorize %||% FALSE),
         stop("unknown penalty: ", penalty))
}


penalty_value <- function(delta, penalty, lambda,
                          gamma = 3, tau = 1, a = 3.7, eps = 1e-6,
                          form = c("effective", "true")) {
  form <- match.arg(form)
  ad <- abs(delta)


  if (form == "effective")
    return(switch(penalty,
                  "AFUSED" = sum(lambda * (ad - eps * log(1 + ad / eps))),
                  "MCP"    = { thr <- gamma * lambda
                  cap <- lambda * thr^2 / 2 - thr^3 / (3 * gamma)
                  sum(ifelse(ad <= thr,
                             lambda * ad^2 / 2 - ad^3 / (3 * gamma), cap)) },
                  "TLP"    = sum(ifelse(ad <= tau, lambda * ad^2 / 2, lambda * tau^2 / 2)),
                  "SCAD"   = { mid <- function(u)
                    lambda^3 / 2 + (1 / (a - 1)) *
                    (a * lambda * (u^2 - lambda^2) / 2 - (u^3 - lambda^3) / 3)
                  sum(ifelse(ad <= lambda, lambda * ad^2 / 2,
                             ifelse(ad <= a * lambda, mid(ad), mid(a * lambda)))) },
                  stop("unknown penalty: ", penalty)))

  switch(penalty,
         "MCP"    = sum(ifelse(ad <= gamma * lambda,
                               lambda * ad - ad^2 / (2 * gamma),
                               gamma * lambda^2 / 2)),
         "TLP"    = sum(lambda * pmin(ad, tau)),
         "SCAD"   = sum(ifelse(ad <= lambda, lambda * ad,
                               ifelse(ad <= a * lambda,
                                      (2 * a * lambda * ad - ad^2 - lambda^2) / (2 * (a - 1)),
                                      lambda^2 * (a + 1) / 2))),
         "AFUSED" = sum(lambda * log(1 + ad / eps)),
         stop("unknown penalty: ", penalty))
}
