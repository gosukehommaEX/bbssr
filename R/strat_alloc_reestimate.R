#' Re-Estimated Total of a Stratified Design with Random Sizes
#'
#' Internal helper returning the final total sample size for each row of the interim
#' pooled numbers of responders \code{s}, when the strata hold \code{nS} interim patients.
#' The nominal total is that of \code{strat_reestimate}, with the fractions of the strata
#' observed at the interim analysis (\code{ss.frac = 'observed'}) or those of the plan.
#' The final total is the total of the group sizes that \code{strat_sizes} gives for the
#' nominal total, at least the interim total, at least the planned total \code{N.plan}
#' when \code{restricted} is \code{TRUE}, at least \code{N.min} and at most \code{N.max}.
#' An interim outcome without a sample size keeps the planned total, raised to the interim
#' total. The allocation of the remaining patients to the strata and the groups is left
#' to the allocation model.
#'
#' @param nS Interim sizes of the strata, all positive
#' @param re.args List of the arguments of \code{strat_reestimate}, with the planned
#'   fractions in \code{frac}
#' @param ss.frac \code{'observed'} or \code{'planned'}
#' @param N.plan Planned total sample size
#' @param restricted Logical. Whether the final total is at least \code{N.plan}
#' @param N.min Smallest final total, or \code{NULL}
#' @param N.max Largest final total, or \code{NULL}
#' @param s Interim pooled numbers of responders, a matrix with one column per stratum, or
#'   \code{NULL} (default) for every interim outcome of two strata, \code{s1} varying
#'   fastest
#' @inheritParams StratRR
#'
#' @return A list with the integer vector \code{N} of final totals, one per row of
#'   \code{s}, and the result \code{re} of \code{strat_reestimate}
#'
#' @keywords internal
#' @noRd
strat_alloc_reestimate <- function(nS, re.args, ss.frac, N.plan, restricted, N.min, N.max,
                                   alternative, tsmethod, s = NULL) {
  n <- sum(nS)
  if (is.null(s)) s <- as.matrix(expand.grid(0:nS[1], 0:nS[2]))
  fr <- if (ss.frac == 'observed') nS / n else re.args$frac
  # Only the totals of the strata enter the blinded estimates; the split between the
  # groups enters only the bounds of strat_reestimate, which are not used here
  n2 <- floor(nS / (1 + re.args$r))
  n1 <- nS - n2
  args <- re.args
  args$frac <- fr
  args[c('N1.plan', 'N2.plan', 'N.min', 'N.max')] <- list(n1, n2, NULL, NULL)
  args$restricted <- FALSE
  re <- do.call(strat_reestimate, c(list(s = s, n1 = n1, n2 = n2,
                                         alternative = alternative, tsmethod = tsmethod),
                                    args))
  ok <- !is.na(re$n.re)
  N <- rep(as.double(max(N.plan, n)), nrow(s))
  if (any(ok)) {
    sz <- strat_sizes(re$n.re[ok], fr, re.args$r)
    N[ok] <- rowSums(sz$N1) + rowSums(sz$N2)
  }
  N <- pmax(N, n, if (restricted) N.plan else 0, if (is.null(N.min)) 0 else N.min)
  if (!is.null(N.max)) N <- pmin(N, max(N.max, n))
  list(N = as.integer(N), re = re)
}
