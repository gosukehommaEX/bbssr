#' Range of the Nominal Total Sample Size Allowed by the Bounds
#'
#' Internal helper returning the smallest nominal total sample size whose group sizes of
#' \code{strat_sizes} are at least \code{low1} and \code{low2} in every stratum and whose
#' total is at least \code{N.min}, and the largest nominal size whose total does not exceed
#' \code{N.max}. The group sizes do not decrease with the nominal size, so both are found
#' by bisection.
#'
#' @param frac Fractions of the strata, summing to 1
#' @param r Allocation ratios of the strata
#' @param low1 Smallest sizes of group 1 in the strata
#' @param low2 Smallest sizes of group 2 in the strata
#' @param N.min Smallest final total sample size, or \code{NULL}
#' @param N.max Largest final total sample size, or \code{NULL}
#'
#' @return A numeric vector \code{c(lo, hi)}, with \code{hi = Inf} when \code{N.max} is
#'   \code{NULL}
#'
#' @keywords internal
#' @noRd
strat_bounds <- function(frac, r, low1, low2, N.min, N.max) {
  if (!is.null(N.min) && (length(N.min) != 1 || !is.numeric(N.min) || is.na(N.min))) {
    stop('N.min must be a single value')
  }
  if (!is.null(N.max) && (length(N.max) != 1 || !is.numeric(N.max) || is.na(N.max))) {
    stop('N.max must be a single value')
  }
  total <- function(n) {
    sz <- strat_sizes(n, frac, r)
    sum(sz$N1) + sum(sz$N2)
  }
  meets <- function(n) {
    sz <- strat_sizes(n, frac, r)
    all(sz$N1 >= low1) && all(sz$N2 >= low2) &&
      (is.null(N.min) || sum(sz$N1) + sum(sz$N2) >= N.min)
  }
  # Smallest n that meets the lower bounds: double until met, then bisect
  hi <- 1
  while (!meets(hi)) hi <- 2 * hi
  lo <- hi / 2
  if (hi == 1) {
    n.lo <- 1
  } else {
    while (hi - lo > 1) {
      mid <- floor((lo + hi) / 2)
      if (meets(mid)) hi <- mid else lo <- mid
    }
    n.lo <- hi
  }
  n.hi <- Inf
  if (!is.null(N.max)) {
    # The total of a nominal size n is at least n, so the largest admissible n is at most
    # N.max
    if (total(1) > N.max) {
      n.hi <- 0
    } else {
      a <- 1
      b <- floor(N.max) + 1
      while (b - a > 1) {
        mid <- floor((a + b) / 2)
        if (total(mid) <= N.max) a <- mid else b <- mid
      }
      n.hi <- a
    }
  }
  if (n.hi < n.lo) {
    stop('N.max is smaller than the smallest final sample size the design allows')
  }
  c(n.lo, n.hi)
}
