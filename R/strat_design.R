#' Validate the Settings of a Stratified Re-Estimation Design
#'
#' Internal helper shared by \code{\link{StratPowerBSSR}}, \code{\link{StratTypeIErrorBSSR}}
#' and \code{\link{StratAlphaAdjBSSR}}. It checks the planned group sizes, the test, the
#' fractions of the strata and the allocation ratios, derives the interim group sizes, and
#' checks the settings that depend on the number of strata.
#'
#' @param N1 Planned sizes of group 1 in the strata
#' @param N2 Planned sizes of group 2 in the strata
#' @param r Allocation ratios, a single value or one per stratum
#' @param frac Fractions of the strata, or \code{NULL} for those of the planned design
#' @param omega Fraction of the planned sizes at the interim analysis, or \code{NULL}
#' @param n.interim Interim group sizes, or \code{NULL}
#' @param ss.method Re-estimation rule
#' @param nsim Number of simulated trials, used with three or more strata
#' @inheritParams StratRR
#'
#' @return A list with the result \code{a} of \code{check_strat_args}, the number \code{J}
#'   of strata, \code{frac}, \code{r}, the interim sizes \code{n1} and \code{n2}, and
#'   \code{nsim}
#'
#' @keywords internal
#' @noRd
strat_design <- function(N1, N2, r, frac, omega, n.interim, alpha, Test, alternative,
                         tsmethod, ss.method, nsim) {
  a <- check_strat_args(N1, N2, alpha, Test, alternative, tsmethod)
  J <- length(a$N1)
  fr <- check_strat_frac(if (is.null(frac)) a$N1 + a$N2 else frac, r, J)
  if (any(a$N1 != ceil_tol(fr$r * a$N2))) {
    warning('N1 differs from ceiling(r * N2) in some stratum, so the planned design is ',
            'not allocated in the ratio r to 1')
  }
  ni <- strat_n1(a$N2, fr$r, omega, n.interim)
  if (J > 2) {
    if (ss.method != 'standard') {
      stop("with three or more strata the re-estimation requires ss.method = 'standard'")
    }
    if (length(nsim) != 1 || !is.numeric(nsim) || is.na(nsim) || nsim < 100 ||
        nsim != round(nsim)) {
      stop('nsim must be a single integer of at least 100')
    }
  }
  list(a = a, J = J, frac = fr$frac, r = fr$r, n1 = ni$n1, n2 = ni$n2,
       nsim = as.integer(nsim))
}
