#' Sample Size of a Stratified Test
#'
#' Calculates the total sample size of a two-arm trial with a binary endpoint stratified
#' into two or more strata, for a test of \code{\link{StratRR}}, either exactly or by the
#' normal approximation of Woolson et al. (1986) and Nam (1992).
#'
#' @param p1 Response probabilities of group 1 in the strata, one element per stratum
#' @param p2 Response probabilities of group 2 in the strata, one element per stratum
#' @param frac Fractions of the strata, one positive element per stratum. They are rescaled
#'   to sum to 1, so they may also be given as numbers of patients
#' @param r Allocation ratio of group 1 to group 2 within the strata, a single value or one
#'   value per stratum. Default is 1
#' @param tar.power Target power
#' @param method Calculation of the sample size. \code{'exact'} (default) searches for the
#'   total sample size at which the exact power of \code{\link{StratPower}} attains
#'   \code{tar.power}. \code{'standard'} uses the normal approximation, see Details
#' @param search How the exact calculation chooses among the sizes at which the power
#'   fluctuates around \code{tar.power}, as in \code{\link{BinarySampleSize}}. Options:
#'   \code{'crossing'} (default), \code{'smallest'} or \code{'stable'}
#' @param search.limit Limit of the search \code{'stable'}: the nominal sizes are examined
#'   downwards from \code{max(ceiling(search.limit[1] n), ceiling(n + search.limit[2]))},
#'   where \code{n} is the nominal size of the normal approximation. Default is
#'   \code{c(2, 50)}
#' @inheritParams StratRR
#'
#' @return An object of class \code{bbssr_ssamplesize}, a list containing the arguments
#'   \code{p1}, \code{p2}, \code{frac}, \code{r}, \code{alpha}, \code{tar.power},
#'   \code{Test} and \code{alternative}, the group sizes \code{N1} and \code{N2} of the
#'   strata, their total \code{N}, and the exact power \code{Power} of the design. The
#'   attribute \code{n.raw} holds the unrounded total of the normal approximation and
#'   \code{n} the nominal total of the design.
#'
#' @details
#' A design is described by the fractions \code{frac} of the strata and the allocation
#' ratios \code{r}. For a nominal total sample size \code{n}, the size of group 2 in
#' stratum \code{j} is \code{ceiling(frac[j] n / (1 + r[j]))} and that of group 1 is
#' \code{ceiling(r[j] N2[j])}, so the total \code{N} is at least \code{n}. The exact
#' calculation searches over \code{n} within this family of designs.
#'
#' The normal approximation is that of Cochran's test. With
#' \code{s[j] = r[j] / (1 + r[j])}, \code{w[j] = frac[j] s[j] (1 - s[j])},
#' \code{D = sum_j w[j] (p1[j] - p2[j])},
#' \code{v0 = sum_j w[j] pbar[j] (1 - pbar[j])} with
#' \code{pbar[j] = s[j] p1[j] + (1 - s[j]) p2[j]}, and
#' \code{v1 = sum_j w[j] ((1 - s[j]) p1[j] (1 - p1[j]) + s[j] p2[j] (1 - p2[j]))}, the
#' total sample size is \code{N* = (z.a sqrt(v0) + z.b sqrt(v1))^2 / D^2}, formula (2.8) of
#' Woolson et al. (1986), where \code{z.a} is the upper \code{alpha} quantile
#' (\code{alpha / 2} for a two-sided alternative) and \code{z.b} the \code{tar.power}
#' quantile of the standard normal distribution. For the tests with a continuity
#' correction it is \code{N = (1 + sqrt(1 + 2 / (D N*)))^2 N* / 4}, formula (2.3) of Nam
#' (1992). The same approximation is used for every test, since the
#' Cochran-Mantel-Haenszel statistic differs from Cochran's only by the factors
#' \code{(N_j - 1) / N_j} of the variance, and it serves as the starting value of the exact
#' search. The nominal total is \code{ceiling(N*)}, or \code{ceiling(N)}; the published
#' tables of Woolson et al. (1986) and Nam (1992) round \code{N*} and \code{N} to the
#' nearest integer instead.
#'
#' The exact power is not monotone in the sample size, so the result depends on
#' \code{search}, as described in \code{\link{BinarySampleSize}}.
#'
#' @references
#' Woolson RF, Bean JA, Rojas PB (1986). Sample size for case-control studies using
#' Cochran's statistic. \emph{Biometrics}, 42(4), 927-932.
#'
#' Nam J (1992). Sample size determination for case-control studies and the comparison of
#' stratified and unstratified analyses. \emph{Biometrics}, 48(2), 389-395.
#'
#' Jung SH (2014). Stratified Fisher's exact test and its sample size calculation.
#' \emph{Biometrical Journal}, 56(1), 129-140.
#'
#' @examples
#' # Normal approximation and exact sample size with two strata
#' StratSampleSize(p1 = c(0.5, 0.7), p2 = c(0.3, 0.5), frac = c(0.5, 0.5),
#'                 alpha = 0.025, tar.power = 0.8, Test = 'CMH', method = 'standard')
#' StratSampleSize(p1 = c(0.5, 0.7), p2 = c(0.3, 0.5), frac = c(0.5, 0.5),
#'                 alpha = 0.025, tar.power = 0.8, Test = 'CMH')
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{StratPower}}, \code{\link{StratPowerBSSR}}
#' @export
StratSampleSize <- function(p1, p2, frac, r = 1, alpha, tar.power, Test,
                            alternative = c('greater', 'less', 'two.sided'),
                            tsmethod = c('minlike', 'central', 'blaker'),
                            method = c('exact', 'standard'),
                            search = c('crossing', 'smallest', 'stable'),
                            search.limit = c(2, 50)) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  method <- match.arg(method)
  search <- match.arg(search)
  J <- length(frac)
  fr <- check_strat_frac(frac, r, J)
  a <- check_strat_args(rep(1, J), rep(1, J), alpha, Test, alternative, tsmethod)
  q1 <- check_strat_probs(p1, J, 'p1')
  q2 <- check_strat_probs(p2, J, 'p2')
  if (nrow(q1) != 1 || nrow(q2) != 1) {
    stop('p1 and p2 must each have one element per stratum')
  }
  if (length(tar.power) != 1 || !is.numeric(tar.power) || is.na(tar.power) ||
      tar.power <= 0 || tar.power >= 1) {
    stop('tar.power must be a single value in (0, 1)')
  }
  n.raw <- strat_raw_n(q1, q2, fr$frac, fr$r, alpha, tar.power, a$alternative, a$cc)
  if (!is.finite(n.raw)) {
    stop('the response probabilities show no effect in the direction of the ',
         'alternative: sum_j w[j] (p1[j] - p2[j]) must be ',
         switch(a$alternative, greater = 'positive', less = 'negative',
                two.sided = 'non-zero'))
  }
  limit <- NA_integer_
  if (method == 'standard') {
    n <- as.integer(max(1, ceil_tol(n.raw)))
  } else {
    s <- strat_ss_search(q1, q2, fr$frac, fr$r, alpha, tar.power, a$Test, a$alternative,
                         a$tsmethod, search, search.limit)
    n <- s$n
    limit <- s$limit
  }
  sz <- strat_sizes(n, fr$frac, fr$r)
  N1 <- as.integer(sz$N1[1, ])
  N2 <- as.integer(sz$N2[1, ])
  Power <- strat_fixed_power(N1, N2, q1, q2, alpha, a$Test, a$alternative, a$tsmethod)
  out <- list(p1 = as.vector(q1), p2 = as.vector(q2), frac = fr$frac, r = fr$r,
              alpha = alpha, tar.power = tar.power, Test = a$Test,
              alternative = a$alternative, N1 = N1, N2 = N2,
              N = as.integer(sum(N1) + sum(N2)), Power = Power)
  attr(out, 'n.raw') <- n.raw
  attr(out, 'n') <- n
  attr(out, 'method') <- method
  attr(out, 'search') <- search
  attr(out, 'search.limit') <- limit
  attr(out, 'tsmethod') <- a$tsmethod
  class(out) <- 'bbssr_ssamplesize'
  out
}
