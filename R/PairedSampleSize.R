#' Sample Size for Paired Binary Data
#'
#' Calculates the number of pairs required for a test comparing two treatments on a binary
#' endpoint in a matched-pairs design to attain a target power. The number can be obtained
#' from the exact power of the selected test or from the normal approximation.
#'
#' @param p10 Probability that only treatment 1 responds in a pair
#' @param p01 Probability that only treatment 2 responds in a pair
#' @param alpha Level of significance for the alternative specified by \code{alternative}
#' @param tar.power Target power
#' @param method How the number of pairs is obtained. \code{'exact'} (default) steps from
#'   the normal approximation to the number of pairs at which the exact power of
#'   \code{Test} attains the target, see Details. \code{'standard'} rounds up the normal
#'   approximation, formula (9) of Nam (1997), which is formula (3) of Connor (1987) when
#'   the margin is 0
#' @param search How the search of \code{method = 'exact'} chooses the number of pairs.
#'   \code{'crossing'} (default) steps from the normal approximation to a number that
#'   attains the target power while the number one smaller does not. \code{'smallest'}
#'   returns the smallest number that attains the target power. \code{'stable'} returns
#'   the smallest number from which every number up to a limit attains it. Ignored by
#'   \code{method = 'standard'}
#' @param search.limit Two numbers \code{a} and \code{b} giving the limit
#'   \code{max(ceiling(a n0), ceiling(n0 + b))} of the number of pairs examined by
#'   \code{search = 'stable'}, where \code{n0} is the normal approximation rounded up.
#'   \code{a} must be at least 1 and \code{b} non-negative. Default is \code{c(2, 50)}
#' @inheritParams PairedRR
#'
#' @return An object of class \code{bbssr_psamplesize}, a data frame with one row
#'   containing:
#' \describe{
#'   \item{p10, p01}{Probabilities of the two kinds of discordant pairs}
#'   \item{alpha}{Level of significance}
#'   \item{tar.power}{Target power}
#'   \item{Test}{Name of the statistical test}
#'   \item{alternative}{Direction of the alternative hypothesis}
#'   \item{Power}{Exact power of \code{Test} with the required number of pairs}
#'   \item{N}{Required number of pairs}
#' }
#' The attributes hold the remaining settings, among them the unrounded normal
#' approximation \code{n.raw} and the limit \code{search.limit} of the search
#' \code{'stable'} (otherwise \code{NA}).
#'
#' @details
#' The normal approximation is formula (9) of Nam (1997),
#' \code{n = (z.a sqrt(v0) + z.b sqrt(v1))^2 / (d1 - d0)^2}, where \code{d1 = p10 - p01}
#' is the true difference, \code{d0 = -margin} the difference under the null hypothesis,
#' \code{v1 = p10 + p01 - d1^2} the variance under the alternative, and
#' \code{v0 = 2 p01.bar + d0 - d0^2} the variance under the null hypothesis evaluated at
#' the large sample value \code{p01.bar} of the restricted maximum likelihood estimate of
#' \code{p01}, formula (5) of that article. \code{z.a} is the upper \code{alpha} quantile
#' of the standard normal distribution, the upper \code{alpha / 2} quantile for a
#' two-sided alternative, and \code{z.b} the \code{tar.power} quantile. With
#' \code{margin = 0} the formula is formula (3) of Connor (1987). For
#' \code{alternative = 'less'} the probabilities \code{p10} and \code{p01} are exchanged.
#' For \code{alternative = 'equivalence'} the approximate power of each one-sided test is
#' obtained in the same way for its own margin, at the level \code{alpha}, and the normal
#' approximation is the smallest number of pairs at which the sum of the two approximate
#' powers minus 1 attains the target power. The sum minus 1 is a lower bound of the
#' approximate probability that both tests reject.
#'
#' The normal approximation does not depend on the test, and for \code{'McNemar-cc'} it
#' leaves out the continuity correction. It is the starting point of the exact search for
#' every test. The exact power is not monotone in the number of pairs, so the three
#' searches can return different numbers, see \code{\link{BinarySampleSize}}.
#'
#' @references
#' Connor RJ (1987). Sample size for testing differences in proportions for the
#' paired-sample design. \emph{Biometrics}, 43(1), 207-211.
#'
#' Nam J (1997). Establishing equivalence of two treatments and sample size requirements in
#' matched-pairs design. \emph{Biometrics}, 53(4), 1422-1430.
#'
#' @examples
#' # Example of Connor (1987, p. 210): theta1 = 0.8, theta2 = 0.7 and psi = 0.3
#' PairedSampleSize(p10 = 0.2, p01 = 0.1, alpha = 0.05, tar.power = 0.9,
#'                  Test = 'McNemar', method = 'standard')
#'
#' \donttest{
#' # Exact number of pairs of the exact conditional test
#' PairedSampleSize(p10 = 0.2, p01 = 0.1, alpha = 0.05, tar.power = 0.9,
#'                  Test = 'McNemar-exact')
#' }
#'
#' @export
PairedSampleSize <- function(p10, p01, alpha, tar.power, Test,
                             alternative = c('greater', 'less', 'two.sided', 'equivalence'),
                             n.grid = 100, bb.gamma = 0, method = c('exact', 'standard'),
                             search = c('crossing', 'smallest', 'stable'),
                             search.limit = c(2, 50), margin = 0, ref.pvalue = FALSE) {
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  search <- match.arg(search)
  if (length(p10) != 1 || length(p01) != 1) stop('p10 and p01 must each be a single value')
  check_paired_probs(p10, p01)
  if (length(tar.power) != 1 || is.na(tar.power) || tar.power <= 0 || tar.power >= 1) {
    stop('tar.power must be a single value in (0, 1)')
  }
  a <- check_paired_args(1, alpha, Test, n.grid, bb.gamma, ref.pvalue, alternative, margin)
  margin <- a$margin
  # The Berger-Boos parameter matters only for the unconditional test; it is set to 0 for
  # the other tests so that the warning is issued once
  if (a$Test != 'McNemar-uncond') bb.gamma <- 0
  d1 <- p10 - p01
  if (alternative == 'equivalence') {
    if (!((d1 + margin[1]) %>>% 0 && (margin[2] - d1) %>>% 0)) {
      stop('p10 - p01 must lie strictly between -lower and upper, the margins of ',
           'equivalence')
    }
  } else if (alternative == 'two.sided') {
    if (d1 %==% 0) stop('p10 and p01 must differ for a sample size to exist')
  } else if (alternative == 'greater' && !((d1 + margin) %>>% 0)) {
    stop("p10 - p01 must exceed -margin when alternative is 'greater'")
  } else if (alternative == 'less' && !((margin - d1) %>>% 0)) {
    stop("p10 - p01 must fall below margin when alternative is 'less'")
  }
  n.raw <- paired_raw_n(p10, p01, alpha, tar.power, alternative, margin)
  limit <- NA_integer_
  if (method == 'standard') {
    N <- as.integer(max(1, ceil_tol(n.raw)))
  } else {
    s <- paired_ss_search(p10, p01, alpha, tar.power, a$Test, alternative, a$n.grid,
                          bb.gamma, a$ref.pvalue, margin, search, search.limit)
    N <- s$N
    limit <- s$limit
  }
  rr <- get_paired_rr(N, alpha, a$Test, alternative, a$n.grid, bb.gamma, a$ref.pvalue,
                      margin)
  Power <- paired_power(rr, p10, p01)
  out <- data.frame(p10 = p10, p01 = p01, alpha = alpha, tar.power = tar.power,
                    Test = a$Test, alternative = alternative, Power = Power, N = N,
                    stringsAsFactors = FALSE)
  attr(out, 'n.grid') <- a$n.grid
  attr(out, 'bb.gamma') <- bb.gamma
  attr(out, 'method') <- method
  attr(out, 'search') <- search
  attr(out, 'search.limit') <- limit
  attr(out, 'n.raw') <- n.raw
  attr(out, 'margin') <- margin
  attr(out, 'ref.pvalue') <- a$ref.pvalue
  class(out) <- c('bbssr_psamplesize', 'data.frame')
  out
}
