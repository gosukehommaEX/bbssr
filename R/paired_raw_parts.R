#' Ingredients of the Normal Approximation for Paired Binary Data
#'
#' Internal helper returning the distance of the true difference from the null boundary
#' and the variances under the null and the alternative hypotheses that enter the sample
#' size formula (9) of Nam (1997) for the score test of \code{\link{PairedRR}}. For
#' \code{alternative = 'greater'} the null difference is \code{d0 = -margin} and the true
#' difference \code{d1 = p10 - p01}. The asymptotic limit \code{p01.bar} of the restricted
#' maximum likelihood estimate of \code{p01} is the larger root of
#' \code{2 x^2 + b x + c = 0} with \code{b = (2 - d1) (1 + d0) - 2 (1 + p01)} and
#' \code{c = -d0 (1 - d0) p01}, formula (5) of Nam (1997). The variance under the null
#' hypothesis is \code{2 p01.bar + d0 - d0^2} and that under the alternative
#' \code{p10 + p01 - d1^2}. With \code{margin = 0}, \code{p01.bar = (p10 + p01) / 2}, the
#' variances become \code{p10 + p01} and \code{p10 + p01 - d1^2}, and the formula is
#' formula (3) of Connor (1987) for the McNemar test. The alternative \code{'less'} is
#' obtained by exchanging \code{p10} and \code{p01}.
#'
#' @param p10 Probability that only treatment 1 responds
#' @param p01 Probability that only treatment 2 responds
#' @param alternative \code{'greater'} or \code{'less'}
#' @param margin Margin, so that the null hypothesis is \code{p10 - p01 <= -margin} for
#'   \code{'greater'} and \code{p10 - p01 >= margin} for \code{'less'}
#'
#' @return A list with the distance \code{d = d1 - d0} and the variances \code{v0} and
#'   \code{v1}, vectors of the length of \code{p10}
#'
#' @references
#' Connor RJ (1987). Sample size for testing differences in proportions for the
#' paired-sample design. \emph{Biometrics}, 43(1), 207-211.
#'
#' Nam J (1997). Establishing equivalence of two treatments and sample size requirements in
#' matched-pairs design. \emph{Biometrics}, 53(4), 1422-1430.
#'
#' @keywords internal
#' @noRd
paired_raw_parts <- function(p10, p01, alternative, margin) {
  if (alternative == 'less') {
    tmp <- p10
    p10 <- p01
    p01 <- tmp
  }
  d0 <- -margin
  d1 <- p10 - p01
  b <- (2 - d1) * (1 + d0) - 2 * (1 + p01)
  c <- -d0 * (1 - d0) * p01
  p01.bar <- (sqrt(pmax(0, b ^ 2 - 8 * c)) - b) / 4
  list(d = d1 - d0, v0 = pmax(0, 2 * p01.bar + d0 - d0 ^ 2),
       v1 = pmax(0, p10 + p01 - d1 ^ 2))
}
