#' Sample Size Calculation for Two-Arm Trials with Binary Endpoints
#'
#' Calculates the required sample size for two-arm trials with binary endpoints using
#' exact statistical tests. Five tests are supported, each of which can be applied with a
#' one-sided or a two-sided alternative.
#'
#' @param p1 True probability of responders for group 1
#' @param p2 True probability of responders for group 2
#' @param r Allocation ratio to group 1 (i.e., allocation ratio of group 1:group 2 = r:1,
#'   r > 0)
#' @param alpha Level of significance for the alternative specified by \code{alternative}
#' @param tar.power Target power
#' @param Test Type of statistical test. Options: \code{'Chisq'}, \code{'Fisher'},
#'   \code{'Fisher-midP'}, \code{'Z-pool'}, or \code{'Boschloo'}
#' @param alternative Direction of the alternative hypothesis. Options: \code{'greater'}
#'   (default) or \code{'two.sided'}
#' @param tsmethod Convention used to construct the two-sided version of the conditional
#'   tests. Options: \code{'minlike'} (default) or \code{'central'}
#' @param n.grid Number of grid points used to search over the nuisance parameter of the
#'   unconditional tests. Default is 100
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure. The default of
#'   0 disables the procedure
#'
#' @return An object of class \code{bbssr_samplesize}, a data frame with one row
#'   containing:
#' \describe{
#'   \item{p1}{True probability of responders for group 1}
#'   \item{p2}{True probability of responders for group 2}
#'   \item{r}{Allocation ratio to group 1}
#'   \item{alpha}{Level of significance}
#'   \item{tar.power}{Target power}
#'   \item{Test}{Name of the statistical test}
#'   \item{alternative}{Direction of the alternative hypothesis}
#'   \item{Power}{Exact power at the selected sample size}
#'   \item{N1}{Required sample size of group 1}
#'   \item{N2}{Required sample size of group 2}
#'   \item{N}{Total required sample size}
#' }
#'
#' @details
#' The calculation uses a three-step approach:
#' \enumerate{
#'   \item Calculate an initial sample size from the normal approximation to the
#'     chi-squared test
#'   \item Evaluate the exact power at the initial sample size
#'   \item Move the sample size up or down one unit at a time until the smallest sample
#'     size attaining the target power is found
#' }
#'
#' The normal approximation of the first step uses \code{alpha} for a one-sided
#' alternative and \code{alpha / 2} for a two-sided alternative. Only the starting value of
#' the search is affected, so the returned sample size is exact in either case.
#'
#' @examples
#' # One-sided chi-squared test
#' BinarySampleSize(p1 = 0.4, p2 = 0.2, r = 1, alpha = 0.025,
#'                  tar.power = 0.8, Test = 'Chisq')
#'
#' \donttest{
#' # Two-sided Fisher exact test
#' BinarySampleSize(p1 = 0.5, p2 = 0.2, r = 2, alpha = 0.05,
#'                  tar.power = 0.9, Test = 'Fisher', alternative = 'two.sided')
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @export
#' @import fpCompare
#' @importFrom stats qnorm
BinarySampleSize <- function(p1, p2, r, alpha, tar.power, Test,
                             alternative = c('greater', 'two.sided'),
                             tsmethod = c('minlike', 'central'),
                             n.grid = 100, bb.gamma = 0) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  if (length(p1) != 1 || length(p2) != 1) stop('p1 and p2 must each be a single value')
  if (length(r) != 1 || is.na(r) || r <= 0) stop('r must be a single positive value')
  if (length(tar.power) != 1 || tar.power <= 0 || tar.power >= 1) {
    stop('tar.power must be a single value in (0, 1)')
  }
  if (p1 %==% p2) stop('p1 and p2 must differ for a sample size to exist')
  power_at <- function(N2) {
    N1 <- ceiling(r * N2)
    BinaryPower(p1, p2, N1, N2, alpha, Test, alternative, tsmethod, n.grid, bb.gamma)$Power
  }
  # Step 0 (initial sample size from the normal approximation to the chi-squared test)
  alpha.eff <- if (alternative == 'two.sided') alpha / 2 else alpha
  p <- (r * p1 + p2) / (1 + r)
  init.N2 <- '*'(
    (1 + 1 / r) / ((p1 - p2) ^ 2),
    (qnorm(alpha.eff) * sqrt(p * (1 - p)) +
       qnorm(1 - tar.power) * sqrt((p1 * (1 - p1) / r + p2 * (1 - p2)) / (1 + 1 / r))) ^ 2
  )
  # Step 1 (power calculation given the initial sample size)
  N2 <- max(1, ceiling(init.N2))
  Power <- power_at(N2)
  # Step 2 (sample size calculation via a grid search algorithm)
  if (Power %>=% tar.power) {
    while ((Power %>=% tar.power) && (N2 > 1)) {
      N2 <- N2 - 1
      Power <- power_at(N2)
    }
    if (Power %<<% tar.power) N2 <- N2 + 1
  } else {
    while (Power %<<% tar.power) {
      N2 <- N2 + 1
      Power <- power_at(N2)
    }
  }
  # Step 3 (determine the final sample size)
  N2 <- as.integer(N2)
  N1 <- as.integer(ceiling(r * N2))
  N <- N1 + N2
  Power <- power_at(N2)
  out <- data.frame(
    p1 = p1, p2 = p2, r = r, alpha = alpha, tar.power = tar.power,
    Test = Test, alternative = alternative, Power = Power,
    N1 = N1, N2 = N2, N = N,
    stringsAsFactors = FALSE
  )
  attr(out, 'tsmethod') <- tsmethod
  attr(out, 'n.grid') <- n.grid
  attr(out, 'bb.gamma') <- bb.gamma
  class(out) <- c('bbssr_samplesize', 'data.frame')
  out
}
