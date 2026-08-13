#' Power Calculation for Two-Arm Trials with Binary Endpoints
#'
#' Calculates power for two-arm trials with binary endpoints using exact statistical tests.
#' Five tests are supported, each of which can be applied with a one-sided or a two-sided
#' alternative, and vectors of response probabilities are accepted.
#'
#' @param p1 True probability of responders for group 1 (can be a vector)
#' @param p2 True probability of responders for group 2 (can be a vector of the same length
#'   as \code{p1})
#' @param N1 Sample size for group 1
#' @param N2 Sample size for group 2
#' @param alpha Level of significance for the alternative specified by \code{alternative}
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
#' @return An object of class \code{bbssr_power}, a data frame with one row per element of
#'   \code{p1} containing:
#' \describe{
#'   \item{p1}{True probability of responders for group 1}
#'   \item{p2}{True probability of responders for group 2}
#'   \item{N1}{Sample size for group 1}
#'   \item{N2}{Sample size for group 2}
#'   \item{alpha}{Level of significance}
#'   \item{Test}{Name of the statistical test}
#'   \item{alternative}{Direction of the alternative hypothesis}
#'   \item{Power}{Exact power}
#' }
#'
#' @details
#' The power is obtained by summing the joint probability mass function of the two
#' independent binomial counts over the rejection region returned by \code{\link{BinaryRR}}.
#' The summation covers the whole rejection region rather than a row-wise tail, so it
#' remains valid for two-sided tests, whose rejection regions are not contiguous within a
#' row of the outcome grid.
#'
#' @examples
#' # Power of the one-sided chi-squared test
#' BinaryPower(p1 = 0.5, p2 = 0.2, N1 = 5, N2 = 5, alpha = 0.025, Test = 'Chisq')
#'
#' \donttest{
#' # Power over a range of response probabilities for the two-sided Boschloo test
#' pw <- BinaryPower(p1 = c(0.5, 0.6, 0.7, 0.8), p2 = rep(0.2, 4),
#'                   N1 = 20, N2 = 20, alpha = 0.05, Test = 'Boschloo',
#'                   alternative = 'two.sided')
#' print(pw)
#' plot(pw)
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @export
#' @importFrom stats dbinom
BinaryPower <- function(p1, p2, N1, N2, alpha, Test,
                        alternative = c('greater', 'two.sided'),
                        tsmethod = c('minlike', 'central'),
                        n.grid = 100, bb.gamma = 0) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  if (length(p1) != length(p2)) stop('p1 and p2 should be the same length')
  if (any(p1 < 0 | p1 > 1 | p2 < 0 | p2 > 1)) stop('p1 and p2 must lie in [0, 1]')
  RR <- BinaryRR(N1, N2, alpha, Test, alternative, tsmethod, n.grid, bb.gamma)
  Test <- attr(RR, 'Test')
  N1 <- attr(RR, 'N1')
  N2 <- attr(RR, 'N2')
  rr <- matrix(as.vector(RR), nrow = N1 + 1L, ncol = N2 + 1L)
  Power <- vapply(
    seq_along(p1),
    function(i) power_from_rr(rr, dbinom(0:N1, N1, p1[i]), dbinom(0:N2, N2, p2[i])),
    numeric(1)
  )
  out <- data.frame(
    p1 = p1, p2 = p2, N1 = N1, N2 = N2, alpha = alpha,
    Test = Test, alternative = alternative, Power = Power,
    stringsAsFactors = FALSE
  )
  class(out) <- c('bbssr_power', 'data.frame')
  out
}
