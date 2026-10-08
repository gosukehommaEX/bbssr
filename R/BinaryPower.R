#' Power Calculation for Two-Arm Trials with Binary Endpoints
#'
#' Calculates the exact power for two-arm trials with binary endpoints. Seven tests are
#' supported (see \code{\link{BinaryRR}}), each of which can be applied with a one-sided or
#' a two-sided alternative, and two of them also test non-inferiority with a margin on the
#' scale of the risk difference or the risk ratio. Vectors of response probabilities are
#' accepted.
#'
#' @param p1 True probability of responders for group 1 (can be a vector)
#' @param p2 True probability of responders for group 2 (can be a vector of the same length
#'   as \code{p1})
#' @param N1 Sample size for group 1
#' @param N2 Sample size for group 2
#' @param alpha Level of significance for the alternative specified by \code{alternative}
#' @param Test Type of statistical test. Options: \code{'Chisq'}, \code{'Fisher'},
#'   \code{'Fisher-midP'}, \code{'Z-pool'}, \code{'Boschloo'}, \code{'Blackwelder'} or
#'   \code{'Farrington-Manning'}
#' @param alternative Direction of the alternative hypothesis. Options: \code{'greater'}
#'   (default), \code{'less'} or \code{'two.sided'}
#' @param tsmethod Convention used to construct the two-sided version of the conditional
#'   tests, see \code{\link{BinaryRR}}. Options: \code{'minlike'} (default),
#'   \code{'central'} or \code{'blaker'}
#' @param n.grid Number of grid points used to search over the nuisance parameter of the
#'   unconditional tests. Default is 100
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure. The default of
#'   0 disables the procedure
#' @param margin Non-inferiority margin, on the scale given by \code{margin.scale}. On the
#'   scale of the risk difference the default of 0 gives a test of superiority, and a
#'   value other than 0 tests the null hypothesis \code{p1 - p2 <= -margin} against
#'   \code{p1 - p2 > -margin} when \code{alternative} is \code{'greater'}, and
#'   \code{p1 - p2 >= margin} against \code{p1 - p2 < margin} when it is \code{'less'}. A
#'   negative value tests for superiority by more than its absolute value. On the scale
#'   of the risk ratio the margin is a positive ratio \code{R0}, and the null hypothesis
#'   is \code{p1 / p2 <= R0} against \code{p1 / p2 > R0} when \code{alternative} is
#'   \code{'greater'}, and \code{p1 / p2 >= R0} against \code{p1 / p2 < R0} when it is
#'   \code{'less'}. A margin other than 0, and any margin on the scale of the risk ratio,
#'   requires \code{Test = 'Blackwelder'} or \code{'Farrington-Manning'} and a one-sided
#'   alternative, see \code{\link{BinaryRR}}
#' @param margin.scale Scale of \code{margin}. Options: \code{'RD'} (default) for the
#'   risk difference \code{p1 - p2} or \code{'RR'} for the risk ratio \code{p1 / p2}
#' @param ref.pvalue Logical. If \code{TRUE}, the maximization over the nuisance parameter
#'   of the unconditional tests is refined between the grid points, see
#'   \code{\link{BinaryRR}}. Default is \code{FALSE}
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
                        alternative = c('greater', 'less', 'two.sided'),
                        tsmethod = c('minlike', 'central', 'blaker'),
                        n.grid = 100, bb.gamma = 0, margin = 0, ref.pvalue = FALSE,
                        margin.scale = c('RD', 'RR')) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  margin.scale <- match.arg(margin.scale)
  if (length(p1) != length(p2)) stop('p1 and p2 should be the same length')
  if (any(p1 < 0 | p1 > 1 | p2 < 0 | p2 > 1)) stop('p1 and p2 must lie in [0, 1]')
  Test <- match.arg(Test, c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo',
                            'Blackwelder', 'Farrington-Manning'))
  rr <- get_rr(N1, N2, alpha, Test, alternative, tsmethod, n.grid, bb.gamma, ref.pvalue,
               margin, margin.scale)
  N1 <- nrow(rr) - 1L
  N2 <- ncol(rr) - 1L
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
  attr(out, 'margin') <- margin
  attr(out, 'margin.scale') <- margin.scale
  class(out) <- c('bbssr_power', 'data.frame')
  out
}
