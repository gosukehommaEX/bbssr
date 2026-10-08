#' Restricted Maximum Likelihood Estimates under a Risk Ratio
#'
#' Internal helper returning the maximum likelihood estimates of two response
#' probabilities under the restriction \code{p1 = R0 p2}, formula (13) of Farrington and
#' Manning (1990). The estimate of \code{p1} is the smaller root of the quadratic
#' \code{a x^2 + b x + c = 0} with \code{a = 1 + theta},
#' \code{b = -(R0 (1 + theta p2) + theta + p1)} and \code{c = R0 (p1 + theta p2)}. The
#' root \code{(-b - sqrt(b^2 - 4 a c)) / (2 a)} of the article is evaluated as
#' \code{2 c / (-b + sqrt(b^2 - 4 a c))}, which is the same value without the cancellation
#' of the two terms when \code{c} is small, since \code{b} is negative. With observed
#' proportions as input the function gives the estimates used by the test statistic, and
#' with the response probabilities expected in the trial it gives the large sample values
#' used by the sample size formula. The estimates are moved onto the range allowed by the
#' restriction, which absorbs rounding error.
#'
#' @param p1 Proportions of responders of group 1, in the unit interval
#' @param p2 Proportions of responders of group 2, of the same length as \code{p1}
#' @param theta Ratio \code{N2 / N1} of the group sizes
#' @param R0 Ratio \code{p1 / p2} under the null hypothesis, a single positive value
#'
#' @return A list with numeric vectors \code{p1} and \code{p2}
#'
#' @references
#' Farrington CP, Manning G (1990). Test statistics and sample size formulae for
#' comparative binomial trials with null hypothesis of non-zero risk difference or
#' non-unity relative risk. \emph{Statistics in Medicine}, 9(12), 1447-1454.
#'
#' @keywords internal
#' @noRd
fm_restricted_rr <- function(p1, p2, theta, R0) {
  p1 <- as.vector(p1)
  p2 <- as.vector(p2)
  a <- 1 + theta
  b <- -(R0 * (1 + theta * p2) + theta + p1)
  c <- R0 * (p1 + theta * p2)
  t1 <- 2 * c / (-b + sqrt(pmax(b^2 - 4 * a * c, 0)))
  t1 <- pmin(min(1, R0), pmax(0, t1))
  list(p1 = t1, p2 = t1 / R0)
}
