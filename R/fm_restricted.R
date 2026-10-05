#' Restricted Maximum Likelihood Estimates of Farrington and Manning
#'
#' Internal helper returning the maximum likelihood estimates of two response
#' probabilities under the restriction \code{p1 - p2 = s0}, formula (12) of Farrington and
#' Manning (1990). The estimate of \code{p1} is the root in the admissible range of a cubic
#' equation, given in closed form. With observed proportions as input the function gives
#' the estimates used by the test statistic, and with the response probabilities expected
#' in the trial it gives the large sample values used by the sample size formula. The
#' estimates are moved onto the range allowed by the restriction, which absorbs rounding
#' error.
#'
#' @param p1 Proportions of responders of group 1, in the unit interval
#' @param p2 Proportions of responders of group 2, of the same length as \code{p1}
#' @param theta Ratio \code{N2 / N1} of the group sizes
#' @param s0 Difference \code{p1 - p2} under the null hypothesis, a single value in
#'   \code{(-1, 1)}
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
fm_restricted <- function(p1, p2, theta, s0) {
  p1 <- as.vector(p1)
  p2 <- as.vector(p2)
  a <- 1 + theta
  b <- -(1 + theta + p1 + theta * p2 + s0 * (theta + 2))
  c <- s0^2 + s0 * (2 * p1 + theta + 1) + p1 + theta * p2
  d <- -p1 * s0 * (1 + s0)
  v <- b^3 / (3 * a)^3 - b * c / (6 * a^2) + d / (2 * a)
  u <- sign(v) * sqrt(pmax(b^2 / (3 * a)^2 - c / (3 * a), 0))
  # When u is 0, w is pi / 2 for any value of the ratio and the cosine term vanishes
  ratio <- ifelse(u == 0, 0, v / u^3)
  w <- (pi + acos(pmin(1, pmax(-1, ratio)))) / 3
  t1 <- 2 * u * cos(w) - b / (3 * a)
  t1 <- pmin(min(1, 1 + s0), pmax(max(0, s0), t1))
  list(p1 = t1, p2 = t1 - s0)
}
