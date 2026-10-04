#' Unrounded Sample Size of Group 2 from the Normal Approximation
#'
#' Internal helper returning the unrounded size of group 2 required by the normal
#' approximation to the pooled Z test, which is equivalent to the Pearson chi-squared test.
#' Group 1 receives \code{r} times as many patients. The function is vectorized over
#' \code{p1} and \code{p2}.
#'
#' With \code{p = (r p1 + p2) / (1 + r)} and \code{d = p1 - p2}, the size is
#' \code{(z.a sqrt(v0) + z.b sqrt(v1))^2 / d^2}, where
#' \code{v0 = p (1 - p) (1 + 1 / r)} is the variance under the null hypothesis and
#' \code{z.a} and \code{z.b} are standard normal quantiles of the significance level and
#' the target power. Under \code{'standard'}, \code{v1 = p1 (1 - p1) / r + p2 (1 - p2)} is
#' the variance under the alternative, as in formula (21.3) of Kieser (2020). Under
#' \code{'null.variance'}, \code{v1 = v0}, as in formula (1) of Friede and Kieser (2004)
#' and formula (21.9) of Kieser (2020).
#'
#' The probabilities recovered from a pooled rate with an assumed risk difference can fall
#' outside the unit interval. Each Bernoulli variance is truncated at zero, so that such
#' probabilities remain usable and the assumed difference \code{d} is kept.
#'
#' @param p1 Response probability of group 1
#' @param p2 Response probability of group 2
#' @param r Allocation ratio to group 1
#' @param alpha Level of significance for the alternative specified by \code{alternative}
#' @param tar.power Target power
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}. A two-sided
#'   alternative uses \code{alpha / 2}
#' @param ss.method \code{'standard'} or \code{'null.variance'}
#'
#' @return A numeric vector. The value is \code{Inf} or \code{NaN} when \code{p1} equals
#'   \code{p2}
#'
#' @keywords internal
#' @noRd
#' @importFrom stats qnorm
ss_raw_n2 <- function(p1, p2, r, alpha, tar.power, alternative, ss.method) {
  alpha.eff <- if (alternative == 'two.sided') alpha / 2 else alpha
  z.a <- qnorm(1 - alpha.eff)
  z.b <- qnorm(tar.power)
  p <- (r * p1 + p2) / (1 + r)
  d <- p1 - p2
  v0 <- pmax(p * (1 - p), 0) * (1 + 1 / r)
  v1 <- if (ss.method == 'standard') {
    pmax(p1 * (1 - p1), 0) / r + pmax(p2 * (1 - p2), 0)
  } else {
    v0
  }
  (z.a * sqrt(v0) + z.b * sqrt(v1))^2 / d^2
}
