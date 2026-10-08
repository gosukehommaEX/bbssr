#' Unrounded Sample Size of Group 2 from the Normal Approximation
#'
#' Internal helper returning the unrounded size of group 2 required by the normal
#' approximation to the pooled Z test, which is equivalent to the Pearson chi-squared test,
#' or by the normal approximation to a test of non-inferiority. Group 1 receives \code{r}
#' times as many patients. The function is vectorized over \code{p1} and \code{p2}.
#'
#' With \code{p = (r p1 + p2) / (1 + r)} and \code{d = p1 - p2}, the size is
#' \code{(z.a sqrt(v0) + z.b sqrt(v1))^2 / d^2}, where
#' \code{v0 = p (1 - p) (1 + 1 / r)} is the variance under the null hypothesis and
#' \code{z.a} and \code{z.b} are standard normal quantiles of the significance level and
#' the target power. Under \code{'standard'}, \code{v1 = p1 (1 - p1) / r + p2 (1 - p2)} is
#' the variance under the alternative, as in formula (21.3) of Kieser (2020). Under
#' \code{'null.variance'}, \code{v1 = v0}, as in formula (1) of Friede and Kieser (2004)
#' and formula (21.9) of Kieser (2020). Under \code{'alternative.variance'}, \code{v0} is
#' replaced by \code{v1}, as in Blackwelder (1982).
#'
#' With a non-zero \code{margin} and \code{alternative = 'greater'}, the difference is
#' \code{d = p1 - p2 + margin} and \code{v0} is the variance at the large sample values of
#' the restricted estimates of Farrington and Manning (1990) under
#' \code{p1 - p2 = -margin}. The method \code{'standard'} is then formula (4) of
#' Farrington and Manning (1990) and formula (2) of Friede et al. (2007), and
#' \code{'alternative.variance'} is the formula of Blackwelder (1982) and formula (1) of
#' Friede et al. (2007). With
#' \code{alternative = 'less'} the two groups exchange their roles. The restricted
#' estimates are computed from the probabilities truncated to the unit interval.
#'
#' With \code{margin.scale = 'RR'} the margin is the ratio \code{R0} of the boundary
#' \code{p1 = R0 p2}, the difference is \code{d = p1 - R0 p2}, or \code{R0 p2 - p1} for
#' \code{alternative = 'less'}, and the variances are \code{v1 = p1 (1 - p1) / r +
#' R0^2 p2 (1 - p2)} and \code{v0}, the same expression at the large sample values of the
#' restricted estimates of \code{fm_restricted_rr}. The method \code{'standard'} is then
#' formula (8) of Farrington and Manning (1990) divided by \code{r}, which gives the size
#' of group 2, and \code{'alternative.variance'} is Method 1 of that article.
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
#' @param ss.method \code{'standard'}, \code{'null.variance'} or
#'   \code{'alternative.variance'}
#' @param margin Non-inferiority margin, 0 for a test of superiority on the scale of the
#'   risk difference
#' @param margin.scale \code{'RD'} or \code{'RR'}, the scale of the margin
#'
#' @return A numeric vector. The value is \code{Inf} or \code{NaN} when \code{p1} and
#'   \code{p2} lie on the boundary of the null hypothesis
#'
#' @keywords internal
#' @noRd
#' @importFrom stats qnorm
ss_raw_n2 <- function(p1, p2, r, alpha, tar.power, alternative, ss.method, margin,
                      margin.scale) {
  alpha.eff <- if (alternative == 'two.sided') alpha / 2 else alpha
  z.a <- qnorm(1 - alpha.eff)
  z.b <- qnorm(tar.power)
  p <- (r * p1 + p2) / (1 + r)
  d <- p1 - p2
  v0 <- pmax(p * (1 - p), 0) * (1 + 1 / r)
  if (margin.scale == 'RR') {
    q1 <- pmin(1, pmax(0, p1))
    q2 <- pmin(1, pmax(0, p2))
    d <- null_distance(p1, p2, alternative, margin, 'RR')
    est <- fm_restricted_rr(q1, q2, 1 / r, margin)
    v0 <- pmax(est$p1 * (1 - est$p1), 0) / r + margin^2 * pmax(est$p2 * (1 - est$p2), 0)
    v.alt <- pmax(p1 * (1 - p1), 0) / r + margin^2 * pmax(p2 * (1 - p2), 0)
    v1 <- if (ss.method == 'null.variance') v0 else v.alt
    if (ss.method == 'alternative.variance') v0 <- v.alt
    return((z.a * sqrt(v0) + z.b * sqrt(v1))^2 / d^2)
  }
  if (margin != 0) {
    q1 <- pmin(1, pmax(0, p1))
    q2 <- pmin(1, pmax(0, p2))
    if (alternative == 'less') {
      d <- p2 - p1 + margin
      est <- fm_restricted(q2, q1, r, -margin)
      v0 <- pmax(est$p2 * (1 - est$p2), 0) / r + pmax(est$p1 * (1 - est$p1), 0)
    } else {
      d <- p1 - p2 + margin
      est <- fm_restricted(q1, q2, 1 / r, -margin)
      v0 <- pmax(est$p1 * (1 - est$p1), 0) / r + pmax(est$p2 * (1 - est$p2), 0)
    }
  }
  v.alt <- pmax(p1 * (1 - p1), 0) / r + pmax(p2 * (1 - p2), 0)
  v1 <- if (ss.method == 'null.variance') v0 else v.alt
  if (ss.method == 'alternative.variance') v0 <- v.alt
  (z.a * sqrt(v0) + z.b * sqrt(v1))^2 / d^2
}
