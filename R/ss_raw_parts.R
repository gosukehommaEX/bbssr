#' Ingredients of the Normal Approximation to the Sample Size
#'
#' Internal helper returning the difference \code{d} and the variances \code{v0} under
#' the null hypothesis and \code{v1} under the alternative, per unit of the size of
#' group 2, that enter the sample size \code{(z.a sqrt(v0) + z.b sqrt(v1))^2 / d^2} of
#' \code{ss_raw_n2}. The definitions under each \code{ss.method}, margin and scale are
#' those given in \code{ss_raw_n2}. The alternative is \code{'greater'}, \code{'less'} or
#' \code{'two.sided'}; a two-sided alternative gives the same ingredients as
#' \code{'greater'}. The function is vectorized over \code{p1} and \code{p2}.
#'
#' @param p1 Response probability of group 1
#' @param p2 Response probability of group 2
#' @param r Allocation ratio to group 1
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param ss.method \code{'standard'}, \code{'null.variance'} or
#'   \code{'alternative.variance'}
#' @param margin Single margin, 0 for a test of superiority on the scale of the risk
#'   difference
#' @param margin.scale \code{'RD'} or \code{'RR'}, the scale of the margin
#'
#' @return A list with the numeric vectors \code{d}, \code{v0} and \code{v1}
#'
#' @keywords internal
#' @noRd
ss_raw_parts <- function(p1, p2, r, alternative, ss.method, margin, margin.scale) {
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
  } else {
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
  }
  v1 <- if (ss.method == 'null.variance') v0 else v.alt
  if (ss.method == 'alternative.variance') v0 <- v.alt
  list(d = d, v0 = v0, v1 = v1)
}
