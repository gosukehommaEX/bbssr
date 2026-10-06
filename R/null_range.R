#' Interval of Pooled Response Probabilities on the Null Boundary
#'
#' Internal helper returning the interval over which the largest type I error rate is
#' certified: the interval from the smallest to the largest element of \code{theta},
#' restricted to the pooled response probabilities at which both response probabilities
#' on the boundary of the null hypothesis lie in the unit interval. On the boundary,
#' \code{p1 = theta + d0 / (1 + r)} and \code{p2 = theta - r d0 / (1 + r)}, with \code{d0}
#' as in \code{null_boundary}, so the restriction is itself an interval. Its ends are
#' widened to the elements of \code{theta} that \code{null_boundary} accepts, which can
#' lie outside it by rounding error.
#'
#' @param theta Pooled response probabilities, at least one of which is accepted by
#'   \code{null_boundary}
#' @param r Allocation ratio to group 1
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param margin Non-inferiority margin on the scale of the risk difference
#'
#' @return A numeric vector with the two ends of the interval
#'
#' @keywords internal
#' @noRd
null_range <- function(theta, r, alternative, margin) {
  d0 <- if (alternative == 'less') margin else -margin
  lo <- max(min(theta), -d0 / (1 + r), r * d0 / (1 + r))
  hi <- min(max(theta), 1 - d0 / (1 + r), 1 + r * d0 / (1 + r))
  ok <- theta[null_boundary(theta, r, alternative, margin)$ok]
  c(min(lo, ok), max(hi, ok))
}
