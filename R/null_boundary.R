#' Response Probabilities on the Boundary of the Null Hypothesis
#'
#' Internal helper returning the response probabilities of the two groups on the boundary
#' of the null hypothesis, for given pooled response probabilities \code{theta}. The
#' boundary is \code{p1 - p2 = -margin} for \code{alternative = 'greater'} and
#' \code{p1 - p2 = margin} for \code{'less'}, and the pooled probability is
#' \code{(r p1 + p2) / (1 + r)}, so the probabilities follow from \code{split_pooled}.
#' With \code{margin = 0} both probabilities equal \code{theta}. With
#' \code{margin.scale = 'RR'} the boundary is \code{p1 = margin p2} for both one-sided
#' alternatives.
#'
#' @param theta Pooled response probabilities
#' @param r Allocation ratio to group 1
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param margin Non-inferiority margin, a difference for \code{margin.scale = 'RD'} and
#'   a ratio for \code{'RR'}
#' @param margin.scale \code{'RD'} or \code{'RR'}
#'
#' @return A list with the probabilities \code{p1} and \code{p2}, moved onto the unit
#'   interval, and the logical vector \code{ok}, which is \code{FALSE} where they lie
#'   outside the unit interval by more than rounding error
#'
#' @keywords internal
#' @noRd
#' @import fpCompare
null_boundary <- function(theta, r, alternative, margin, margin.scale) {
  if (margin.scale == 'RR') {
    sp <- split_pooled(theta, margin, r, 'RR')
  } else {
    d0 <- if (alternative == 'less') margin else -margin
    sp <- split_pooled(theta, d0, r, 'RD')
  }
  ok <- (sp$p1 %>=% 0) & (sp$p1 %<=% 1) & (sp$p2 %>=% 0) & (sp$p2 %<=% 1)
  list(p1 = pmin(1, pmax(0, sp$p1)), p2 = pmin(1, pmax(0, sp$p2)), ok = ok)
}
