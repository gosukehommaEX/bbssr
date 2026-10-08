#' Distance of Response Probabilities from the Null Boundary
#'
#' Internal helper returning the signed distance of the response probabilities from the
#' boundary of the null hypothesis, positive on the side of the alternative. On the scale
#' of the risk difference the distance is \code{p1 - p2 + margin} for
#' \code{alternative = 'greater'} and \code{p2 - p1 + margin} for \code{'less'}, and on
#' the scale of the risk ratio it is \code{p1 - margin p2} for \code{'greater'} and
#' \code{margin p2 - p1} for \code{'less'}. These are the numerators of the statistics of

#' Farrington and Manning (1990) evaluated at the true probabilities, and the differences
#' that enter their sample size formulae. A two-sided alternative is treated as
#' \code{'greater'}. The function is vectorized over \code{p1} and \code{p2}.
#'
#' @param p1 Response probability of group 1
#' @param p2 Response probability of group 2
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param margin Non-inferiority margin, a difference for \code{margin.scale = 'RD'} and
#'   a ratio for \code{'RR'}
#' @param margin.scale \code{'RD'} or \code{'RR'}
#'
#' @return A numeric vector
#'
#' @keywords internal
#' @noRd
null_distance <- function(p1, p2, alternative, margin, margin.scale) {
  if (margin.scale == 'RR') {
    if (alternative == 'less') margin * p2 - p1 else p1 - margin * p2
  } else if (alternative == 'less') {
    p2 - p1 + margin
  } else {
    p1 - p2 + margin
  }
}
