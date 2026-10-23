#' Points of the Null Hypothesis of a Stratified Design
#'
#' Internal helper turning the argument \code{theta} of \code{\link{StratTypeIErrorBSSR}}
#' and \code{\link{StratAlphaAdjBSSR}} into a matrix of common response probabilities of
#' the strata, one row per point and one column per stratum. A matrix is used as it is. A
#' vector gives, with two strata, every pair of its values, and with more strata the points
#' at which every stratum has the same value.
#'
#' @param theta Vector of probabilities, or matrix with one column per stratum
#' @param J Number of strata
#'
#' @return A numeric matrix with \code{J} columns named \code{theta.S1, theta.S2, ...}
#'
#' @keywords internal
#' @noRd
strat_theta <- function(theta, J) {
  if (!is.numeric(theta) || length(theta) == 0 || anyNA(theta) ||
      any(theta < 0 | theta > 1)) {
    stop('theta must contain probabilities in [0, 1]')
  }
  if (is.matrix(theta)) {
    if (ncol(theta) != J) stop('theta must have one column per stratum')
    out <- matrix(as.double(theta), ncol = J)
  } else {
    th <- sort(unique(as.double(theta)))
    out <- if (J == 2) {
      as.matrix(expand.grid(th, th))
    } else {
      matrix(th, nrow = length(th), ncol = J)
    }
  }
  dimnames(out) <- list(NULL, paste0('theta.S', seq_len(J)))
  out
}
