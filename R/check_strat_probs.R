#' Validate Response Probabilities of the Strata
#'
#' Internal helper turning a vector with one probability per stratum, or a matrix with one
#' row per scenario and one column per stratum, into a matrix, after checking that the
#' probabilities lie in the unit interval.
#'
#' @param p Probabilities, a vector of length \code{J} or a matrix with \code{J} columns
#' @param J Number of strata
#' @param name Name of the argument, used in the error messages
#'
#' @return A numeric matrix with \code{J} columns
#'
#' @keywords internal
#' @noRd
check_strat_probs <- function(p, J, name) {
  if (!is.numeric(p) || length(p) == 0 || anyNA(p) || any(p < 0 | p > 1)) {
    stop(name, ' must contain probabilities in [0, 1]')
  }
  if (is.matrix(p)) {
    if (ncol(p) != J) stop(name, ' must have one column per stratum')
    return(matrix(as.double(p), ncol = J))
  }
  if (length(p) != J) {
    stop(name, ' must have one element per stratum, or be a matrix with one column per ',
         'stratum')
  }
  matrix(as.double(p), nrow = 1)
}
