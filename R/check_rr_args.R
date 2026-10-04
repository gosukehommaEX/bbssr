#' Validate the Arguments of a Rejection Region
#'
#' Internal helper shared by \code{\link{BinaryRR}} and \code{get_rr}. It checks the
#' sample sizes, the level, the grid size and the Berger-Boos parameter, warns when the
#' Berger-Boos parameter is supplied for a conditional test, and returns the arguments in
#' the form used by the computation.
#'
#' @param N1 Sample size for group 1
#' @param N2 Sample size for group 2
#' @param alpha Level of significance
#' @param Test Type of statistical test, possibly abbreviated
#' @param n.grid Number of grid points over the nuisance parameter
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure
#'
#' @return A list with the integer sample sizes \code{N1} and \code{N2}, the integer grid
#'   size \code{n.grid} and the full name of the test \code{Test}
#'
#' @keywords internal
#' @noRd
check_rr_args <- function(N1, N2, alpha, Test, n.grid, bb.gamma) {
  Test <- match.arg(Test, c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo'))
  if (length(N1) != 1 || length(N2) != 1) stop('N1 and N2 must each be a single value')
  if (is.na(N1) || is.na(N2) || N1 != round(N1) || N2 != round(N2) || N1 < 1 || N2 < 1) {
    stop('N1 and N2 must be positive integers')
  }
  if (length(alpha) != 1 || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop('alpha must be a single value in (0, 1)')
  }
  if (length(n.grid) != 1 || is.na(n.grid) || n.grid < 2) {
    stop('n.grid must be a single integer of at least 2')
  }
  if (length(bb.gamma) != 1 || is.na(bb.gamma) || bb.gamma < 0 || bb.gamma >= alpha) {
    stop('bb.gamma must be a single value satisfying 0 <= bb.gamma < alpha')
  }
  if (bb.gamma > 0 && !(Test %in% c('Z-pool', 'Boschloo'))) {
    warning('bb.gamma applies only to the Z-pool and Boschloo tests and is ignored here')
  }
  list(N1 = as.integer(N1), N2 = as.integer(N2), n.grid = as.integer(n.grid), Test = Test)
}
