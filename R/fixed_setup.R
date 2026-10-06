#' Interim Setup of a Fixed-Sample Design
#'
#' Internal helper describing a fixed-sample design in the form returned by
#' \code{bssr_setup}: a single interim cell without patients, from which the sample sizes
#' \code{N1} and \code{N2} are always reached. With this setup the helpers written for
#' re-estimation designs, such as \code{tie_weights}, apply to the fixed-sample design.
#'
#' @param N1 Sample size of group 1
#' @param N2 Sample size of group 2
#'
#' @return A list with the components of the result of \code{bssr_setup}
#'
#' @keywords internal
#' @noRd
fixed_setup <- function(N1, N2) {
  list(n11 = 0L, n12 = 0L, x11 = 0L, x12 = 0L, rr.id = 0L, N1 = N1, N2 = N2,
       n21 = as.integer(N1), n22 = as.integer(N2))
}
