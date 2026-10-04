#' Ceiling with a Tolerance for Rounding Errors
#'
#' Internal helper returning the smallest integer not less than \code{x - tol}. A value
#' that exceeds an integer only through floating-point error, such as
#' \code{0.1 * 3 * 10}, is therefore not rounded up to the next integer. It is used for
#' sample sizes obtained from closed-form formulas.
#'
#' @param x Numeric vector
#' @param tol Absolute tolerance. Default is \code{1e-9}
#'
#' @return A numeric vector of whole numbers (\code{Inf} is returned unchanged)
#'
#' @keywords internal
#' @noRd
ceil_tol <- function(x, tol = 1e-9) {
  ceiling(x - tol)
}
