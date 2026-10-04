#' Maximum of a Smooth Function from a Grid and Local Refinement
#'
#' Internal helper locating the maximum of a function of one variable that has been
#' evaluated on a grid. The local maxima of the grid values with the largest values are
#' refined by \code{stats::optimize} between the neighbouring grid points, and the largest
#' value found is returned. Rejection probabilities are polynomials in the response
#' probability, so the refinement recovers maxima that fall between grid points.
#'
#' @param f Function of one numeric argument, returning one numeric value
#' @param x Increasing grid of arguments
#' @param y Values of \code{f} at \code{x}
#' @param n.peaks Number of local maxima of the grid values that are refined
#' @param tol Tolerance passed to \code{stats::optimize}
#'
#' @return A list with components \code{x}, the location of the maximum, and \code{y},
#'   its value
#'
#' @keywords internal
#' @noRd
#' @importFrom stats optimize
refine_max <- function(f, x, y, n.peaks = 3, tol = 1e-10) {
  n <- length(x)
  best <- list(x = x[which.max(y)], y = max(y))
  if (n < 2) return(best)
  peak <- which(y >= c(-Inf, y[-n]) & y >= c(y[-1], -Inf))
  peak <- peak[order(y[peak], decreasing = TRUE)][seq_len(min(n.peaks, length(peak)))]
  for (i in peak) {
    lo <- x[max(1, i - 1)]
    hi <- x[min(n, i + 1)]
    o <- optimize(f, c(lo, hi), maximum = TRUE, tol = tol)
    if (o$objective > best$y) best <- list(x = o$maximum, y = o$objective)
  }
  best
}
