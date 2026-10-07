#' Adjusted Level Reported between Two p-values
#'
#' Internal helper of \code{BinaryAlphaAdjBSSR} choosing the adjusted level to report from
#' the level found by a bisection. At a level the package rejects the p-values that are
#' below it by more than the tolerance of \code{fpCompare} (about 1.49e-8). The bisection
#' therefore ends just below the smallest p-value not rejected plus this tolerance. That
#' level rejects the smallest p-value not rejected when a p-value is rejected if it is
#' below the level or at most equal to it, and rounding it up can make it do so under the
#' rule of the package as well. The function returns instead the largest value with at
#' most \code{digits} significant digits below the smallest p-value not rejected at which
#' the rejected p-values are the same under the three rules, trying more digits, up to 15,
#' if there is none. The level found is returned unchanged if no p-value is rejected, if
#' every p-value is rejected or if no such value exists.
#'
#' @param pv.list List of matrices of p-values, one for each final sample size of the
#'   design
#' @param level Level found by the bisection
#' @param digits Smallest number of significant digits tried. Default is 6, the number of
#'   digits printed by \code{print.bbssr_alphaadj}
#'
#' @return The level to report
#'
#' @keywords internal
#' @noRd
report_level <- function(pv.list, level, digits = 6L) {
  p <- unlist(lapply(pv.list, as.vector), use.names = FALSE)
  p <- p[is.finite(p)]
  rej <- p %<<% level
  if (!any(rej) || all(rej)) return(level)
  hi <- min(p[!rej])
  # The same p-values are rejected at x under the rule of the package and when a p-value
  # is rejected if it is below x or at most x
  same <- function(x) {
    identical(p %<<% x, rej) && identical(p <= x, rej) && identical(p < x, rej)
  }
  for (d in digits:15) {
    # Largest value below hi with at most d significant digits
    u <- 10^(floor(log10(hi)) - d + 1)
    cand <- signif(floor(hi / u) * u, d)
    if (cand >= hi) cand <- signif(cand - u, d)
    if (same(cand)) return(cand)
  }
  level
}
