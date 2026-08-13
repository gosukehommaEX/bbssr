#' Pooled Z Statistic over the Whole Outcome Grid
#'
#' Internal helper returning the Z statistic of the two-sample test for proportions with a
#' pooled variance estimator, evaluated at every cell of the \code{(N1 + 1)} by
#' \code{(N2 + 1)} outcome grid. Cells for which the pooled proportion is 0 or 1 give a
#' zero denominator, and the statistic is set to 0 there.
#'
#' @param N1 Sample size for group 1
#' @param N2 Sample size for group 2
#'
#' @return A numeric matrix of dimension \code{(N1 + 1)} by \code{(N2 + 1)}
#'
#' @keywords internal
#' @noRd
zstat <- function(N1, N2) {
  x1 <- 0:N1
  x2 <- 0:N2
  diff.hat.p <- outer(x1 / N1, x2 / N2, '-')
  hat.p <- outer(x1, x2, '+') / (N1 + N2)
  se <- sqrt(hat.p * (1 - hat.p) * (1 / N1 + 1 / N2))
  Z <- diff.hat.p / se
  Z[!is.finite(Z)] <- 0
  Z
}
