#' Clopper-Pearson Bounds for the Pooled Response Probability
#'
#' Internal helper used by the Berger-Boos procedure. For every possible value of the
#' total number of responders in a trial of \code{N} patients, the function returns the
#' exact \code{100 (1 - gamma)} percent confidence bounds for the common response
#' probability under the null hypothesis.
#'
#' @param N Total sample size of both groups
#' @param gamma Confidence level parameter of the Berger-Boos procedure
#'
#' @return A numeric matrix with \code{N + 1} rows and two columns named \code{lower} and
#'   \code{upper}, the row \code{s + 1} corresponding to \code{s} responders
#'
#' @keywords internal
#' @noRd
#' @importFrom stats qbeta
cp_bounds <- function(N, gamma) {
  s <- 0:N
  lower <- qbeta(gamma / 2, pmax(s, 1), N - s + 1)
  upper <- qbeta(1 - gamma / 2, s + 1, pmax(N - s, 1))
  lower[s == 0] <- 0
  upper[s == N] <- 1
  cbind(lower = lower, upper = upper)
}
