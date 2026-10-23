#' Score Statistic for Paired Binary Data over the Whole Outcome Grid
#'
#' Internal helper returning the score statistic of Nam (1997) and Tango (1998) for the
#' null hypothesis \code{p10 - p01 = -margin} at every outcome of \code{N} pairs. The
#' outcome is given by the numbers \code{a} of pairs in which only treatment 1 responds
#' and \code{b} of pairs in which only treatment 2 responds. The statistic is
#' \code{Z = (a - b + N margin) / sqrt(N (2 q01 - margin (margin + 1)))}, where \code{q01}
#' is the larger root of
#' \code{2 N x^2 - (a + b + (2 N - a + b) margin) x + b margin (margin + 1) = 0}, the
#' maximum likelihood estimate of \code{p01} under the null hypothesis, formulae (24) to
#' (26) of Tango (1998). Tango (1998, Appendix II) shows that the larger root is the
#' estimate for either sign of the margin and for every outcome, including outcomes with
#' a zero count. With \code{margin = 0} the statistic is the McNemar statistic
#' \code{(a - b) / sqrt(a + b)}, which is then computed directly. The statistic is set to 0
#' at the outcome without discordant pairs when \code{margin = 0}, the only outcome at
#' which its denominator vanishes.
#'
#' @param N Number of pairs
#' @param margin Margin, so that the null hypothesis is \code{p10 - p01 = -margin}
#'
#' @return A numeric matrix of dimension \code{(N + 1)} by \code{(N + 1)} whose element in
#'   row \code{a + 1} and column \code{b + 1} is the statistic of the outcome
#'   \code{(a, b)}, and \code{NA} where \code{a + b > N}
#'
#' @references
#' Nam J (1997). Establishing equivalence of two treatments and sample size requirements in
#' matched-pairs design. \emph{Biometrics}, 53(4), 1422-1430.
#'
#' Tango T (1998). Equivalence test and confidence interval for the difference in
#' proportions for the paired-sample design. \emph{Statistics in Medicine}, 17(8),
#' 891-908.
#'
#' @keywords internal
#' @noRd
paired_score_z <- function(N, margin) {
  Z <- matrix(NA_real_, N + 1, N + 1)
  # Possible outcomes only, so that no square root of a negative number is taken
  ok <- outer(0:N, 0:N, '+') <= N
  a <- row(Z)[ok] - 1
  b <- col(Z)[ok] - 1
  if (margin == 0) {
    z <- (a - b) / sqrt(a + b)
    z[a + b == 0] <- 0
  } else {
    B <- -(a + b) - (2 * N - a + b) * margin
    C <- b * margin * (margin + 1)
    q01 <- (sqrt(pmax(0, B ^ 2 - 8 * N * C)) - B) / (4 * N)
    z <- (a - b + N * margin) / sqrt(N * (2 * q01 - margin * (margin + 1)))
  }
  Z[ok] <- z
  Z
}
