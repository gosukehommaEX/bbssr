#' Z Statistic for a Null Hypothesis of a Non-Zero Risk Difference
#'
#' Internal helper returning, at every cell of the \code{(N1 + 1)} by \code{(N2 + 1)}
#' outcome grid, the statistic
#' \code{z = (hat.p1 - hat.p2 + margin) / SE} for the null hypothesis
#' \code{p1 - p2 <= -margin}. Under \code{se = 'unpooled'} the standard error uses the
#' observed proportions, as in Blackwelder (1982). Under \code{se = 'restricted'} it uses
#' the maximum likelihood estimates under the restriction \code{p1 - p2 = -margin}, as in
#' Farrington and Manning (1990). With \code{margin = 0} the restricted estimates are the
#' pooled proportion, so the statistic of \code{zstat} is returned.
#'
#' A standard error of 0 arises only when both proportions are 0 or 1. The statistic is
#' then set to \code{Inf} or \code{-Inf} according to the sign of the numerator, and to 0
#' when the numerator is 0 as well.
#'
#' @param N1 Sample size for group 1
#' @param N2 Sample size for group 2
#' @param margin Non-inferiority margin, a single value in \code{(-1, 1)}
#' @param se \code{'unpooled'} or \code{'restricted'}
#'
#' @return A numeric matrix of dimension \code{(N1 + 1)} by \code{(N2 + 1)}
#'
#' @references
#' Blackwelder WC (1982). "Proving the null hypothesis" in clinical trials.
#' \emph{Controlled Clinical Trials}, 3(4), 345-353.
#'
#' Farrington CP, Manning G (1990). Test statistics and sample size formulae for
#' comparative binomial trials with null hypothesis of non-zero risk difference or
#' non-unity relative risk. \emph{Statistics in Medicine}, 9(12), 1447-1454.
#'
#' @keywords internal
#' @noRd
zstat_margin <- function(N1, N2, margin, se) {
  if (se == 'restricted' && margin == 0) return(zstat(N1, N2))
  hat.p1 <- (0:N1) / N1
  hat.p2 <- (0:N2) / N2
  num <- outer(hat.p1, hat.p2, '-') + margin
  if (se == 'unpooled') {
    v <- outer(hat.p1 * (1 - hat.p1) / N1, hat.p2 * (1 - hat.p2) / N2, '+')
  } else {
    P1 <- rep(hat.p1, times = N2 + 1)
    P2 <- rep(hat.p2, each = N1 + 1)
    est <- fm_restricted(P1, P2, N2 / N1, -margin)
    v <- pmax(est$p1 * (1 - est$p1), 0) / N1 + pmax(est$p2 * (1 - est$p2), 0) / N2
    dim(v) <- c(N1 + 1L, N2 + 1L)
  }
  Z <- num / sqrt(v)
  zero <- v <= 0
  Z[zero] <- ifelse(num[zero] == 0, 0, sign(num[zero]) * Inf)
  Z
}
