#' Unrounded Total Sample Size of a Stratified Design from the Normal Approximation
#'
#' Internal helper returning the total sample size at which the normal approximation to
#' the power of Cochran's test reaches the target power, formula (2.8) of Woolson et al.
#' (1986), or with the continuity correction, formula (2.3) of Nam (1992). With
#' \code{s[j] = r[j] / (1 + r[j])}, the fraction of group 1 in stratum \code{j}, and
#' \code{w[j] = frac[j] s[j] (1 - s[j])},
#' \code{D = sum_j w[j] (p1[j] - p2[j])},
#' \code{v0 = sum_j w[j] pbar[j] (1 - pbar[j])} with
#' \code{pbar[j] = s[j] p1[j] + (1 - s[j]) p2[j]}, and
#' \code{v1 = sum_j w[j] ((1 - s[j]) p1[j] (1 - p1[j]) + s[j] p2[j] (1 - p2[j]))}, the
#' uncorrected size is \code{N* = (z.a sqrt(v0) + z.b sqrt(v1))^2 / D^2}, where \code{z.a}
#' is the upper \code{alpha} quantile (\code{alpha / 2} for a two-sided alternative) and
#' \code{z.b} the \code{tar.power} quantile of the standard normal distribution. The
#' corrected size is \code{N = (1 + sqrt(1 + 2 / (D N*)))^2 N* / 4}. For the alternative
#' \code{'less'} the difference \code{D} changes its sign, and for \code{'two.sided'} its
#' absolute value is used. The products \code{p (1 - p)} are set to 0 when negative, so
#' that probabilities outside the unit interval, which keep an assumed risk difference,
#' can be used as in \code{ss_raw_parts}.
#'
#' @param p1 Response probabilities of group 1, a matrix with one row per scenario and one
#'   column per stratum, or a vector for a single scenario
#' @param p2 Response probabilities of group 2, of the same form
#' @param frac Fractions of the strata, summing to 1
#' @param r Allocation ratios of the strata
#' @param alpha Level of significance
#' @param tar.power Target power
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param cc Logical. Whether the size of the continuity-corrected test is returned
#'
#' @return A numeric vector with one element per scenario, \code{NaN} when the difference
#'   \code{D} is not positive in the direction of the alternative
#'
#' @references
#' Woolson RF, Bean JA, Rojas PB (1986). Sample size for case-control studies using
#' Cochran's statistic. \emph{Biometrics}, 42(4), 927-932.
#'
#' Nam J (1992). Sample size determination for case-control studies and the comparison of
#' stratified and unstratified analyses. \emph{Biometrics}, 48(2), 389-395.
#'
#' @keywords internal
#' @noRd
#' @importFrom stats qnorm
strat_raw_n <- function(p1, p2, frac, r, alpha, tar.power, alternative, cc) {
  J <- length(frac)
  p1 <- matrix(p1, ncol = J)
  p2 <- matrix(p2, ncol = J)
  s <- r / (1 + r)
  w <- frac * s * (1 - s)
  D <- drop((p1 - p2) %*% w)
  if (alternative == 'less') D <- -D
  if (alternative == 'two.sided') D <- abs(D)
  pbar <- sweep(p1, 2, s, '*') + sweep(p2, 2, 1 - s, '*')
  v0 <- drop(pmax(pbar * (1 - pbar), 0) %*% w)
  v1 <- drop((sweep(pmax(p1 * (1 - p1), 0), 2, 1 - s, '*') +
                sweep(pmax(p2 * (1 - p2), 0), 2, s, '*')) %*% w)
  z.a <- qnorm(1 - if (alternative == 'two.sided') alpha / 2 else alpha)
  z.b <- qnorm(tar.power)
  n <- (z.a * sqrt(v0) + z.b * sqrt(v1))^2 / D^2
  if (cc) n <- (1 + sqrt(1 + 2 / (D * n)))^2 * n / 4
  n[!(D > 0)] <- NaN
  n
}
