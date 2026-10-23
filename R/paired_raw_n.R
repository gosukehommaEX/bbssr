#' Unrounded Number of Pairs from the Normal Approximation
#'
#' Internal helper returning the unrounded number of pairs at which the normal
#' approximation to the power of the score test of \code{\link{PairedRR}} reaches the
#' target power. For a one-sided alternative it is formula (9) of Nam (1997),
#' \code{(z.a sqrt(v0) + z.b sqrt(v1))^2 / d^2} with the ingredients of
#' \code{paired_raw_parts}, \code{z.a} the upper \code{alpha} quantile and \code{z.b} the
#' \code{tar.power} quantile of the standard normal distribution. With \code{margin = 0}
#' it is formula (3) of Connor (1987). For \code{alternative = 'two.sided'}, where the
#' margin is 0, \code{z.a} is the upper \code{alpha / 2} quantile and \code{d} is the
#' absolute difference. For \code{alternative = 'equivalence'} each one-sided test is at
#' the level \code{alpha} and the approximate power at \code{n} pairs is
#' \code{Phi((d_L sqrt(n) - z.a sqrt(v0_L)) / sqrt(v1)) +
#' Phi((d_U sqrt(n) - z.a sqrt(v0_U)) / sqrt(v1)) - 1}, a lower bound for the
#' approximate probability that both tests reject, where \code{L} refers to the test
#' against \code{'greater'} with the margin \code{lower} and \code{U} to the test against
#' \code{'less'} with the margin \code{upper}. The size is the root of the approximation
#' minus the target power, found as in \code{ss_raw_n2_equiv}. A term whose variance
#' \code{v1} is 0 is 1 when its numerator is positive and 0 otherwise. The function is
#' vectorized over \code{p10} and \code{p01}.
#'
#' @param p10 Probability that only treatment 1 responds
#' @param p01 Probability that only treatment 2 responds
#' @param alpha Level of significance
#' @param tar.power Target power
#' @param alternative \code{'greater'}, \code{'less'}, \code{'two.sided'} or
#'   \code{'equivalence'}
#' @param margin Margin, or the pair \code{c(lower, upper)} of margins of a test of
#'   equivalence
#'
#' @return A numeric vector. Probabilities on the wrong side of the null boundary give
#'   \code{NaN} for a test of equivalence and a value that is not meaningful otherwise, so
#'   the callers check them first
#'
#' @references
#' Connor RJ (1987). Sample size for testing differences in proportions for the
#' paired-sample design. \emph{Biometrics}, 43(1), 207-211.
#'
#' Nam J (1997). Establishing equivalence of two treatments and sample size requirements in
#' matched-pairs design. \emph{Biometrics}, 53(4), 1422-1430.
#'
#' @keywords internal
#' @noRd
#' @importFrom stats pnorm qnorm uniroot
paired_raw_n <- function(p10, p01, alpha, tar.power, alternative, margin) {
  if (alternative != 'equivalence') {
    alpha.eff <- if (alternative == 'two.sided') alpha / 2 else alpha
    s <- paired_raw_parts(p10, p01, if (alternative == 'less') 'less' else 'greater',
                          margin)
    d <- if (alternative == 'two.sided') abs(s$d) else s$d
    return((qnorm(1 - alpha.eff) * sqrt(s$v0) + qnorm(tar.power) * sqrt(s$v1)) ^ 2 / d ^ 2)
  }
  m <- norm_margin('equivalence', margin, 'RD')
  z.a <- qnorm(1 - alpha)
  n <- max(length(p10), length(p01))
  lo <- lapply(paired_raw_parts(p10, p01, 'greater', m[1]), rep_len, length.out = n)
  up <- lapply(paired_raw_parts(p10, p01, 'less', m[2]), rep_len, length.out = n)
  term <- function(d, v0, v1, x) {
    num <- d * sqrt(x) - z.a * sqrt(v0)
    if (v1 > 0) pnorm(num / sqrt(v1)) else as.numeric(num > 0)
  }
  vapply(seq_len(n), function(i) {
    v <- c(lo$d[i], lo$v0[i], lo$v1[i], up$d[i], up$v0[i], up$v1[i])
    if (anyNA(v) || any(!is.finite(v))) return(NA_real_)
    if (!(lo$d[i] > 0 && up$d[i] > 0)) return(NaN)
    gap <- function(x) {
      term(lo$d[i], lo$v0[i], lo$v1[i], x) + term(up$d[i], up$v0[i], up$v1[i], x) - 1 -
        tar.power
    }
    # With alpha above 0.5 the approximation can exceed the target without pairs
    if (gap(0) >= 0) return(0)
    hi <- 1
    while (gap(hi) < 0) hi <- 2 * hi
    uniroot(gap, c(0, hi), tol = 1e-12)$root
  }, numeric(1))
}
