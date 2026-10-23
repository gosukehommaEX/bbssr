#' p-Values of a Test for Paired Binary Data over the Whole Outcome Grid
#'
#' Internal helper computing the p-value of every outcome of \code{N} pairs for one of the
#' tests of \code{\link{PairedRR}}. The alternative \code{'less'} is obtained by
#' exchanging the two discordant counts, which turns the null hypothesis
#' \code{p10 - p01 >= margin} into \code{p01 - p10 <= -margin}, and a test of equivalence
#' by the larger of the p-values of the two one-sided tests.
#'
#' @param N Number of pairs
#' @param Test Full name of the test
#' @param alternative \code{'greater'}, \code{'less'}, \code{'two.sided'} or
#'   \code{'equivalence'}
#' @param n.grid Number of grid points over the probability of a discordant pair
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure
#' @param ref.pvalue Logical. Whether the maximum over the probability of a discordant pair
#'   is refined between the grid points
#' @param margin Margin, or the pair \code{c(lower, upper)} of margins of a test of
#'   equivalence
#'
#' @return A numeric matrix of dimension \code{(N + 1)} by \code{(N + 1)} whose element in
#'   row \code{a + 1} and column \code{b + 1} is the p-value of the outcome
#'   \code{(a, b)}, and \code{NA} where \code{a + b > N}
#'
#' @keywords internal
#' @noRd
#' @importFrom stats pnorm
paired_pvalue <- function(N, Test, alternative, n.grid, bb.gamma, ref.pvalue, margin) {
  if (alternative == 'equivalence') {
    # Two one-sided tests, each at the level alpha: the null hypothesis is rejected when
    # both reject, so the p-value is the larger of the two one-sided p-values
    return(pmax(paired_pvalue(N, Test, 'greater', n.grid, bb.gamma, ref.pvalue, margin[1]),
                paired_pvalue(N, Test, 'less', n.grid, bb.gamma, ref.pvalue, margin[2])))
  }
  if (alternative == 'less') {
    return(t(paired_pvalue(N, Test, 'greater', n.grid, bb.gamma, ref.pvalue, margin)))
  }
  if (Test %in% c('McNemar', 'Paired-Score')) {
    Z <- paired_score_z(N, if (Test == 'McNemar') 0 else margin)
    if (alternative == 'greater') {
      pnorm(Z, lower.tail = FALSE)
    } else {
      pmin(2 * pnorm(abs(Z), lower.tail = FALSE), 1)
    }
  } else if (Test %in% c('McNemar-exact', 'McNemar-midP')) {
    paired_cond_pvalue(N, alternative, midp = Test == 'McNemar-midP')
  } else {
    Z <- paired_score_z(N, 0)
    stat <- if (alternative == 'greater') Z else abs(Z)
    paired_uncond_pvalue(stat, N, n.grid, bb.gamma, ref.pvalue)
  }
}
