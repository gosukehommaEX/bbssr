#' Validate the Arguments of a Rejection Region for Paired Data
#'
#' Internal helper shared by \code{\link{PairedRR}} and \code{get_paired_rr}. It checks the
#' number of pairs, the level, the grid size, the Berger-Boos parameter, the refinement
#' flag and the margin, warns when the Berger-Boos parameter is supplied for a test other
#' than \code{'McNemar-uncond'}, and returns the arguments in the form used by the
#' computation. A margin other than 0 requires \code{Test = 'Paired-Score'} and a
#' one-sided alternative, and a test of equivalence requires \code{Test = 'Paired-Score'};
#' its margin is validated and returned as the pair \code{c(lower, upper)} by
#' \code{norm_margin}. The refinement flag is set to \code{FALSE} for the tests without a
#' nuisance parameter, so that their p-values are stored only once.
#'
#' @param N Number of pairs
#' @param alpha Level of significance
#' @param Test Type of statistical test, possibly abbreviated
#' @param n.grid Number of grid points over the nuisance parameter
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure
#' @param ref.pvalue Logical. Whether the maximum over the nuisance parameter is refined
#'   between the grid points
#' @param alternative \code{'greater'}, \code{'less'}, \code{'two.sided'} or
#'   \code{'equivalence'}
#' @param margin Margin on the scale of the difference \code{p10 - p01}, 0 for a test of
#'   superiority, or the margins of a test of equivalence
#'
#' @return A list with the integer number of pairs \code{N}, the integer grid size
#'   \code{n.grid}, the full name of the test \code{Test}, the refinement flag
#'   \code{ref.pvalue} and the margin \code{margin}
#'
#' @keywords internal
#' @noRd
check_paired_args <- function(N, alpha, Test, n.grid, bb.gamma, ref.pvalue, alternative,
                              margin) {
  Test <- match.arg(Test, c('McNemar', 'McNemar-cc', 'McNemar-exact', 'McNemar-midP',
                            'McNemar-uncond', 'Paired-Score'))
  if (length(N) != 1 || is.na(N) || N != round(N) || N < 1) {
    stop('N must be a single positive integer')
  }
  if (length(alpha) != 1 || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop('alpha must be a single value in (0, 1)')
  }
  if (length(n.grid) != 1 || is.na(n.grid) || n.grid < 2) {
    stop('n.grid must be a single integer of at least 2')
  }
  if (length(bb.gamma) != 1 || is.na(bb.gamma) || bb.gamma < 0 || bb.gamma >= alpha) {
    stop('bb.gamma must be a single value satisfying 0 <= bb.gamma < alpha')
  }
  if (length(ref.pvalue) != 1 || !is.logical(ref.pvalue) || is.na(ref.pvalue)) {
    stop('ref.pvalue must be TRUE or FALSE')
  }
  if (alternative == 'equivalence') {
    if (Test != 'Paired-Score') {
      stop("alternative = 'equivalence' requires Test = 'Paired-Score'")
    }
    margin <- norm_margin(alternative, margin, 'RD')
  } else {
    if (length(margin) != 1 || !is.numeric(margin) || is.na(margin) || abs(margin) >= 1) {
      stop('margin must be a single value in (-1, 1)')
    }
    if (margin != 0 && Test != 'Paired-Score') {
      stop("a non-zero margin requires Test = 'Paired-Score'")
    }
    if (margin != 0 && alternative == 'two.sided') {
      stop('a non-zero margin requires a one-sided alternative')
    }
  }
  if (bb.gamma > 0 && Test != 'McNemar-uncond') {
    warning('bb.gamma applies only to the McNemar-uncond test and is ignored here')
  }
  list(N = as.integer(N), n.grid = as.integer(n.grid), Test = Test,
       ref.pvalue = ref.pvalue && Test == 'McNemar-uncond', margin = as.double(margin))
}
