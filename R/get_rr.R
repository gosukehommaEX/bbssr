#' Rejection Region from Stored p-Values
#'
#' Internal helper returning the rejection region of \code{\link{BinaryRR}} as a plain
#' logical matrix. The arguments are validated as in \code{\link{BinaryRR}}, and the
#' p-values are taken from the session store kept by \code{get_pvalue}, so the region is
#' identical to the one returned by \code{\link{BinaryRR}}.
#'
#' @param N1 Sample size for group 1
#' @param N2 Sample size for group 2
#' @param alpha Level of significance
#' @param Test Type of statistical test, possibly abbreviated
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param tsmethod \code{'minlike'} or \code{'central'}
#' @param n.grid Number of grid points over the nuisance parameter
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure
#'
#' @return A logical matrix of dimension \code{(N1 + 1)} by \code{(N2 + 1)} without
#'   attributes other than \code{dim}
#'
#' @keywords internal
#' @noRd
#' @import fpCompare
get_rr <- function(N1, N2, alpha, Test, alternative, tsmethod, n.grid, bb.gamma) {
  a <- check_rr_args(N1, N2, alpha, Test, n.grid, bb.gamma)
  pv <- get_pvalue(a$N1, a$N2, a$Test, alternative, tsmethod, a$n.grid, bb.gamma)
  pv %<<% alpha
}
