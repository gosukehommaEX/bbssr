#' Rejection Region of a Test for Paired Binary Data
#'
#' Internal helper validating the arguments with \code{check_paired_args} and returning
#' the rejection region at the level \code{alpha}, the outcomes whose p-value is smaller
#' than \code{alpha} in the comparison \code{\%<<\%} of the \pkg{fpCompare} package, as in
#' \code{get_rr}.
#'
#' @inheritParams PairedRR
#'
#' @return A logical matrix of dimension \code{(N + 1)} by \code{(N + 1)}, \code{TRUE} at
#'   the rejected outcomes \code{(a, b)} and \code{FALSE} elsewhere, including the
#'   impossible outcomes \code{a + b > N}
#'
#' @keywords internal
#' @noRd
get_paired_rr <- function(N, alpha, Test, alternative, n.grid, bb.gamma, ref.pvalue,
                          margin) {
  a <- check_paired_args(N, alpha, Test, n.grid, bb.gamma, ref.pvalue, alternative, margin)
  pv <- get_paired_pvalue(a$N, a$Test, alternative, a$n.grid, bb.gamma, a$ref.pvalue,
                          a$margin)
  rr <- pv %<<% alpha
  rr[is.na(rr)] <- FALSE
  rr
}
