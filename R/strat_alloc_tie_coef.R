#' Type I Error Rate of a Two-Stratum Design with Random Stratum and Group Sizes
#'
#' Internal helper returning the coefficients of the exact type I error rate of a fixed
#' two-stratum design whose stratum and group sizes are random, in the products of the
#' Bernstein bases of degree \code{N} (or the fixed stratum sizes of \code{alloc =
#' 'group'}) in the common response probabilities of the two strata. It is the sum over
#' the configurations of \code{strat_alloc_power} of their probability times the
#' coefficients of the fixed design they give.
#'
#' @inheritParams strat_alloc_power
#'
#' @return A list with the matrix \code{coef} of coefficients and the total probability
#'   \code{skipped} of the configurations left out
#'
#' @keywords internal
#' @noRd
strat_alloc_tie_coef <- function(N, al, alloc, alpha, Test, alternative, tsmethod,
                                 eps = 1e-14) {
  a <- check_strat_args(c(1, 1), c(1, 1), alpha, Test, alternative, tsmethod)
  n1 <- floor(N * al$frac[1] + 1e-9)
  strat_random_tie_coef(match(alloc, c('stratum', 'group', 'both')), as.integer(N),
                        al$frac[1], as.integer(c(n1, N - n1)), al$b,
                        if (is.null(al$block)) 4L else al$block,
                        if (is.null(al$g)) c(1L, 1L) else al$g, a$type, a$alt, a$ts,
                        alpha, eps)
}
