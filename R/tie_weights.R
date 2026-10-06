#' Rejection Probability of a Design in Terms of the Final Responder Counts
#'
#' Internal helper expressing the rejection probability of a design through the final
#' responder counts. Given that \code{X1} of the \code{N1} patients of group 1 respond, the
#' number of responders among the first \code{n11} is hypergeometric whatever the response
#' probability, and likewise for group 2. The rejection probability is therefore
#' \code{sum_K sum_{X1, X2} C_K[X1, X2] b(X1; N1_K, p1) b(X2; N2_K, p2)}, where \code{K}
#' runs over the distinct final sample sizes, \code{b} is the binomial probability and
#' \code{C_K[X1, X2]} is the indicator of rejection times the sum, over the interim cells
#' that lead to \code{K}, of the hypergeometric probabilities of their two interim counts.
#'
#' @param setup List returned by \code{bssr_setup} or \code{fixed_setup}
#' @param rr.list Rejection regions of the final sample sizes of \code{setup}, as logical
#'   matrices
#'
#' @return A list with one element per final sample size, itself a list with the sizes
#'   \code{N1} and \code{N2} and the matrix \code{C} with \code{N1 + 1} rows and
#'   \code{N2 + 1} columns
#'
#' @keywords internal
#' @noRd
#' @importFrom stats dhyper
tie_weights <- function(setup, rr.list) {
  lapply(seq_along(rr.list), function(k) {
    N1 <- setup$N1[k]
    N2 <- setup$N2[k]
    cells <- which(setup$rr.id == k - 1L)
    # Probability of x11 responders among the first n11 patients given X1 responders among
    # the N1 patients, with one row per value of X1 and one column per interim cell
    a <- outer(0:N1, setup$x11[cells],
               function(X, x) dhyper(x, setup$n11, setup$n21[k], X))
    b <- outer(0:N2, setup$x12[cells],
               function(X, x) dhyper(x, setup$n12, setup$n22[k], X))
    list(N1 = N1, N2 = N2, C = rr.list[[k]] * tcrossprod(a, b))
  })
}
