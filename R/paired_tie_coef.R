#' Bernstein Coefficients of the Rejection Probability on a Null Boundary
#'
#' Internal helper returning the Bernstein coefficients, of degree \code{max(Ns)} in the
#' parameter \code{t} of \code{paired_boundary}, of the rejection probability of a design
#' for paired binary data on a boundary of the null hypothesis. Each final number of pairs
#' contributes the coefficients of \code{paired_bern_coef}, of degree \code{N_k}, which
#' are raised to the common degree by the degree elevation of Bernstein polynomials
#' before they are added.
#'
#' @param rr.list Rejection regions of the distinct final numbers of pairs
#' @param Ns Distinct final numbers of pairs
#' @param W.list Weights of \code{paired_tie_weights}, or a list with the vector of ones
#'   for a fixed design
#' @param type,m Type and value of the boundary, see \code{paired_boundary}
#'
#' @return A numeric vector of length \code{max(Ns) + 1}
#'
#' @keywords internal
#' @noRd
paired_tie_coef <- function(rr.list, Ns, W.list, type, m) {
  coef <- 0
  for (k in order(Ns)) {
    h <- paired_bern_coef(rr.list[[k]], as.double(W.list[[k]]), type, m)
    n <- length(coef) - 1
    while (n < Ns[k]) {
      i <- 0:(n + 1)
      coef <- c(coef, 0) * (1 - i / (n + 1)) + c(0, coef) * (i / (n + 1))
      n <- n + 1
    }
    coef <- coef + h
  }
  coef
}
