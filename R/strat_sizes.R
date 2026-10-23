#' Group Sizes of the Strata for a Nominal Total Sample Size
#'
#' Internal helper returning the group sizes of the strata of a stratified design with
#' the nominal total sample size \code{n}. The size of group 2 in stratum \code{j} is
#' \code{ceiling(frac[j] n / (1 + r[j]))} and that of group 1 is
#' \code{ceiling(r[j] N2[j])}, both at least 1, with the tolerance of \code{ceil_tol}. The
#' total of the group sizes is at least \code{n}, and the sizes do not decrease with
#' \code{n}, so the designs for \code{n = 1, 2, ...} form the family within which the
#' sample size is searched and re-estimated.
#'
#' @param n Nominal total sample sizes, a vector
#' @param frac Fractions of the strata, summing to 1
#' @param r Allocation ratios of the strata
#'
#' @return A list with the integer matrices \code{N1} and \code{N2}, with one row for each
#'   element of \code{n} and one column per stratum
#'
#' @keywords internal
#' @noRd
strat_sizes <- function(n, frac, r) {
  # The matrix comes first in pmax(), which keeps the attributes of its first argument
  N2 <- pmax(ceil_tol(outer(n, frac / (1 + r))), 1)
  N1 <- pmax(ceil_tol(sweep(N2, 2, r, '*')), 1)
  storage.mode(N1) <- 'integer'
  storage.mode(N2) <- 'integer'
  list(N1 = N1, N2 = N2)
}
