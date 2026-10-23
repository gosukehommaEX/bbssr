#' Certified Maximum of the Type I Error Rate of a Two-Stratum Design
#'
#' Internal helper returning the largest type I error rate of a two-stratum design over a
#' rectangle of common response probabilities of the strata, together with an upper bound
#' that the largest rate cannot exceed. The rate is the bivariate polynomial whose
#' Bernstein coefficients \code{coef} are returned by \code{strat_tie_coef}; its maximum is
#' located by \code{strat_bern_max}, which subdivides the rectangle until the bound exceeds
#' the largest value found by at most \code{tol}.
#'
#' @param coef Matrix of Bernstein coefficients over the unit square
#' @param box Matrix with two rows, the smallest and largest probability, and one column
#'   per stratum, with a positive width in each column
#' @param tol Largest difference between the bound and the value. Default is \code{1e-12}
#'
#' @return A list with the location \code{x} (two probabilities), the value \code{y} and
#'   the upper bound \code{bound}
#'
#' @keywords internal
#' @noRd
strat_tie_certify <- function(coef, box, tol = 1e-12) {
  m <- strat_bern_max(coef, box[1, 1], box[2, 1], box[1, 2], box[2, 2], tol, 100000L)
  if (!m$complete) {
    warning('the certified maximization stopped at its limit of subdivisions; ',
            'the upper bound is valid but may exceed the maximum by more than ', tol)
  }
  list(x = c(m$x, m$y), y = m$value, bound = m$bound)
}
