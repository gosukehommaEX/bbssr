#' Certified Largest Type I Error Rate for Paired Binary Data
#'
#' Internal helper returning the largest rejection probability of a design for paired
#' binary data over an interval of the probability \code{psi} of a discordant pair on a
#' boundary of the null hypothesis, together with an upper bound. The rejection
#' probability is a polynomial in the parameter \code{t} of \code{paired_boundary}. Its
#' Bernstein coefficients over \code{[0, 1]} are restricted to the interval by
#' \code{bern_restrict} and passed to \code{bernstein_max}, as \code{tie_certify} does for
#' two independent groups.
#'
#' @param rr.list Rejection regions of the distinct final numbers of pairs
#' @param Ns Distinct final numbers of pairs
#' @param W.list Weights of \code{paired_tie_weights}
#' @param interval Interval of \code{psi}, within the boundary
#' @param alternative,margin Boundary, see \code{paired_boundary}
#' @param tol Tolerance of \code{bernstein_max}
#'
#' @return A list with the location \code{x} (a value of \code{psi}), the largest value
#'   \code{y} and the upper bound \code{bound}
#'
#' @keywords internal
#' @noRd
paired_tie_certify <- function(rr.list, Ns, W.list, interval, alternative, margin,
                               tol = 1e-12) {
  b <- paired_boundary(interval, alternative, margin)
  coef <- paired_tie_coef(rr.list, Ns, W.list, b$type, b$m)
  coef <- bern_restrict(coef, b$t[1], b$t[2])
  mx <- bernstein_max(coef, tol, 100000L)
  if (!mx$complete) {
    warning('the certified maximization stopped at its limit of subdivisions; ',
            'the upper bound is valid but may exceed the maximum by more than ', tol)
  }
  t <- b$t[1] + (b$t[2] - b$t[1]) * mx$x
  list(x = b$m + (1 - b$m) * t, y = mx$value, bound = mx$bound)
}
