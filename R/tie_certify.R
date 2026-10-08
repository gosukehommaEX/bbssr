#' Certified Maximum of the Type I Error Rate over an Interval
#'
#' Internal helper returning the largest rejection probability of a design over an
#' interval of pooled response probabilities on the boundary of the null hypothesis,
#' together with an upper bound that the largest value cannot exceed. On the boundary both
#' response probabilities are linear in the pooled probability, so the rejection
#' probability is a polynomial whose coefficients in the Bernstein basis over the interval
#' follow from \code{tie_bernstein}. The maximum is located by \code{bernstein_max}, which
#' stops once the bound exceeds the largest value found by at most \code{tol}.
#'
#' @param weights List returned by \code{tie_weights}
#' @param interval The two ends of the interval of pooled response probabilities, as
#'   returned by \code{null_range}
#' @param r Allocation ratio to group 1
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param margin Non-inferiority margin, a difference for \code{margin.scale = 'RD'} and
#'   a ratio for \code{'RR'}
#' @param tol Largest difference between the bound and the value. Default is \code{1e-12}
#' @param margin.scale \code{'RD'} or \code{'RR'}. On both scales the response
#'   probabilities on the boundary are linear in the pooled probability
#'
#' @return A list with the pooled response probability \code{x} of the largest value
#'   found, the value \code{y} and the upper bound \code{bound}
#'
#' @keywords internal
#' @noRd
tie_certify <- function(weights, interval, r, alternative, margin, tol = 1e-12,
                        margin.scale) {
  b <- null_boundary(interval, r, alternative, margin, margin.scale)
  m <- bernstein_max(tie_bernstein(weights, b$p1, b$p2), tol, 100000L)
  if (!m$complete) {
    warning('the certified maximization stopped at its limit of subdivisions; ',
            'the upper bound is valid but may exceed the maximum by more than ', tol)
  }
  list(x = interval[1] + (interval[2] - interval[1]) * m$x, y = m$value, bound = m$bound)
}
