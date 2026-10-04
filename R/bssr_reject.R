#' Rejection Probability of a Re-Estimation Design
#'
#' Internal helper evaluating the probability that a design with blinded sample size
#' re-estimation rejects the null hypothesis, for each pair of response probabilities.
#'
#' @param setup List returned by \code{bssr_setup}
#' @param rr.list Rejection regions of the distinct final sample sizes of \code{setup},
#'   as logical matrices
#' @param p1 Response probabilities of group 1
#' @param p2 Response probabilities of group 2, of the same length as \code{p1}
#'
#' @return A numeric vector of rejection probabilities
#'
#' @keywords internal
#' @noRd
bssr_reject <- function(setup, rr.list, p1, p2) {
  bssr_power(rr.list, setup$rr.id, setup$x11, setup$x12, setup$n21, setup$n22,
             as.double(p1), as.double(p2), setup$n11, setup$n12)
}
