#' Boundary of the Null Hypothesis for Paired Binary Data
#'
#' Internal helper returning the probabilities of the two kinds of discordant pairs on the
#' boundary of the null hypothesis, as functions of the probability \code{psi} of a
#' discordant pair. The difference on the boundary is \code{d0 = -margin} for
#' \code{alternative = 'greater'}, \code{d0 = margin} for \code{'less'} and \code{d0 = 0}
#' for \code{'two.sided'}, and \code{p10 = (psi + d0) / 2} and \code{p01 = (psi - d0) / 2}.
#' A value of \code{psi} lies on the boundary when \code{|d0| <= psi <= 1}. The boundary is
#' also parametrized by \code{t = (psi - |d0|) / (1 - |d0|)} in \code{[0, 1]}, in which the
#' probabilities \code{p10}, \code{p01} and \code{1 - psi} are linear, as required by the
#' Bernstein coefficients of \code{paired_bern_coef}: \code{type} is 0 when \code{d0 = 0},
#' 1 when \code{d0 < 0} and 2 when \code{d0 > 0}, and \code{m = |d0|}.
#'
#' @param psi Values of the probability of a discordant pair
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param margin Margin of the one-sided test
#'
#' @return A list with the probabilities \code{p10} and \code{p01}, the logical vector
#'   \code{ok} of the values of \code{psi} on the boundary, the parameter \code{t}, the
#'   type \code{type}, the value \code{m} and the interval \code{range} of \code{psi}
#'   between the smallest and the largest value on the boundary within the range of
#'   \code{psi}, extended to the values of \code{psi} accepted as on the boundary
#'
#' @keywords internal
#' @noRd
paired_boundary <- function(psi, alternative, margin) {
  d0 <- switch(alternative, greater = -margin, less = margin, 0)
  m <- abs(d0)
  type <- if (d0 == 0) 0L else if (d0 < 0) 1L else 2L
  ok <- (psi %>=% m) & (psi %<=% 1)
  p10 <- pmin(1, pmax(0, (psi + d0) / 2))
  p01 <- pmin(1, pmax(0, (psi - d0) / 2))
  t <- pmin(1, pmax(0, (psi - m) / (1 - m)))
  range <- if (any(ok)) {
    c(min(max(min(psi), m), psi[ok]), max(min(max(psi), 1), psi[ok]))
  } else {
    c(NA_real_, NA_real_)
  }
  list(p10 = p10, p01 = p01, ok = ok, t = t, type = type, m = m, range = range)
}
