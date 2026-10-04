#' Group-Specific Response Probabilities from a Pooled Probability
#'
#' Internal helper returning the response probabilities of the two groups that have a
#' given pooled probability \code{p = (r p1 + p2) / (1 + r)} and a given treatment effect.
#' The function is vectorized over \code{p}.
#'
#' For the risk difference \code{Delta = p1 - p2},
#' \code{p1 = p + Delta / (1 + r)} and \code{p2 = p - r Delta / (1 + r)}. For the risk
#' ratio \code{Delta = p1 / p2}, \code{p2 = (1 + r) p / (1 + r Delta)} and
#' \code{p1 = Delta p2}, formula (21.11) of Kieser (2020). For the odds ratio
#' \code{Delta = p1 (1 - p2) / (p2 (1 - p1))}, \code{p2} is the root in the unit interval
#' of \code{(Delta - 1) x^2 + b x - (1 + r) p = 0} with
#' \code{b = 1 + r Delta - (1 + r) p (Delta - 1)}, evaluated as
#' \code{2 (1 + r) p / (b + sqrt(b^2 + 4 (Delta - 1) (1 + r) p))} to avoid cancellation,
#' and \code{p1 = Delta p2 / (1 - p2 + Delta p2)}.
#'
#' The risk difference and the risk ratio can give probabilities outside the unit
#' interval, which are returned unchanged so that the caller can decide how to treat them.
#'
#' @param p Pooled response probability
#' @param Delta Treatment effect on the scale given by \code{effect}
#' @param r Allocation ratio to group 1
#' @param effect \code{'RD'}, \code{'RR'} or \code{'OR'}
#'
#' @return A list with components \code{p1} and \code{p2}
#'
#' @keywords internal
#' @noRd
split_pooled <- function(p, Delta, r, effect) {
  if (effect == 'RD') {
    # The expressions match those of version 2.0.0 term by term, so that the recovered
    # probabilities, and with them the re-estimated sample sizes, are unchanged
    p1 <- p + (1 / (1 + r)) * Delta
    p2 <- p - (r / (1 + r)) * Delta
  } else if (effect == 'RR') {
    p2 <- (1 + r) * p / (1 + r * Delta)
    p1 <- Delta * p2
  } else {
    b <- 1 + r * Delta - (1 + r) * p * (Delta - 1)
    disc <- pmax(b^2 + 4 * (Delta - 1) * (1 + r) * p, 0)
    p2 <- 2 * (1 + r) * p / (b + sqrt(disc))
    p1 <- Delta * p2 / (1 - p2 + Delta * p2)
  }
  list(p1 = p1, p2 = p2)
}
