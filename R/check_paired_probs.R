#' Validate the Probabilities of Discordant Pairs
#'
#' Internal helper checking that \code{p10} and \code{p01} are vectors of the same length
#' with elements in \code{[0, 1]} and sums at most 1. A sum that exceeds 1 only by
#' rounding error is accepted.
#'
#' @param p10 Probability that only treatment 1 responds
#' @param p01 Probability that only treatment 2 responds
#'
#' @return \code{NULL}, invisibly. An error is raised when a check fails
#'
#' @keywords internal
#' @noRd
check_paired_probs <- function(p10, p01) {
  if (!is.numeric(p10) || !is.numeric(p01) || length(p10) != length(p01) ||
      length(p10) == 0) {
    stop('p10 and p01 must be numeric vectors of the same length')
  }
  if (anyNA(p10) || anyNA(p01) || any(p10 < 0 | p10 > 1 | p01 < 0 | p01 > 1)) {
    stop('p10 and p01 must lie in [0, 1]')
  }
  if (any((p10 + p01) %>>% 1)) stop('p10 + p01 must not exceed 1')
  invisible(NULL)
}
