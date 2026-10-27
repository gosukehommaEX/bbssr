#' Largest Size of Group 2 within an Upper Bound on the Total
#'
#' Internal helper returning the largest size \code{N2} of group 2 whose total
#' \code{N2 + ceiling(r N2)} does not exceed \code{N.max}, the upper bound of the rounding
#' rule \code{'group'} of \code{final_sizes}. The exact re-estimation does not search
#' beyond it.
#'
#' @param r Allocation ratio to group 1
#' @param N.max Upper bound on the final total, or \code{NULL}
#'
#' @return The largest size of group 2, \code{Inf} when \code{N.max} is \code{NULL}
#'
#' @keywords internal
#' @noRd
group2_max <- function(r, N.max) {
  if (is.null(N.max)) return(Inf)
  # An infinite bound would make the steps below run for ever
  if (length(N.max) != 1 || !is.numeric(N.max) || !is.finite(N.max)) {
    stop('N.max must be a single finite number')
  }
  hi2 <- floor(N.max / (1 + r))
  while ((hi2 + 1) + ceil_tol(r * (hi2 + 1)) <= N.max) hi2 <- hi2 + 1
  while (hi2 > 0 && hi2 + ceil_tol(r * hi2) > N.max) hi2 <- hi2 - 1
  hi2
}
