#' Final Sample Sizes from Re-Estimated Sample Sizes
#'
#' Internal helper turning re-estimated sample sizes into the final sample sizes of a
#' design with blinded sample size re-estimation, applying the lower and the upper bound
#' on the final size and the rounding rule.
#'
#' The final total is never below the interim total. Under the restricted rule it is
#' also never below the planned total, and \code{N.min} and \code{N.max} add a further
#' lower and upper bound.
#' \describe{
#'   \item{group}{The size of group 2 is bounded from below by the interim size of group 2,
#'     by the planned size of group 2 under the restricted rule and by
#'     \code{ceiling(N.min / (1 + r))}, and from above by the largest size whose total
#'     \code{N2 + ceiling(r N2)} does not exceed \code{N.max}. Group 1 receives
#'     \code{ceiling(r N2)} patients, and at least its interim size.}
#'   \item{friede-kieser}{The unrounded total is truncated to the bounds. The second stage
#'     receives the rounded-up excess of this total over the interim total, split into
#'     \code{ceiling(n2 / (1 + r))} patients of group 2 and \code{ceiling(r n2 / (1 + r))}
#'     patients of group 1, as in Friede and Kieser (2004). The total can therefore exceed
#'     \code{N.max} by one patient.}
#'   \item{total}{The unrounded total is truncated to the bounds and rounded up. Group 2
#'     receives \code{floor(N / (1 + r))} patients and group 1 the remainder, each at least
#'     its interim size.}
#' }
#'
#' @param re Data frame returned by \code{reestimate}
#' @param n11 Interim sample size of group 1
#' @param n12 Interim sample size of group 2
#' @param r Allocation ratio to group 1
#' @param rounding \code{'group'}, \code{'friede-kieser'} or \code{'total'}
#' @param restricted Logical. If \code{TRUE}, the planned sizes are lower bounds
#' @param N1.plan Planned sample size of group 1, or \code{NULL}
#' @param N2.plan Planned sample size of group 2, or \code{NULL}
#' @param N.min Lower bound on the final total, or \code{NULL}
#' @param N.max Upper bound on the final total, or \code{NULL}
#'
#' @return A list with integer vectors \code{N1} and \code{N2}
#'
#' @keywords internal
#' @noRd
final_sizes <- function(re, n11, n12, r, rounding, restricted, N1.plan, N2.plan,
                        N.min, N.max) {
  n1 <- n11 + n12
  if (rounding == 'group') {
    lo2 <- max(n12, if (restricted) N2.plan else 0,
               if (!is.null(N.min)) ceil_tol(N.min / (1 + r)) else 0)
    hi2 <- Inf
    if (!is.null(N.max)) {
      hi2 <- floor(N.max / (1 + r))
      while ((hi2 + 1) + ceiling(r * (hi2 + 1)) <= N.max) hi2 <- hi2 + 1
      while (hi2 > 0 && hi2 + ceiling(r * hi2) > N.max) hi2 <- hi2 - 1
    }
    if (hi2 < lo2) {
      stop('N.max is smaller than the smallest final sample size the design allows')
    }
    N2 <- pmax(lo2, pmin(hi2, re$N2.re))
    N1 <- pmax(n11, ceiling(r * N2))
    return(list(N1 = as.integer(N1), N2 = as.integer(N2)))
  }
  lo <- max(n1, if (restricted) N1.plan + N2.plan else 0,
            if (!is.null(N.min)) N.min else 0)
  hi <- if (!is.null(N.max)) N.max else Inf
  if (hi < lo) stop('N.max is smaller than the smallest final sample size the design allows')
  n.c <- pmin(hi, pmax(lo, re$n.raw))
  if (rounding == 'friede-kieser') {
    n2 <- pmax(ceil_tol(n.c - n1), 0)
    N2 <- n12 + ceil_tol(n2 / (1 + r))
    N1 <- n11 + ceil_tol(r * n2 / (1 + r))
  } else {
    N <- ceil_tol(n.c)
    N2 <- pmax(n12, floor(N / (1 + r)))
    N1 <- pmax(n11, N - N2)
  }
  list(N1 = as.integer(N1), N2 = as.integer(N2))
}
