#' Bernstein Coefficients on a Subinterval
#'
#' Internal helper returning the Bernstein coefficients, on the subinterval
#' \code{[t0, t1]} of \code{[0, 1]} mapped onto \code{[0, 1]}, of a polynomial given by its
#' Bernstein coefficients on \code{[0, 1]}. The coefficients are obtained by two
#' subdivisions of de Casteljau, first at \code{t1}, keeping the left part, and then at
#' \code{t0 / t1}, keeping the right part. Every step is a convex combination, so the
#' computation is numerically stable. When \code{t0 = t1} the polynomial is constant on
#' the subinterval and all coefficients equal its value there.
#'
#' @param coef Bernstein coefficients on \code{[0, 1]}
#' @param t0,t1 End points of the subinterval, \code{0 <= t0 <= t1 <= 1}
#'
#' @return A numeric vector of the length of \code{coef}
#'
#' @keywords internal
#' @noRd
bern_restrict <- function(coef, t0, t1) {
  # Control points of the two parts of a split at t
  split <- function(b, t) {
    n <- length(b) - 1L
    left <- numeric(n + 1L)
    right <- numeric(n + 1L)
    left[1] <- b[1]
    right[n + 1L] <- b[n + 1L]
    for (j in seq_len(n)) {
      b <- (1 - t) * b[-length(b)] + t * b[-1]
      left[j + 1L] <- b[1]
      right[n + 1L - j] <- b[length(b)]
    }
    list(left = left, right = right)
  }
  if (t1 <= t0) return(rep(split(coef, t0)$left[length(coef)], length(coef)))
  if (t1 < 1) coef <- split(coef, t1)$left
  if (t0 > 0) coef <- split(coef, t0 / t1)$right
  coef
}
