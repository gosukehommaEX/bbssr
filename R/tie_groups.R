#' Identify Groups of Tied Values in a Monotone Sequence
#'
#' Internal helper. Given a vector that has already been sorted, either increasing or
#' decreasing, the function returns the position of the first and the last element of the
#' tie group that each element belongs to. Two adjacent values are treated as tied when
#' their absolute difference is below a relative tolerance.
#'
#' A relative tolerance is used rather than the absolute tolerance of \code{fpCompare},
#' because the sequences passed to this function are often exact test p-values that can be
#' smaller than 1e-12, and an absolute tolerance of 1.5e-08 would merge values that are
#' genuinely distinct.
#'
#' @param x A numeric vector sorted in either increasing or decreasing order
#'
#' @return An integer matrix with \code{length(x)} rows and two columns named
#'   \code{first} and \code{last}
#'
#' @keywords internal
#' @noRd
tie_groups <- function(x) {
  n <- length(x)
  if (n == 0) {
    return(matrix(integer(0), nrow = 0, ncol = 2, dimnames = list(NULL, c('first', 'last'))))
  }
  if (n == 1) {
    return(matrix(1L, nrow = 1, ncol = 2, dimnames = list(NULL, c('first', 'last'))))
  }
  tol <- 1e-10
  ref <- pmax(abs(x[-n]), abs(x[-1]))
  new.grp <- abs(diff(x)) > tol * ref
  grp.end <- c(which(new.grp), n)
  grp.size <- diff(c(0L, grp.end))
  grp.start <- grp.end - grp.size + 1L
  cbind(
    first = rep(as.integer(grp.start), times = grp.size),
    last  = rep(as.integer(grp.end),   times = grp.size)
  )
}
