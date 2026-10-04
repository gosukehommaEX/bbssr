#' Interim Cells and Final Sample Sizes of a Re-Estimation Design
#'
#' Internal helper that lists every interim outcome of a design described by
#' \code{bssr_map}, together with the distinct pairs of final sample sizes and the index of
#' the pair reached from each interim outcome. The result is passed to \code{bssr_power}
#' with the rejection regions of the distinct final sample sizes.
#'
#' @param map Data frame returned by \code{bssr_map}
#'
#' @return A list with the interim sizes \code{n11} and \code{n12}, the interim responder
#'   counts \code{x11} and \code{x12} of every interim outcome, the zero based index
#'   \code{rr.id} of its final sample sizes, the distinct final sizes \code{N1} and
#'   \code{N2}, and the corresponding second-stage sizes \code{n21} and \code{n22}
#'
#' @keywords internal
#' @noRd
bssr_setup <- function(map) {
  n11 <- attr(map, 'n11')
  n12 <- attr(map, 'n12')
  key <- paste(map$N1, map$N2, sep = '_')
  first <- which(!duplicated(key))
  id.s <- match(key, key[first]) - 1L
  # Interim responder counts, with the count of group 1 varying fastest
  x11 <- rep(0:n11, times = n12 + 1L)
  x12 <- rep(0:n12, each = n11 + 1L)
  list(n11 = n11, n12 = n12, x11 = as.integer(x11), x12 = as.integer(x12),
       rr.id = as.integer(id.s[x11 + x12 + 1L]),
       N1 = map$N1[first], N2 = map$N2[first],
       n21 = as.integer(map$N1[first] - n11), n22 = as.integer(map$N2[first] - n12))
}
