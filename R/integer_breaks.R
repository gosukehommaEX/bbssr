#' Axis Breaks Restricted to Integers
#'
#' Internal helper returning a break function for a continuous scale whose values are
#' counts of patients or of responders. Every integer of the range is placed when the range
#' is short enough, and a subset chosen by \code{pretty} is placed otherwise.
#'
#' A scale drawn with \code{geom_tile} extends half a tile beyond the outermost value and
#' is then expanded further, so the limits passed to a break function reach below zero and
#' above the largest count. Supplying \code{range} keeps the breaks inside the counts that
#' can actually occur. The limits are accepted in either order, so the function can be used
#' with a reversed scale.
#'
#' @param range Optional numeric vector of length two giving the smallest and the largest
#'   value that a break may take. The limits of the scale are used when it is \code{NULL}
#' @param max.breaks Largest number of integer breaks to place before thinning
#'
#' @return A function suitable for the \code{breaks} argument of a continuous scale
#'
#' @keywords internal
#' @noRd
integer_breaks <- function(range = NULL, max.breaks = 21) {
  function(limits) {
    if (is.null(range)) {
      lower <- ceiling(min(limits))
      upper <- floor(max(limits))
    } else {
      lower <- ceiling(min(range))
      upper <- floor(max(range))
    }
    if (!is.finite(lower) || !is.finite(upper) || upper < lower) return(numeric(0))
    if (upper - lower + 1 <= max.breaks) return(seq(lower, upper))
    candidates <- unique(round(pretty(c(lower, upper))))
    candidates[candidates >= lower & candidates <= upper]
  }
}
