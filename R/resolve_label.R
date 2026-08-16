#' Resolve an Optional Plot Element
#'
#' Internal helper implementing the convention shared by the plotting methods of the
#' package. An argument left at its default of \code{NULL} takes the value built from the
#' object, a single missing value removes the element from the plot, and anything else is
#' passed through unchanged.
#'
#' @param value Value supplied through the plotting method
#' @param default Value used when \code{value} is \code{NULL}
#'
#' @return \code{default} when \code{value} is \code{NULL}, \code{NULL} when \code{value}
#'   is a single missing value, and \code{value} otherwise
#'
#' @keywords internal
#' @noRd
resolve_label <- function(value, default) {
  if (is.null(value)) return(default)
  if (length(value) == 1 && is.na(value)) return(NULL)
  value
}
