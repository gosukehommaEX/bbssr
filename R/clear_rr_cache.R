#' Empty the Store of p-Values
#'
#' Internal helper that removes every p-value matrix kept by \code{get_pvalue}, from which
#' the rejection regions are derived, and resets the count of stored cells.
#'
#' @return \code{NULL}, invisibly
#'
#' @keywords internal
#' @noRd
clear_rr_cache <- function() {
  .bbssr_cache$pv <- new.env(parent = emptyenv())
  .bbssr_cache$cells <- 0
  invisible(NULL)
}
