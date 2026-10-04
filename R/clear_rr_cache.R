#' Empty the Store of Rejection Regions
#'
#' Internal helper that removes every rejection region kept by \code{get_rr} and resets
#' the count of stored cells.
#'
#' @return \code{NULL}, invisibly
#'
#' @keywords internal
#' @noRd
clear_rr_cache <- function() {
  .bbssr_cache$rr <- new.env(parent = emptyenv())
  .bbssr_cache$cells <- 0
  invisible(NULL)
}
