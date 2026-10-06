#' Description of How the Largest Type I Error Rate Was Located
#'
#' Internal helper of the print methods for \code{bbssr_tie} and \code{bbssr_alphaadj}
#' objects, describing the attribute \code{maximize} and, for a certified maximum, the
#' attribute \code{interval}.
#'
#' @param x Object returned by \code{BinaryTypeIErrorBSSR} or \code{BinaryAlphaAdjBSSR}
#'
#' @return A character string
#'
#' @keywords internal
#' @noRd
maximize_label <- function(x) {
  maximize <- attr(x, 'maximize')
  interval <- attr(x, 'interval')
  if (is.null(maximize)) return('not recorded')
  if (maximize == 'certified' && length(interval) == 2) {
    return(sprintf('certified over theta in [%s, %s]', format(signif(interval[1], 6)),
                   format(signif(interval[2], 6))))
  }
  switch(maximize,
         certified = 'certified',
         refined = 'largest values on the grid refined',
         grid = 'largest value on the grid',
         maximize)
}
