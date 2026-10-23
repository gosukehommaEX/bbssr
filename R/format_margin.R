#' Label of a Non-Inferiority Margin for Printing
#'
#' Internal helper of the print methods returning the margin stored with a result as a
#' character string, or \code{NULL} when the result has no margin. A margin on the scale
#' of the risk difference is shown as its value, and is absent when it is 0. A margin on
#' the scale of the risk ratio is followed by \code{(risk ratio)}. Results without the
#' attribute \code{margin.scale} have a margin on the scale of the risk difference. The
#' two margins of a test of equivalence are shown as the equivalence region, from
#' \code{-lower} to \code{upper} for the risk difference and from \code{lower} to
#' \code{upper} for the risk ratio.
#'
#' @param margin The attribute \code{margin} of the result, or \code{NULL}
#' @param margin.scale The attribute \code{margin.scale} of the result, or \code{NULL}
#'
#' @return A character string, or \code{NULL}
#'
#' @keywords internal
#' @noRd
format_margin <- function(margin, margin.scale) {
  if (is.null(margin)) return(NULL)
  if (length(margin) == 2) {
    if (identical(margin.scale, 'RR')) {
      return(sprintf('%s to %s (risk ratio)', format(margin[1]), format(margin[2])))
    }
    return(sprintf('%s to %s (risk difference)', format(-margin[1]), format(margin[2])))
  }
  if (identical(margin.scale, 'RR')) return(paste(format(margin), '(risk ratio)'))
  if (margin == 0) return(NULL)
  format(margin)
}
