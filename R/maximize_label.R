#' Description of How the Largest Type I Error Rate Was Located
#'
#' Internal helper of the print methods for \code{bbssr_tie}, \code{bbssr_ptie} and
#' \code{bbssr_alphaadj} objects, describing the attribute \code{maximize} and, for a
#' certified maximum, the attribute \code{interval}. The interval is one of the nuisance
#' parameter named by the attribute \code{nuisance}, \code{theta} when it is absent.
#'
#' @param x Object returned by \code{BinaryTypeIErrorBSSR}, \code{BinaryAlphaAdjBSSR},
#'   \code{PairedTypeIErrorBSSR} or \code{PairedAlphaAdjBSSR}
#'
#' @return A character string
#'
#' @keywords internal
#' @noRd
maximize_label <- function(x) {
  maximize <- attr(x, 'maximize')
  interval <- attr(x, 'interval')
  if (is.null(maximize)) return('not recorded')
  var <- attr(x, 'nuisance')
  if (is.null(var)) var <- 'theta'
  if (maximize == 'certified' && is.matrix(interval)) {
    return(sprintf(paste('certified over %s in [%s, %s] on the lower boundary and',
                         '[%s, %s] on the upper boundary'), var,
                   format(signif(interval[1, 1], 6)), format(signif(interval[1, 2], 6)),
                   format(signif(interval[2, 1], 6)), format(signif(interval[2, 2], 6))))
  }
  if (maximize == 'certified' && length(interval) == 2) {
    return(sprintf('certified over %s in [%s, %s]', var, format(signif(interval[1], 6)),
                   format(signif(interval[2], 6))))
  }
  switch(maximize,
         certified = 'certified',
         refined = 'largest values on the grid refined',
         grid = 'largest value on the grid',
         maximize)
}
