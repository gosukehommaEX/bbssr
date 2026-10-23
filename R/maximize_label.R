#' Description of How the Largest Type I Error Rate Was Located
#'
#' Internal helper of the print methods for \code{bbssr_tie}, \code{bbssr_ptie},
#' \code{bbssr_stie} and \code{bbssr_alphaadj} objects, describing the attribute
#' \code{maximize} and, for a certified maximum, the attribute \code{interval}. The
#' interval is one of the nuisance parameter named by the attribute \code{nuisance},
#' \code{theta} when it is absent. For a stratified design it is a matrix with one row per
#' stratum, named \code{theta.S1} and \code{theta.S2}, and describes a rectangle.
#'
#' @param x Object returned by \code{BinaryTypeIErrorBSSR}, \code{BinaryAlphaAdjBSSR},
#'   \code{PairedTypeIErrorBSSR}, \code{PairedAlphaAdjBSSR},
#'   \code{StratTypeIErrorBSSR} or \code{StratAlphaAdjBSSR}
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
  if (maximize == 'certified' && is.matrix(interval) &&
      identical(rownames(interval), c('theta.S1', 'theta.S2'))) {
    return(sprintf('certified over (theta.S1, theta.S2) in [%s, %s] x [%s, %s]',
                   format(signif(interval[1, 1], 6)), format(signif(interval[1, 2], 6)),
                   format(signif(interval[2, 1], 6)), format(signif(interval[2, 2], 6))))
  }
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
