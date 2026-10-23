#' Print Method for bbssr_ptie Objects
#'
#' Prints the settings and the largest type I error rates of a matched-pairs design with
#' blinded sample size re-estimation and of the fixed design.
#'
#' @param x An object of class \code{bbssr_ptie} returned by
#'   \code{\link{PairedTypeIErrorBSSR}}
#' @param digits Number of significant digits of the rates. Default is 4
#' @param ... Further arguments, ignored
#'
#' @return \code{x}, invisibly
#'
#' @examples
#' tie <- PairedTypeIErrorBSSR(Delta.A = 0.2, N = 40, n.interim = 20, alpha = 0.05,
#'                             tar.power = 0.8, Test = 'McNemar',
#'                             alternative = 'two.sided', ss.method = 'standard',
#'                             psi = seq(0.1, 0.9, by = 0.1))
#' print(tie)
#'
#' @export
print.bbssr_ptie <- function(x, digits = 4, ...) {
  m <- attr(x, 'max')
  if (is.null(m) || !all(c('psi', 'TIE.BSSR', 'TIE.TRAD') %in% names(x))) {
    return(NextMethod())
  }
  cat('Type I error rate of blinded sample size re-estimation for paired binary data\n\n')
  cat(sprintf('  Test            : %s\n', attr(x, 'Test')))
  cat(sprintf('  Alternative     : %s\n', attr(x, 'alternative')))
  cat(sprintf('  Pairs           : N = %d planned, %d at the interim analysis\n',
              attr(x, 'N'), attr(x, 'n.interim')))
  cat(sprintf('  Assumed effect  : %s\n', format(attr(x, 'Delta.A'))))
  margin <- format_margin(attr(x, 'margin'), 'RD')
  if (!is.null(margin)) {
    cat(sprintf('  Margin          : %s\n', margin))
  }
  cat(sprintf('  Nominal level   : %s\n', format(attr(x, 'alpha'))))
  if ('boundary' %in% names(x)) {
    cat(sprintf(paste('  Grid            : %d values of psi on the lower boundary and',
                      '%d on the upper boundary\n'),
                sum(x$boundary == 'lower'), sum(x$boundary == 'upper')))
  } else {
    cat(sprintf('  Grid            : %d values of psi in [%s, %s]\n', nrow(x),
                format(min(x$psi)), format(max(x$psi))))
  }
  cat(sprintf('  Maximum         : %s\n\n', maximize_label(x)))
  cat('Largest type I error rate\n')
  tab <- data.frame(Design = c('BSSR', 'Fixed sample'), psi = signif(m$psi, 6),
                    TIE = signif(m$TIE, digits))
  if (!is.null(m$boundary)) tab$Boundary <- m$boundary
  print(tab, row.names = FALSE)
  invisible(x)
}
