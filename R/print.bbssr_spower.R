#' Print Method for bbssr_spower Objects
#'
#' Prints the exact power of a stratified test returned by \code{\link{StratPower}}.
#'
#' An object whose attributes were removed, for example by selecting some of its
#' columns, is printed as a data frame.
#'
#' @param x An object of class \code{bbssr_spower}
#' @param digits Number of decimal places. Default is 4
#' @param ... Further arguments, ignored
#'
#' @return \code{x}, invisibly
#'
#' @examples
#' pw <- StratPower(p1 = c(0.5, 0.7), p2 = c(0.3, 0.5), N1 = c(40, 40), N2 = c(40, 40),
#'                  alpha = 0.025, Test = 'CMH')
#' print(pw)
#'
#' @export
print.bbssr_spower <- function(x, digits = 4, ...) {
  if (is.null(attr(x, 'Test')) || !('Power' %in% names(x))) return(NextMethod())
  cat('Exact power of a stratified test\n\n')
  cat(sprintf('  Test             : %s\n', attr(x, 'Test')))
  cat(sprintf('  Alternative      : %s\n', attr(x, 'alternative')))
  if (attr(x, 'Test') %in% c('Strat-Fisher', 'Unstrat-Fisher') &&
      attr(x, 'alternative') == 'two.sided') {
    cat(sprintf('  Two-sided method : %s\n', attr(x, 'tsmethod')))
  }
  alloc <- attr(x, 'alloc')
  if (is.null(alloc) || alloc == 'fixed') {
    cat(sprintf('  Group sizes      : N1 = (%s), N2 = (%s)\n',
                paste(attr(x, 'N1'), collapse = ', '),
                paste(attr(x, 'N2'), collapse = ', ')))
  } else {
    cat(sprintf('  Allocation       : %s\n', strat_alloc_label(alloc, attr(x, 'block'),
                                                           attr(x, 'n.strata'))))
    cat(sprintf('  Total size       : N = %d, fractions (%s), ratio %s\n', attr(x, 'N'),
                paste(format(round(attr(x, 'frac'), 4)), collapse = ', '),
                paste(format(unique(attr(x, 'r'))), collapse = ', ')))
    cat(sprintf('  Left out         : configurations of total probability %s\n',
                format(signif(attr(x, 'skipped'), 2))))
  }
  cat(sprintf('  Alpha            : %s\n\n', format(attr(x, 'alpha'))))
  tab <- data.frame(lapply(unclass(x), identity), check.names = FALSE)
  tab$Power <- round(tab$Power, digits)
  print(tab, row.names = FALSE)
  invisible(x)
}
