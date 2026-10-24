#' Print Method for bbssr_srr Objects
#'
#' Prints the critical values of a stratified test returned by \code{\link{StratRR}}.
#'
#' An object whose attributes were removed, for example by selecting some of its
#' columns, is printed as a data frame.
#'
#' @param x An object of class \code{bbssr_srr}
#' @param n Number of rows of the table to print. Default is 20
#' @param ... Further arguments, ignored
#'
#' @return \code{x}, invisibly
#'
#' @examples
#' rr <- StratRR(N1 = c(10, 12), N2 = c(10, 12), alpha = 0.025, Test = 'CMH')
#' print(rr, n = 10)
#'
#' @export
print.bbssr_srr <- function(x, n = 20, ...) {
  if (is.null(attr(x, 'Test')) ||
      !all(c('S.min', 'S.max', 'lower', 'upper') %in% names(x))) {
    return(NextMethod())
  }
  cat('Rejection region of a stratified test\n\n')
  cat(sprintf('  Test             : %s\n', attr(x, 'Test')))
  cat(sprintf('  Alternative      : %s\n', attr(x, 'alternative')))
  if (attr(x, 'Test') %in% c('Strat-Fisher', 'Unstrat-Fisher') &&
      attr(x, 'alternative') == 'two.sided') {
    cat(sprintf('  Two-sided method : %s\n', attr(x, 'tsmethod')))
  }
  cat(sprintf('  Group sizes      : N1 = (%s), N2 = (%s)\n',
              paste(attr(x, 'N1'), collapse = ', '), paste(attr(x, 'N2'), collapse = ', ')))
  cat(sprintf('  Alpha            : %s\n', format(attr(x, 'alpha'))))
  rej <- !is.na(x$lower) | !is.na(x$upper)
  cat(sprintf('  Combinations     : %d, of which %d have a rejected value of S\n\n',
              nrow(x), sum(rej)))
  tab <- data.frame(lapply(unclass(x), identity), check.names = FALSE)
  tab <- tab[rej, , drop = FALSE]
  if (nrow(tab) > 0) {
    print(utils::head(tab, n), row.names = FALSE)
    if (nrow(tab) > n) cat(sprintf('  ... %d more rows\n', nrow(tab) - n))
  }
  invisible(x)
}
