#' Print Method for bbssr_ppower Objects
#'
#' Prints the exact power of a test for paired binary data.
#'
#' @param x An object of class \code{bbssr_ppower} returned by \code{\link{PairedPower}}
#' @param digits Number of decimal places of the power. Default is 4
#' @param ... Further arguments, ignored
#'
#' @return \code{x}, invisibly
#'
#' @examples
#' pw <- PairedPower(p10 = c(0.3, 0.35), p01 = 0.1, N = 40, alpha = 0.025,
#'                   Test = 'McNemar-exact')
#' print(pw)
#'
#' @export
print.bbssr_ppower <- function(x, digits = 4, ...) {
  if (!all(c('p10', 'p01', 'N', 'alpha', 'Test', 'alternative', 'Power') %in% names(x))) {
    return(NextMethod())
  }
  cat('Exact power for paired binary data\n\n')
  cat(sprintf('  Test        : %s\n', x$Test[1]))
  cat(sprintf('  Alternative : %s\n', x$alternative[1]))
  cat(sprintf('  Pairs       : N = %d\n', x$N[1]))
  margin <- format_margin(attr(x, 'margin'), 'RD')
  if (!is.null(margin)) {
    cat(sprintf('  Margin      : %s\n', margin))
  }
  cat(sprintf('  Alpha       : %s\n\n', format(x$alpha[1])))
  tab <- data.frame(p10 = x$p10, p01 = x$p01, Power = round(x$Power, digits))
  print(tab, row.names = FALSE)
  invisible(x)
}
