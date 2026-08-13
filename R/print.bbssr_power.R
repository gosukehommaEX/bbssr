#' Print Exact Power Results
#'
#' @param x An object of class \code{bbssr_power} returned by \code{\link{BinaryPower}}
#' @param digits Number of significant digits used for the power
#' @param ... Further arguments, currently ignored
#'
#' @return The object \code{x}, returned invisibly
#'
#' @examples
#' pw <- BinaryPower(p1 = 0.5, p2 = 0.2, N1 = 5, N2 = 5, alpha = 0.025, Test = 'Chisq')
#' print(pw)
#'
#' @export
print.bbssr_power <- function(x, digits = 4, ...) {
  if (!all(c('p1', 'p2', 'N1', 'N2', 'alpha', 'Test', 'alternative', 'Power') %in%
           names(x))) return(NextMethod())
  cat('Exact power for a two-arm trial with a binary endpoint\n\n')
  cat(sprintf('  Test         : %s\n', x$Test[1]))
  cat(sprintf('  Alternative  : %s\n', x$alternative[1]))
  cat(sprintf('  Sample sizes : N1 = %d, N2 = %d\n', x$N1[1], x$N2[1]))
  cat(sprintf('  Alpha        : %s\n\n', format(x$alpha[1])))
  tab <- data.frame(p1 = x$p1, p2 = x$p2, Power = round(x$Power, digits))
  print(tab, row.names = FALSE)
  invisible(x)
}
