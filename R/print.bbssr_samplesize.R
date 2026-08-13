#' Print a Sample Size Calculation
#'
#' @param x An object of class \code{bbssr_samplesize} returned by
#'   \code{\link{BinarySampleSize}}
#' @param digits Number of significant digits used for the power
#' @param ... Further arguments, currently ignored
#'
#' @return The object \code{x}, returned invisibly
#'
#' @examples
#' ss <- BinarySampleSize(p1 = 0.4, p2 = 0.2, r = 1, alpha = 0.025,
#'                        tar.power = 0.8, Test = 'Chisq')
#' print(ss)
#'
#' @export
print.bbssr_samplesize <- function(x, digits = 4, ...) {
  if (!all(c('p1', 'p2', 'r', 'alpha', 'tar.power', 'Test', 'alternative', 'Power',
             'N1', 'N2', 'N') %in% names(x))) return(NextMethod())
  cat('Sample size for a two-arm trial with a binary endpoint\n\n')
  cat(sprintf('  Test             : %s\n', x$Test[1]))
  cat(sprintf('  Alternative      : %s\n', x$alternative[1]))
  cat(sprintf('  Response rates   : p1 = %s, p2 = %s\n', format(x$p1[1]), format(x$p2[1])))
  cat(sprintf('  Allocation ratio : %s to 1\n', format(x$r[1])))
  cat(sprintf('  Alpha            : %s\n', format(x$alpha[1])))
  cat(sprintf('  Target power     : %s\n\n', format(x$tar.power[1])))
  cat(sprintf('  Required sample size: N1 = %d, N2 = %d, total N = %d\n',
              x$N1[1], x$N2[1], x$N[1]))
  cat(sprintf('  Attained power      : %s\n', format(round(x$Power[1], digits))))
  invisible(x)
}
