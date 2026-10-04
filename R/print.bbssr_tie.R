#' Print Method for bbssr_tie Objects
#'
#' Prints the design settings and the largest type I error rate of a design with blinded
#' sample size re-estimation and of the corresponding fixed-sample design.
#'
#' @param x An object of class \code{bbssr_tie}
#' @param digits Number of significant digits of the type I error rates. Default is 4
#' @param ... Further arguments, ignored
#'
#' @return \code{x}, invisibly
#'
#' @examples
#' tie <- BinaryTypeIErrorBSSR(
#'   Delta.A = 0.3, N1 = 20, N2 = 20, omega = 0.5, r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Chisq', ss.method = 'standard',
#'   theta = seq(0.05, 0.95, by = 0.05)
#' )
#' print(tie)
#'
#' @export
print.bbssr_tie <- function(x, digits = 4, ...) {
  m <- attr(x, 'max')
  if (is.null(m) || !all(c('theta', 'TIE.BSSR', 'TIE.TRAD') %in% names(x))) {
    return(NextMethod())
  }
  cat('Type I error rate of blinded sample size re-estimation for a binary endpoint\n\n')
  cat(sprintf('  Test            : %s\n', attr(x, 'Test')))
  cat(sprintf('  Alternative     : %s\n', attr(x, 'alternative')))
  cat(sprintf('  Initial size    : N1 = %s, N2 = %s\n',
              format(attr(x, 'N1')), format(attr(x, 'N2'))))
  cat(sprintf('  Interim size    : n1 = %d, n2 = %d\n',
              attr(x, 'n1.interim'), attr(x, 'n2.interim')))
  cat(sprintf('  Assumed effect  : %s (%s)\n', format(attr(x, 'Delta.A')), attr(x, 'effect')))
  cat(sprintf('  Nominal level   : %s\n', format(attr(x, 'alpha'))))
  cat(sprintf('  Grid            : %d values of theta in [%s, %s]%s\n\n', nrow(x),
              format(min(x$theta)), format(max(x$theta)),
              if (isTRUE(attr(x, 'refine'))) ', maxima refined' else ''))
  cat('Largest type I error rate\n')
  tab <- data.frame(Design = c('BSSR', 'Fixed sample'),
                    theta = signif(m$theta, 6), TIE = signif(m$TIE, digits))
  print(tab, row.names = FALSE)
  invisible(x)
}
