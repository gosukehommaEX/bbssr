#' Print a Sample Size Re-estimation from Interim Data
#'
#' @param x An object of class \code{bbssr_bssr} returned by \code{\link{BinaryBSSR}}
#' @param digits Number of significant digits used for the proportions and the power
#' @param ... Further arguments, currently ignored
#'
#' @return The object \code{x}, returned invisibly
#'
#' @examples
#' res <- BinaryBSSR(n1 = 20, n2 = 20, S = 11, Delta.A = 0.3, r = 1,
#'                   alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
#' print(res)
#'
#' @export
print.bbssr_bssr <- function(x, digits = 4, ...) {
  if (!all(c('n1', 'n2', 'n', 'S', 'hat.p', 'hat.p1', 'hat.p2', 'N1.re', 'N2.re',
             'N.re', 'n1.stage2', 'n2.stage2', 'n.stage2', 'N1.final', 'N2.final',
             'N.final', 'Power') %in% names(x))) return(NextMethod())
  cat('Blinded sample size re-estimation from interim data\n\n')
  cat(sprintf('  Test             : %s\n', attr(x, 'Test')))
  cat(sprintf('  Alternative      : %s\n', attr(x, 'alternative')))
  cat(sprintf('  Design rule      : %s\n',
              if (isTRUE(attr(x, 'restricted'))) 'restricted' else 'unrestricted'))
  cat(sprintf('  Assumed effect   : %s\n', format(attr(x, 'Delta.A'))))
  cat(sprintf('  Alpha            : %s, target power %s\n\n',
              format(attr(x, 'alpha')), format(attr(x, 'tar.power'))))
  cat('Interim data\n')
  cat(sprintf('  Patients         : n1 = %d, n2 = %d, total n = %d\n', x$n1, x$n2, x$n))
  cat(sprintf('  Responders       : S = %d, blinded pooled rate = %s\n',
              x$S, format(round(x$hat.p, digits))))
  cat(sprintf('  Recovered rates  : hat.p1 = %s, hat.p2 = %s\n\n',
              format(round(x$hat.p1, digits)), format(round(x$hat.p2, digits))))
  cat('Re-estimation\n')
  cat(sprintf('  Required total   : N1 = %d, N2 = %d, total N = %d\n',
              x$N1.re, x$N2.re, x$N.re))
  cat(sprintf('  Still to enrol   : group 1 = %d, group 2 = %d, total = %d\n',
              x$n1.stage2, x$n2.stage2, x$n.stage2))
  cat(sprintf('  Final size       : N1 = %d, N2 = %d, total N = %d\n',
              x$N1.final, x$N2.final, x$N.final))
  cat(sprintf('  Power at final N : %s\n', format(round(x$Power, digits))))
  invisible(x)
}
