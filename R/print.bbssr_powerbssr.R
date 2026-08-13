#' Print the Operating Characteristics of a BSSR Design
#'
#' @param x An object of class \code{bbssr_powerbssr} returned by
#'   \code{\link{BinaryPowerBSSR}}
#' @param digits Number of significant digits used for the power
#' @param ... Further arguments, currently ignored
#'
#' @return The object \code{x}, returned invisibly
#'
#' @examples
#' res <- BinaryPowerBSSR(
#'   asmd.p1 = 0.6, asmd.p2 = 0.3, p = 0.45,
#'   Delta.A = 0.3, Delta.T = 0.3,
#'   N1 = 5, N2 = 5, omega = 0.5, r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
#' )
#' print(res)
#'
#' @export
print.bbssr_powerbssr <- function(x, digits = 4, ...) {
  if (!all(c('p1', 'p2', 'p', 'power.BSSR', 'power.TRAD', 'E.N') %in% names(x))) {
    return(NextMethod())
  }
  cat('Blinded sample size re-estimation for a binary endpoint\n\n')
  cat(sprintf('  Test            : %s\n', attr(x, 'Test')))
  cat(sprintf('  Alternative     : %s\n', attr(x, 'alternative')))
  cat(sprintf('  Design rule     : %s\n',
              if (isTRUE(attr(x, 'restricted'))) 'restricted' else 'unrestricted'))
  cat(sprintf('  Initial size    : N1 = %s, N2 = %s\n',
              format(attr(x, 'N1')), format(attr(x, 'N2'))))
  cat(sprintf('  Interim fraction: %s\n', format(attr(x, 'omega'))))
  cat(sprintf('  Treatment effect: assumed %s, true %s\n',
              format(attr(x, 'Delta.A')), format(attr(x, 'Delta.T'))))
  cat(sprintf('  Alpha           : %s, target power %s\n\n',
              format(attr(x, 'alpha')), format(attr(x, 'tar.power'))))
  tab <- data.frame(
    p = x$p, p1 = x$p1, p2 = x$p2,
    power.BSSR = round(x$power.BSSR, digits),
    power.TRAD = round(x$power.TRAD, digits),
    E.N = round(x$E.N, 1)
  )
  print(tab, row.names = FALSE)
  invisible(x)
}
