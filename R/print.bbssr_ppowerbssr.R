#' Print Method for bbssr_ppowerbssr Objects
#'
#' Prints the power and the expected final number of pairs of a matched-pairs design with
#' blinded sample size re-estimation.
#'
#' An object whose attributes were removed, for example by selecting some of its
#' columns, is printed as a data frame.
#'
#' @param x An object of class \code{bbssr_ppowerbssr} returned by
#'   \code{\link{PairedPowerBSSR}}
#' @param digits Number of decimal places. Default is 4
#' @param ... Further arguments, ignored
#'
#' @return \code{x}, invisibly
#'
#' @examples
#' res <- PairedPowerBSSR(psi = c(0.2, 0.3), Delta.T = 0.15, Delta.A = 0.15, N = 72,
#'                        n.interim = 30, alpha = 0.05, tar.power = 0.8,
#'                        Test = 'McNemar', alternative = 'two.sided',
#'                        ss.method = 'standard')
#' print(res)
#'
#' @export
print.bbssr_ppowerbssr <- function(x, digits = 4, ...) {
  if (is.null(attr(x, 'ss.method')) ||
      !all(c('psi', 'p10', 'p01', 'power.BSSR', 'power.TRAD', 'E.N') %in% names(x))) {
    return(NextMethod())
  }
  rule <- c(exact = 'exact sample size', standard = 'normal approximation',
            conditional.power = 'conditional power')[[attr(x, 'ss.method')]]
  cat('Blinded sample size re-estimation for paired binary data\n\n')
  cat(sprintf('  Test             : %s\n', attr(x, 'Test')))
  cat(sprintf('  Alternative      : %s\n', attr(x, 'alternative')))
  cat(sprintf('  Re-estimation    : %s (%s)\n', rule,
              if (isTRUE(attr(x, 'restricted'))) 'restricted' else 'unrestricted'))
  cat(sprintf('  Pairs            : N = %d planned, %d at the interim analysis\n',
              attr(x, 'N'), attr(x, 'n.interim')))
  cat(sprintf('  Effect           : assumed %s, true %s\n', format(attr(x, 'Delta.A')),
              format(attr(x, 'Delta.T'))))
  margin <- format_margin(attr(x, 'margin'), 'RD')
  if (!is.null(margin)) {
    cat(sprintf('  Margin           : %s\n', margin))
  }
  cat(sprintf('  Alpha            : %s, target power %s\n\n', format(attr(x, 'alpha')),
              format(attr(x, 'tar.power'))))
  tab <- data.frame(psi = x$psi, p10 = x$p10, p01 = x$p01,
                    power.BSSR = round(x$power.BSSR, digits),
                    power.TRAD = round(x$power.TRAD, digits), E.N = round(x$E.N, 2))
  print(tab, row.names = FALSE)
  invisible(x)
}
