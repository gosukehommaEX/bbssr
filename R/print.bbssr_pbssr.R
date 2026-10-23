#' Print Method for bbssr_pbssr Objects
#'
#' Prints the result of a blinded sample size re-estimation from paired interim data.
#'
#' An object whose attributes were removed, for example by selecting some of its
#' columns, is printed as a data frame.
#'
#' @param x An object of class \code{bbssr_pbssr} returned by \code{\link{PairedBSSR}}
#' @param digits Number of decimal places. Default is 4
#' @param ... Further arguments, ignored
#'
#' @return \code{x}, invisibly
#'
#' @examples
#' res <- PairedBSSR(n1 = 25, D1 = 5, Delta.A = 0.15, alpha = 0.05, tar.power = 0.8,
#'                   Test = 'McNemar', alternative = 'two.sided', ss.method = 'standard')
#' print(res)
#'
#' @export
print.bbssr_pbssr <- function(x, digits = 4, ...) {
  if (is.null(attr(x, 'ss.method')) ||
      !all(c('n1', 'D1', 'psi.hat', 'psi.plan', 'p10.plan', 'p01.plan', 'N.re',
             'n.stage2', 'N.final', 'Power', 'CP') %in% names(x))) return(NextMethod())
  rule <- c(exact = 'exact sample size', standard = 'normal approximation',
            conditional.power = 'conditional power')[[attr(x, 'ss.method')]]
  cat('Blinded sample size re-estimation from paired interim data\n\n')
  cat(sprintf('  Test             : %s\n', attr(x, 'Test')))
  cat(sprintf('  Alternative      : %s\n', attr(x, 'alternative')))
  cat(sprintf('  Re-estimation    : %s (%s)\n', rule,
              if (isTRUE(attr(x, 'restricted'))) 'restricted' else 'unrestricted'))
  cat(sprintf('  Assumed effect   : %s\n', format(attr(x, 'Delta.A'))))
  margin <- format_margin(attr(x, 'margin'), 'RD')
  if (!is.null(margin)) {
    cat(sprintf('  Margin           : %s\n', margin))
  }
  cat(sprintf('  Alpha            : %s, target power %s\n\n', format(attr(x, 'alpha')),
              format(attr(x, 'tar.power'))))
  cat('Interim data\n')
  cat(sprintf('  Pairs            : n1 = %d, discordant D1 = %d\n', x$n1, x$D1))
  cat(sprintf('  Discordance      : estimate %s, used %s\n',
              format(round(x$psi.hat, digits)), format(round(x$psi.plan, digits))))
  cat(sprintf('  Planning values  : p10 = %s, p01 = %s\n\n',
              format(round(x$p10.plan, digits)), format(round(x$p01.plan, digits))))
  cat('Re-estimation\n')
  cat(sprintf('  Re-estimated N   : %d\n', x$N.re))
  cat(sprintf('  Still to enrol   : %d pairs\n', x$n.stage2))
  cat(sprintf('  Final N          : %d\n', x$N.final))
  cat(sprintf('  Power at final N : %s\n', format(round(x$Power, digits))))
  cat(sprintf('  Conditional power: %s\n', format(round(x$CP, digits))))
  invisible(x)
}
