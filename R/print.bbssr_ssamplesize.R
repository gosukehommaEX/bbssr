#' Print Method for bbssr_ssamplesize Objects
#'
#' Prints the sample size of a stratified test returned by \code{\link{StratSampleSize}}.
#'
#' @param x An object of class \code{bbssr_ssamplesize}
#' @param digits Number of decimal places of the power. Default is 4
#' @param ... Further arguments, ignored
#'
#' @return \code{x}, invisibly
#'
#' @examples
#' ss <- StratSampleSize(p1 = c(0.5, 0.7), p2 = c(0.3, 0.5), frac = c(0.5, 0.5),
#'                       alpha = 0.025, tar.power = 0.8, Test = 'CMH',
#'                       method = 'standard')
#' print(ss)
#'
#' @export
print.bbssr_ssamplesize <- function(x, digits = 4, ...) {
  cat('Sample size of a stratified test\n\n')
  cat(sprintf('  Test             : %s\n', x$Test))
  cat(sprintf('  Alternative      : %s\n', x$alternative))
  cat(sprintf('  Method           : %s\n',
              if (attr(x, 'method') == 'exact') {
                sprintf('exact (search %s)', attr(x, 'search'))
              } else {
                'normal approximation'
              }))
  cat(sprintf('  Alpha            : %s, target power %s\n', format(x$alpha),
              format(x$tar.power)))
  cat(sprintf('  Normal approx.   : %s (unrounded total)\n',
              format(round(attr(x, 'n.raw'), 2))))
  alloc <- attr(x, 'alloc')
  random <- !is.null(alloc) && alloc != 'fixed'
  if (random) {
    cat(sprintf('  Allocation       : %s\n', strat_alloc_label(alloc, attr(x, 'block'),
                                                           attr(x, 'n.strata'))))
  }
  cat(sprintf('  Total size       : N = %d, power %s\n\n', x$N,
              format(round(x$Power, digits))))
  tab <- data.frame(Stratum = seq_along(x$frac), frac = round(x$frac, 4), r = x$r,
                    p1 = x$p1, p2 = x$p2)
  if (!random) {
    tab$N1 <- x$N1
    tab$N2 <- x$N2
  }
  print(tab, row.names = FALSE)
  invisible(x)
}
