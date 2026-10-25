#' Print Method for bbssr_sbssr Objects
#'
#' Prints the re-estimated group sizes of a stratified trial returned by
#' \code{\link{StratBSSR}}.
#'
#' An object whose attributes were removed, for example by selecting some of its
#' columns, is printed as a data frame.
#'
#' @param x An object of class \code{bbssr_sbssr}
#' @param digits Number of decimal places. Default is 4
#' @param ... Further arguments, ignored
#'
#' @return \code{x}, invisibly
#'
#' @examples
#' re <- StratBSSR(n1 = c(15, 15), n2 = c(15, 15), S = c(12, 18), Delta.A = 0.2, r = 1,
#'                 alpha = 0.025, tar.power = 0.8, Test = 'CMH', N1 = c(30, 30),
#'                 N2 = c(30, 30), ss.method = 'standard')
#' print(re)
#'
#' @export
print.bbssr_sbssr <- function(x, digits = 4, ...) {
  if (is.null(attr(x, 'ss.method')) || !all(c('N1.final', 'N2.final') %in% names(x))) {
    return(NextMethod())
  }
  rule <- c(exact = 'exact sample size', standard = 'normal approximation')[[
    attr(x, 'ss.method')]]
  cat('Blinded sample size re-estimation for a stratified trial\n\n')
  cat(sprintf('  Test             : %s\n', attr(x, 'Test')))
  cat(sprintf('  Alternative      : %s\n', attr(x, 'alternative')))
  cat(sprintf('  Re-estimation    : %s (%s%s)\n', rule,
              if (isTRUE(attr(x, 'restricted'))) 'restricted' else 'unrestricted',
              if (identical(attr(x, 'ss.strata'), 'pooled')) ', strata pooled' else ''))
  cat(sprintf('  Effect           : %s, assumed %s\n', attr(x, 'effect'),
              format(attr(x, 'Delta.A'))))
  alloc <- attr(x, 'alloc')
  random <- !is.null(alloc) && alloc != 'fixed'
  if (random) {
    cat(sprintf('  Allocation       : %s\n', strat_alloc_label(alloc, attr(x, 'block'))))
  }
  cat(sprintf('  Final total      : %d, power %s at the planning probabilities\n\n',
              attr(x, 'N.final'), format(round(attr(x, 'Power'), digits))))
  tab <- data.frame(lapply(unclass(x), identity), check.names = FALSE)
  # With a random allocation only the total is re-estimated
  if (random) tab[c('N1.final', 'N2.final', 'n1.stage2', 'n2.stage2')] <- NULL
  for (v in c('hat.p', 'hat.p1', 'hat.p2')) tab[[v]] <- round(tab[[v]], digits)
  print(tab, row.names = FALSE)
  invisible(x)
}
