#' Print Method for bbssr_spowerbssr Objects
#'
#' Prints the power and the expected final total sample size of a stratified design with
#' blinded sample size re-estimation returned by \code{\link{StratPowerBSSR}}.
#'
#' An object whose attributes were removed, for example by selecting some of its
#' columns, is printed as a data frame.
#'
#' @param x An object of class \code{bbssr_spowerbssr}
#' @param digits Number of decimal places. Default is 4
#' @param ... Further arguments, ignored
#'
#' @return \code{x}, invisibly
#'
#' @examples
#' res <- StratPowerBSSR(p = c(0.3, 0.5), Delta.A = 0.2, Delta.T = 0.2,
#'                       N1 = c(30, 30), N2 = c(30, 30), omega = 0.5, r = 1,
#'                       alpha = 0.025, tar.power = 0.8, Test = 'CMH',
#'                       ss.method = 'standard')
#' print(res)
#'
#' @export
print.bbssr_spowerbssr <- function(x, digits = 4, ...) {
  if (is.null(attr(x, 'ss.method')) ||
      !all(c('power.BSSR', 'power.TRAD', 'E.N') %in% names(x))) {
    return(NextMethod())
  }
  rule <- c(exact = 'exact sample size', standard = 'normal approximation')[[
    attr(x, 'ss.method')]]
  sizes <- function(n1, n2) {
    sprintf('(%s) and (%s)', paste(n1, collapse = ', '), paste(n2, collapse = ', '))
  }
  cat('Blinded sample size re-estimation for a stratified design\n\n')
  cat(sprintf('  Test             : %s\n', attr(x, 'Test')))
  cat(sprintf('  Alternative      : %s\n', attr(x, 'alternative')))
  cat(sprintf('  Re-estimation    : %s (%s%s)\n', rule,
              if (isTRUE(attr(x, 'restricted'))) 'restricted' else 'unrestricted',
              if (identical(attr(x, 'ss.strata'), 'pooled')) ', strata pooled' else ''))
  cat(sprintf('  Planned sizes    : %s\n', sizes(attr(x, 'N1'), attr(x, 'N2'))))
  alloc <- attr(x, 'alloc')
  if (!is.null(alloc) && alloc != 'fixed') {
    cat(sprintf('  Allocation       : %s\n', strat_alloc_label(alloc, attr(x, 'block'))))
    cat(sprintf('  Interim analysis : first %d patients%s\n', attr(x, 'n.interim'),
                if (is.null(attr(x, 'nseq'))) ''
                else sprintf(', %d allocation sequences', attr(x, 'nseq'))))
  } else {
    cat(sprintf('  Interim sizes    : %s\n', sizes(attr(x, 'n1.interim'),
                                                 attr(x, 'n2.interim'))))
  }
  cat(sprintf('  Effect           : %s, assumed %s, true %s\n', attr(x, 'effect'),
              format(attr(x, 'Delta.A')),
              paste(format(attr(x, 'Delta.T')), collapse = ', ')))
  cat(sprintf('  Alpha            : %s, target power %s\n', format(attr(x, 'alpha')),
              format(attr(x, 'tar.power'))))
  if (!is.null(attr(x, 'nsim'))) {
    cat(sprintf('  Simulation       : %d trials per scenario\n', attr(x, 'nsim')))
  }
  cat('\n')
  tab <- data.frame(lapply(unclass(x), identity), check.names = FALSE)
  for (v in c('power.BSSR', 'power.TRAD', 'se.BSSR', 'se.TRAD')) {
    if (!is.null(tab[[v]])) tab[[v]] <- round(tab[[v]], digits)
  }
  tab$E.N <- round(tab$E.N, 2)
  print(tab, row.names = FALSE)
  invisible(x)
}
