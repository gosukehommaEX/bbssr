#' Print Method for bbssr_scrp Objects
#'
#' Prints the design settings, the number of outcomes at which the conditional rejection
#' probabilities of a stratified design exceed the nominal level, their largest values
#' and, when the type I error rate has been decomposed, its largest value over the
#' supplied response probabilities, as returned by \code{\link{StratCondRejectBSSR}}.
#'
#' An object whose attributes were removed, for example by selecting some of its
#' columns, is printed as a data frame.
#'
#' @param x An object of class \code{bbssr_scrp}
#' @param digits Number of significant digits of the probabilities. Default is 4
#' @param ... Further arguments, ignored
#'
#' @return \code{x}, invisibly
#'
#' @examples
#' crp <- StratCondRejectBSSR(Delta.A = 0.3, N1 = c(8, 8), N2 = c(8, 8), omega = 0.5,
#'                            r = 1, alpha = 0.025, tar.power = 0.8, Test = 'CMH',
#'                            ss.method = 'standard')
#' print(crp)
#'
#' @export
print.bbssr_scrp <- function(x, digits = 4, ...) {
  m <- attr(x, 'max')
  if (is.null(m) || !all(c('s1', 's2', 'CRP.max', 'CRP.total.max') %in% names(x))) {
    return(NextMethod())
  }
  ex <- attr(x, 'exceed')
  sizes <- function(n1, n2) {
    sprintf('(%s) and (%s)', paste(n1, collapse = ', '), paste(n2, collapse = ', '))
  }
  cat('Conditional rejection probabilities of blinded sample size re-estimation\n',
      'for a stratified design, under equal response probabilities in every stratum\n\n',
      sep = '')
  cat(sprintf('  Test             : %s\n', attr(x, 'Test')))
  cat(sprintf('  Alternative      : %s\n', attr(x, 'alternative')))
  cat(sprintf('  Planned sizes    : %s\n', sizes(attr(x, 'N1'), attr(x, 'N2'))))
  cat(sprintf('  Interim sizes    : %s\n', sizes(attr(x, 'n1.interim'),
                                               attr(x, 'n2.interim'))))
  cat(sprintf('  Assumed effect   : %s (%s)\n', format(attr(x, 'Delta.A')),
              attr(x, 'effect')))
  cat(sprintf('  Nominal level    : %s\n', format(attr(x, 'alpha'))))
  cat(sprintf('  Interim outcomes : %d pairs (s1, s2)\n', nrow(x)))
  cat(sprintf('  Above the level  : %d outcomes for CRP, %d for CRP.total\n\n',
              ex[['CRP']], ex[['CRP.total']]))
  cat('Largest conditional rejection probability\n')
  tab <- m
  tab$value <- signif(tab$value, digits)
  print(tab, row.names = FALSE)
  tie <- attr(x, 'TIE')
  if (!is.null(tie)) {
    i <- which.max(tie$TIE)
    cat(sprintf(paste0('\nType I error rate over %d point%s: largest %s at ',
                       'theta = (%s, %s)\n'),
                nrow(tie), if (nrow(tie) == 1) '' else 's',
                format(signif(tie$TIE[i], digits)), format(signif(tie$theta.S1[i], 6)),
                format(signif(tie$theta.S2[i], 6))))
  }
  invisible(x)
}
