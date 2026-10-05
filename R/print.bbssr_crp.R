#' Print Method for bbssr_crp Objects
#'
#' Prints the design settings, the number of outcomes at which the conditional rejection
#' probabilities exceed the nominal level, their largest values and, when the type I error
#' rate has been decomposed, its largest value over the supplied response probabilities.
#'
#' @param x An object of class \code{bbssr_crp}
#' @param digits Number of significant digits of the probabilities. Default is 4
#' @param ... Further arguments, ignored
#'
#' @return \code{x}, invisibly
#'
#' @examples
#' crp <- BinaryCondRejectBSSR(
#'   Delta.A = 0.3, N1 = 12, N2 = 12, n.interim = c(6, 6), r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Fisher', ss.method = 'standard',
#'   theta = c(0.3, 0.5)
#' )
#' print(crp)
#'
#' @export
print.bbssr_crp <- function(x, digits = 4, ...) {
  m <- attr(x, 'max')
  if (is.null(m) || !all(c('s', 's2', 'CRP', 'CRP.total') %in% names(x))) {
    return(NextMethod())
  }
  ex <- attr(x, 'exceed')
  cat('Conditional rejection probabilities of blinded sample size re-estimation\n',
      'for a binary endpoint, under equal response probabilities\n\n', sep = '')
  cat(sprintf('  Test            : %s\n', attr(x, 'Test')))
  cat(sprintf('  Alternative     : %s\n', attr(x, 'alternative')))
  cat(sprintf('  Initial size    : N1 = %s, N2 = %s\n',
              format(attr(x, 'N1')), format(attr(x, 'N2'))))
  cat(sprintf('  Interim size    : n1 = %d, n2 = %d\n',
              attr(x, 'n1.interim'), attr(x, 'n2.interim')))
  cat(sprintf('  Assumed effect  : %s (%s)\n', format(attr(x, 'Delta.A')),
              attr(x, 'effect')))
  cat(sprintf('  Nominal level   : %s\n', format(attr(x, 'alpha'))))
  cat(sprintf('  Outcomes        : %d pairs (s, s2)\n', nrow(x)))
  cat(sprintf('  Above the level : %d for CRP, %d for CRP.total\n\n',
              ex[['CRP']], ex[['CRP.total']]))
  cat('Largest conditional rejection probability\n')
  tab <- m
  tab$value <- signif(tab$value, digits)
  print(tab, row.names = FALSE)
  tie <- attr(x, 'TIE')
  if (!is.null(tie)) {
    i <- which.max(tie$TIE)
    cat(sprintf(paste0('\nType I error rate over %d value%s of theta: ',
                       'largest %s at theta = %s\n'),
                nrow(tie), if (nrow(tie) == 1) '' else 's',
                format(signif(tie$TIE[i], digits)), format(signif(tie$theta[i], 6))))
  }
  invisible(x)
}
