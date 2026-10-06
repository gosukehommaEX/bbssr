#' Print Method for bbssr_alphaadj Objects
#'
#' Prints the adjusted nominal significance levels of a design with blinded sample size
#' re-estimation and of the corresponding fixed-sample design.
#'
#' @param x An object of class \code{bbssr_alphaadj}
#' @param digits Number of significant digits. Default is 6
#' @param ... Further arguments, ignored
#'
#' @return \code{x}, invisibly
#'
#' @examples
#' \donttest{
#' adj <- BinaryAlphaAdjBSSR(
#'   Delta.A = 0.3, N1 = 39, N2 = 39, n.interim = c(20, 20), r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Chisq', ss.method = 'standard',
#'   theta = seq(0.01, 0.99, by = 0.01)
#' )
#' print(adj)
#' }
#'
#' @export
print.bbssr_alphaadj <- function(x, digits = 6, ...) {
  if (!all(c('Design', 'alpha', 'max.TIE', 'alpha.adj', 'max.TIE.adj') %in% names(x))) {
    return(NextMethod())
  }
  cat('Adjusted significance level for blinded sample size re-estimation\n\n')
  cat(sprintf('  Test            : %s\n', attr(x, 'Test')))
  cat(sprintf('  Alternative     : %s\n', attr(x, 'alternative')))
  cat(sprintf('  Adjusted part   : %s\n',
              if (attr(x, 'adjust') == 'test') 'final analysis only'
              else 'final analysis and re-estimation'))
  margin <- attr(x, 'margin')
  if (!is.null(margin) && margin != 0) {
    cat(sprintf('  Margin          : %s\n', format(margin)))
  }
  cat(sprintf('  Maximum         : %s\n', maximize_label(x)))
  cat(sprintf('  Target level    : %s\n\n', format(x$alpha[1])))
  tab <- data.frame(
    Design = ifelse(x$Design == 'BSSR', 'BSSR', 'Fixed sample'),
    max.TIE = signif(x$max.TIE, digits),
    alpha.adj = signif(x$alpha.adj, digits),
    max.TIE.adj = signif(x$max.TIE.adj, digits)
  )
  print(tab, row.names = FALSE)
  invisible(x)
}
