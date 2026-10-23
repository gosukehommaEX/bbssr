#' Print Method for bbssr_psamplesize Objects
#'
#' Prints the number of pairs required for a test for paired binary data.
#'
#' @param x An object of class \code{bbssr_psamplesize} returned by
#'   \code{\link{PairedSampleSize}}
#' @param digits Number of decimal places of the power. Default is 4
#' @param ... Further arguments, ignored
#'
#' @return \code{x}, invisibly
#'
#' @examples
#' ss <- PairedSampleSize(p10 = 0.2, p01 = 0.1, alpha = 0.05, tar.power = 0.9,
#'                        Test = 'McNemar', method = 'standard')
#' print(ss)
#'
#' @export
print.bbssr_psamplesize <- function(x, digits = 4, ...) {
  if (!all(c('p10', 'p01', 'alpha', 'tar.power', 'Test', 'alternative', 'Power', 'N') %in%
           names(x))) return(NextMethod())
  cat('Number of pairs for paired binary data\n\n')
  cat(sprintf('  Test           : %s\n', x$Test[1]))
  cat(sprintf('  Alternative    : %s\n', x$alternative[1]))
  cat(sprintf('  Discordance    : p10 = %s, p01 = %s\n', format(x$p10[1]),
              format(x$p01[1])))
  margin <- format_margin(attr(x, 'margin'), 'RD')
  if (!is.null(margin)) {
    cat(sprintf('  Margin         : %s\n', margin))
  }
  cat(sprintf('  Alpha          : %s\n', format(x$alpha[1])))
  cat(sprintf('  Target power   : %s\n', format(x$tar.power[1])))
  method <- attr(x, 'method')
  if (identical(method, 'exact')) {
    search <- attr(x, 'search')
    limit <- attr(x, 'search.limit')
    cat(sprintf('  Exact search   : %s\n',
                if (identical(search, 'stable') && !is.null(limit) && !is.na(limit)) {
                  sprintf('stable up to N = %d', as.integer(limit))
                } else {
                  search
                }))
  } else if (identical(method, 'standard')) {
    cat('  Method         : normal approximation\n')
  }
  cat('\n')
  cat(sprintf('  Required number of pairs: N = %d\n', x$N[1]))
  cat(sprintf('  Attained power          : %s\n', format(round(x$Power[1], digits))))
  invisible(x)
}
