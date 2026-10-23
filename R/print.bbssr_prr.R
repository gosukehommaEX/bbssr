#' Print Method for bbssr_prr Objects
#'
#' Prints the settings of a rejection region for paired binary data and, for a small
#' number of pairs, a map of the outcomes.
#'
#' @param x An object of class \code{bbssr_prr} returned by \code{\link{PairedRR}}
#' @param show.map Logical. Whether the map of the outcomes is printed. The default of
#'   \code{NULL} prints it when there are at most 20 pairs
#' @param ... Further arguments, ignored
#'
#' @details
#' The map has one row for each number \code{a} of pairs in which only treatment 1
#' responds and one column for each number \code{b} of pairs in which only treatment 2
#' responds. \code{X} marks a rejected outcome, \code{.} a retained outcome, and the
#' impossible outcomes \code{a + b > N} are left blank.
#'
#' @return \code{x}, invisibly
#'
#' @examples
#' RR <- PairedRR(N = 15, alpha = 0.05, Test = 'McNemar')
#' print(RR)
#'
#' @export
print.bbssr_prr <- function(x, show.map = NULL, ...) {
  N <- attr(x, 'N')
  Test <- attr(x, 'Test')
  cat('Rejection region for paired binary data\n\n')
  cat(sprintf('  Test          : %s\n', Test))
  cat(sprintf('  Alternative   : %s\n', attr(x, 'alternative')))
  cat(sprintf('  Pairs         : N = %d\n', N))
  cat(sprintf('  Alpha         : %s\n', format(attr(x, 'alpha'))))
  margin <- format_margin(attr(x, 'margin'), 'RD')
  if (!is.null(margin)) {
    cat(sprintf('  Margin        : %s\n', margin))
  }
  if (Test == 'McNemar-uncond') {
    cat(sprintf('  Grid points   : %d\n', attr(x, 'n.grid')))
    bb <- attr(x, 'bb.gamma')
    cat(sprintf('  Berger-Boos   : %s\n', if (bb > 0) format(bb) else 'not used'))
    ref <- if (isTRUE(attr(x, 'ref.pvalue'))) 'local maxima refined' else 'not used'
    cat(sprintf('  Refinement    : %s\n', ref))
  }
  m <- matrix(as.vector(x), nrow = N + 1L, ncol = N + 1L)
  possible <- outer(0:N, 0:N, '+') <= N
  cat(sprintf('\n  Rejected outcomes: %d of %d\n\n', sum(m), sum(possible)))
  if (is.null(show.map)) show.map <- N <= 20
  if (isTRUE(show.map)) {
    chr <- ifelse(m, 'X', '.')
    chr[!possible] <- ' '
    dimnames(chr) <- list(paste0('a=', 0:N), paste0('b=', 0:N))
    print(noquote(chr))
    cat('\n  a: pairs in which only treatment 1 responds\n')
    cat('  b: pairs in which only treatment 2 responds\n')
    cat('  X marks rejection of the null hypothesis\n')
  } else {
    cat('  Outcome map suppressed, use print(x, show.map = TRUE) to display it\n')
  }
  invisible(x)
}
