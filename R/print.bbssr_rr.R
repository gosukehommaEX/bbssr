#' Print a Rejection Region
#'
#' @param x An object of class \code{bbssr_rr} returned by \code{\link{BinaryRR}}
#' @param show.map Logical. If \code{TRUE}, the outcome grid is printed as a map in which
#'   \code{X} marks rejection of the null hypothesis. The default prints the map only when
#'   the grid has at most 400 cells
#' @param ... Further arguments, currently ignored
#'
#' @return The object \code{x}, returned invisibly
#'
#' @examples
#' RR <- BinaryRR(N1 = 5, N2 = 5, alpha = 0.025, Test = 'Chisq')
#' print(RR)
#'
#' @export
print.bbssr_rr <- function(x, show.map = NULL, ...) {
  N1 <- attr(x, 'N1')
  N2 <- attr(x, 'N2')
  Test <- attr(x, 'Test')
  alternative <- attr(x, 'alternative')
  cat('Rejection region for a two-arm trial with a binary endpoint\n\n')
  cat(sprintf('  Test          : %s\n', Test))
  cat(sprintf('  Alternative   : %s\n', alternative))
  if (alternative == 'two.sided' && Test %in% c('Fisher', 'Fisher-midP', 'Boschloo')) {
    cat(sprintf('  Two-sided rule: %s\n', attr(x, 'tsmethod')))
  }
  cat(sprintf('  Sample sizes  : N1 = %d, N2 = %d\n', N1, N2))
  cat(sprintf('  Alpha         : %s\n', format(attr(x, 'alpha'))))
  if (Test %in% c('Z-pool', 'Boschloo')) {
    cat(sprintf('  Grid points   : %d\n', attr(x, 'n.grid')))
    bb <- attr(x, 'bb.gamma')
    cat(sprintf('  Berger-Boos   : %s\n', if (bb > 0) format(bb) else 'not used'))
  }
  m <- matrix(as.vector(x), nrow = N1 + 1L, ncol = N2 + 1L)
  cat(sprintf('\n  Rejected outcomes: %d of %d\n\n', sum(m), length(m)))
  if (is.null(show.map)) show.map <- (N1 + 1) * (N2 + 1) <= 400
  if (isTRUE(show.map)) {
    chr <- ifelse(m, 'X', '.')
    dimnames(chr) <- list(paste0('x1=', 0:N1), paste0('x2=', 0:N2))
    print(noquote(chr))
    cat('\n  X marks rejection of the null hypothesis\n')
  } else {
    cat('  Outcome grid suppressed, use print(x, show.map = TRUE) to display it\n')
  }
  invisible(x)
}
