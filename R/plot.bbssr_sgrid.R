#' Plot Method for bbssr_sgrid Objects
#'
#' Draws the power or the expected total sample size of the stratified designs of
#' \code{\link{StratGridBSSR}} against a chosen column, with one colour per value of
#' another column and, optionally, one panel per value of a third, as
#' \code{\link{plot.bbssr_grid}} does for two groups.
#'
#' @param x An object of class \code{bbssr_sgrid}
#' @param x.by Name of the numeric column of \code{x} on the horizontal axis. The default
#'   of \code{NULL} takes the first of the true pooled response probabilities
#'   \code{p.S1, p.S2, ...} that varies between the scenarios, or the index
#'   \code{scenario} when none varies
#' @param xlab Label of the horizontal axis. \code{NULL} makes it from \code{x.by}
#' @param ... Further arguments of \code{\link{plot.bbssr_grid}}: \code{what},
#'   \code{colour.by}, \code{facet.by}, \code{main}, \code{sub}, \code{ylab},
#'   \code{ref.line}, \code{legend.title} and \code{base_size}
#'
#' @return A \code{ggplot} object
#'
#' @examples
#' design <- data.frame(Test = c('CMH', 'Strat-Fisher'))
#' grid <- StratGridBSSR(design, p = rbind(c(0.3, 0.5), c(0.4, 0.6)), Delta.A = 0.3,
#'                       N1 = c(10, 10), N2 = c(10, 10), omega = 0.5, r = 1,
#'                       alpha = 0.05, tar.power = 0.8, ss.method = 'standard')
#' plot(grid, colour.by = 'Test')
#' plot(grid, what = 'E.N', x.by = 'p.S2', colour.by = 'Test')
#'
#' @export
plot.bbssr_sgrid <- function(x, x.by = NULL, xlab = NULL, ...) {
  ax <- strat_x_axis(x, x.by, xlab)
  y <- x
  y$p <- ax$x
  class(y) <- c('bbssr_grid', 'data.frame')
  plot.bbssr_grid(y, xlab = ax$xlab, ...)
}
