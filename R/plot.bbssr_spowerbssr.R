#' Plot the Power of a Stratified BSSR Design against the Fixed Design
#'
#' Draws the power of the design with blinded sample size re-estimation of
#' \code{\link{StratPowerBSSR}} and that of the fixed design against a chosen column of
#' the scenarios, as \code{\link{plot.bbssr_powerbssr}} does for two groups.
#'
#' @param x An object of class \code{bbssr_spowerbssr} returned by
#'   \code{\link{StratPowerBSSR}}
#' @param x.by Name of the numeric column of \code{x} on the horizontal axis, or
#'   \code{'scenario'} for the rows of \code{p} that the scenarios come from. The default
#'   of \code{NULL} takes the first of the true pooled response probabilities
#'   \code{p.S1, p.S2, ...} that varies between the scenarios, or the scenario when none
#'   varies
#' @param xlab Label of the horizontal axis. \code{NULL} makes it from \code{x.by} and a
#'   single \code{NA} drops it
#' @param ... Further arguments of \code{\link{plot.bbssr_powerbssr}}: \code{main},
#'   \code{sub}, \code{ylab}, \code{ylim}, \code{ref.line}, \code{legend.title},
#'   \code{legend.labels}, \code{show.points}, \code{colours} and \code{base_size}
#'
#' @return A \code{ggplot} object
#'
#' @details
#' The lines join the scenarios in the order of the horizontal axis, so the scenarios
#' should differ in a single probability, or \code{x.by} should name a column in which
#' they all differ; a warning is given when two scenarios share a value on the axis.
#'
#' The Monte Carlo standard errors of a simulated or averaged power (three or more strata,
#' or a random allocation) are not drawn. With \code{Delta.T = 0} in
#' \code{\link{StratPowerBSSR}} the curves are rejection probabilities under the null
#' hypothesis, and \code{ref.line} can be set to the level of significance.
#'
#' @examples
#' res <- StratPowerBSSR(p = cbind(0.4, seq(0.3, 0.7, by = 0.1)), Delta.A = 0.3,
#'                       Delta.T = 0.3, N1 = c(10, 10), N2 = c(10, 10), omega = 0.5,
#'                       r = 1, alpha = 0.05, tar.power = 0.8, Test = 'CMH',
#'                       ss.method = 'standard')
#' plot(res)
#'
#' @export
plot.bbssr_spowerbssr <- function(x, x.by = NULL, xlab = NULL, ...) {
  ax <- strat_x_axis(x, x.by, xlab, attr(x, 'scenario'))
  y <- x
  y$p <- ax$x
  class(y) <- c('bbssr_powerbssr', 'data.frame')
  plot.bbssr_powerbssr(y, xlab = ax$xlab, ...)
}
