#' Plot Method for bbssr_stie Objects
#'
#' Plots the type I error rate of a two-stratum design with blinded sample size
#' re-estimation and of the fixed design over the grid of common response probabilities
#' of the two strata, as two heat maps.
#'
#' @param x An object of class \code{bbssr_stie} with two strata, computed on a grid
#' @param main Title, or \code{NULL} for the default
#' @param xlab Label of the horizontal axis, or \code{NULL} for the default
#' @param ylab Label of the vertical axis, or \code{NULL} for the default
#' @param low Colour of the smallest rate. Default is \code{'white'}
#' @param high Colour of the largest rate. Default is \code{'firebrick'}
#' @param base_size Base font size. Default is 11
#' @param ... Further arguments, ignored
#'
#' @return A \code{ggplot} object
#'
#' @examples
#' tie <- StratTypeIErrorBSSR(Delta.A = 0.2, N1 = c(20, 20), N2 = c(20, 20), omega = 0.5,
#'                            r = 1, alpha = 0.025, tar.power = 0.8, Test = 'CMH',
#'                            ss.method = 'standard', theta = seq(0.1, 0.9, by = 0.1))
#' plot(tie)
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_tile facet_wrap scale_fill_gradient labs
#' @importFrom ggplot2 coord_equal theme_bw theme element_text
plot.bbssr_stie <- function(x, main = NULL, xlab = NULL, ylab = NULL, low = 'white',
                            high = 'firebrick', base_size = 11, ...) {
  if (!all(c('theta.S1', 'theta.S2') %in% names(x)) || 'theta.S3' %in% names(x)) {
    stop('the plot is available for designs with two strata')
  }
  design.levels <- c('BSSR', 'Fixed sample')
  df <- data.frame(theta1 = rep(x$theta.S1, 2), theta2 = rep(x$theta.S2, 2),
                   TIE = c(x$TIE.BSSR, x$TIE.TRAD),
                   Design = factor(rep(design.levels, each = nrow(x)),
                                   levels = design.levels))
  ggplot(df, aes(x = theta1, y = theta2, fill = TIE)) +
    geom_tile() +
    facet_wrap(~ Design) +
    scale_fill_gradient(low = low, high = high) +
    coord_equal() +
    labs(title = resolve_label(main, sprintf('Type I error rate, %s test',
                                             attr(x, 'Test'))),
         x = resolve_label(xlab, 'Response probability of stratum 1'),
         y = resolve_label(ylab, 'Response probability of stratum 2'),
         fill = 'Type I error') +
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(face = 'bold'))
}
