#' Plot Exact Power against the Response Probability of Group 1
#'
#' @param x An object of class \code{bbssr_power} returned by \code{\link{BinaryPower}}
#' @param main Title of the plot. The default of \code{NULL} builds the title from the
#'   object and a single \code{NA} drops it
#' @param sub Subtitle of the plot, following the convention of \code{main}
#' @param xlab Label of the horizontal axis, following the convention of \code{main}
#' @param ylab Label of the vertical axis, following the convention of \code{main}
#' @param ylim Numeric vector of length two giving the range of the vertical axis, or
#'   \code{NULL} for a range chosen from the data
#' @param show.points Logical. If \code{TRUE}, the default, a marker is drawn at every
#'   evaluated response probability. Set it to \code{FALSE} to draw the curve alone, which
#'   is easier to read on a fine grid
#' @param colours Colour of the curve
#' @param base_size Base font size of the theme, in points
#' @param ... Further arguments, currently ignored
#'
#' @return A \code{ggplot} object
#'
#' @details
#' The vertical range is imposed with \code{\link[ggplot2]{coord_cartesian}}, so points
#' outside \code{ylim} are hidden rather than removed from the data.
#'
#' @examples
#' pw <- BinaryPower(p1 = seq(0.3, 0.8, by = 0.1), p2 = rep(0.2, 6),
#'                   N1 = 20, N2 = 20, alpha = 0.025, Test = 'Chisq')
#' plot(pw)
#'
#' # Enlarge the type, rescale the vertical axis, rename the axes and drop the markers
#' plot(pw, base_size = 14, ylim = c(0.2, 1), sub = NA, show.points = FALSE,
#'      xlab = 'Response probability, experimental group', ylab = 'Exact power')
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs coord_cartesian
#' @importFrom ggplot2 theme_bw theme element_text
plot.bbssr_power <- function(x, main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
                             ylim = c(0, 1), show.points = TRUE, colours = NULL,
                             base_size = 11, ...) {
  colour <- if (is.null(colours)) 'steelblue' else colours[1]
  df <- data.frame(p1 = x$p1, Power = x$Power)
  out <- ggplot(df, aes(x = p1, y = Power)) +
    geom_line(colour = colour, linewidth = 0.8)
  if (isTRUE(show.points)) {
    out <- out + geom_point(size = 2, colour = colour)
  }
  out +
    coord_cartesian(ylim = ylim) +
    labs(
      title = resolve_label(
        main,
        sprintf('Exact power, %s test (%s)', x$Test[1], x$alternative[1])
      ),
      subtitle = resolve_label(
        sub,
        sprintf('N1 = %d, N2 = %d, alpha = %s', x$N1[1], x$N2[1], format(x$alpha[1]))
      ),
      x = resolve_label(xlab, 'Response probability of group 1'),
      y = resolve_label(ylab, 'Power')
    ) +
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(face = 'bold'))
}
