#' Plot the Exact Power Curve around a Sample Size Solution
#'
#' Recomputes the exact power over a range of sample sizes of group 2 and marks the
#' selected sample size and the target power.
#'
#' @param x An object of class \code{bbssr_samplesize} returned by
#'   \code{\link{BinarySampleSize}}
#' @param N2.range Optional integer vector of sample sizes of group 2 at which the power is
#'   evaluated. By default the selected sample size plus or minus five is used
#' @param main Title of the plot. The default of \code{NULL} builds the title from the
#'   object and a single \code{NA} drops it
#' @param sub Subtitle of the plot, following the convention of \code{main}
#' @param xlab Label of the horizontal axis, following the convention of \code{main}
#' @param ylab Label of the vertical axis, following the convention of \code{main}
#' @param ylim Numeric vector of length two giving the range of the vertical axis, or
#'   \code{NULL} for a range chosen from the data
#' @param ref.line Numeric vector of heights at which a dashed horizontal line is drawn.
#'   The default of \code{NULL} places one line at the target power and a single \code{NA}
#'   draws no line
#' @param ref.line.N2 Numeric vector of positions at which a dotted vertical line is drawn.
#'   The default of \code{NULL} places one line at the selected sample size and a single
#'   \code{NA} draws no line
#' @param show.points Logical. If \code{TRUE}, the default, a marker is drawn at every
#'   evaluated sample size. Set it to \code{FALSE} to draw the curve alone, which is easier
#'   to read over a wide range
#' @param colours Colour of the curve
#' @param base_size Base font size of the theme, in points
#' @param ... Further arguments, currently ignored
#'
#' @return A \code{ggplot} object
#'
#' @details
#' The power is recomputed at every point of \code{N2.range}, so a wide range combined with
#' one of the unconditional tests can take a long time to evaluate.
#'
#' The vertical range is imposed with \code{\link[ggplot2]{coord_cartesian}}, so points
#' outside \code{ylim} are hidden rather than removed from the data.
#'
#' @examples
#' \donttest{
#' ss <- BinarySampleSize(p1 = 0.4, p2 = 0.2, r = 1, alpha = 0.025,
#'                        tar.power = 0.8, Test = 'Chisq')
#' plot(ss)
#'
#' # Enlarge the type, rename the axes and draw neither reference line
#' plot(ss, base_size = 14, ref.line = NA, ref.line.N2 = NA,
#'      xlab = 'Sample size of the control group', ylab = 'Exact power')
#' }
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline geom_vline
#' @importFrom ggplot2 labs theme_bw theme element_text scale_x_continuous
#' @importFrom ggplot2 coord_cartesian
plot.bbssr_samplesize <- function(x, N2.range = NULL, main = NULL, sub = NULL,
                                  xlab = NULL, ylab = NULL, ylim = NULL,
                                  ref.line = NULL, ref.line.N2 = NULL,
                                  show.points = TRUE, colours = NULL,
                                  base_size = 11, ...) {
  if (is.null(N2.range)) N2.range <- seq(max(1, x$N2[1] - 5), x$N2[1] + 5)
  N2.range <- sort(unique(as.integer(N2.range)))
  colour <- if (is.null(colours)) 'steelblue' else colours[1]
  line.y <- resolve_label(ref.line, x$tar.power[1])
  line.x <- resolve_label(ref.line.N2, x$N2[1])
  default.sub <- if (is.null(line.y)) {
    sprintf('Selected N2 = %d', x$N2[1])
  } else if (is.null(ref.line)) {
    sprintf('Selected N2 = %d, dashed line is the target power of %s',
            x$N2[1], format(x$tar.power[1]))
  } else {
    sprintf('Selected N2 = %d, dashed line at %s',
            x$N2[1], paste(format(line.y), collapse = ', '))
  }
  tsmethod <- attr(x, 'tsmethod')
  n.grid <- attr(x, 'n.grid')
  bb.gamma <- attr(x, 'bb.gamma')
  Power <- vapply(N2.range, function(n2) {
    BinaryPower(x$p1[1], x$p2[1], ceiling(x$r[1] * n2), n2, x$alpha[1], x$Test[1],
                x$alternative[1], tsmethod, n.grid, bb.gamma)$Power
  }, numeric(1))
  df <- data.frame(N2 = N2.range, Power = Power)
  out <- ggplot(df, aes(x = N2, y = Power)) +
    geom_line(colour = colour, linewidth = 0.8)
  if (isTRUE(show.points)) {
    out <- out + geom_point(size = 2, colour = colour)
  }
  if (!is.null(line.y)) {
    out <- out + geom_hline(yintercept = line.y, linetype = 'dashed', colour = 'grey40')
  }
  if (!is.null(line.x)) {
    out <- out + geom_vline(xintercept = line.x, linetype = 'dotted', colour = 'firebrick')
  }
  out +
    scale_x_continuous(breaks = integer_breaks(range(N2.range))) +
    coord_cartesian(ylim = ylim) +
    labs(
      title = resolve_label(
        main,
        sprintf('Exact power by sample size, %s test (%s)', x$Test[1], x$alternative[1])
      ),
      subtitle = resolve_label(sub, default.sub),
      x = resolve_label(xlab, 'Sample size of group 2'),
      y = resolve_label(ylab, 'Power')
    ) +
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(face = 'bold'))
}
