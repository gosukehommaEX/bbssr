#' Plot the Power of a BSSR Design against the Fixed-Sample Design
#'
#' @param x An object of class \code{bbssr_powerbssr} returned by
#'   \code{\link{BinaryPowerBSSR}}
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
#' @param legend.title Title of the legend. The default of \code{NULL} leaves the legend
#'   untitled
#' @param legend.labels Character vector of length two replacing the entries of the
#'   legend, in the order BSSR then fixed sample
#' @param show.points Logical. If \code{TRUE}, the default, a marker is drawn at every
#'   evaluated pooled probability. Set it to \code{FALSE} to draw the curves alone, which
#'   is easier to read on a fine grid
#' @param colours Character vector of length two giving the colour of the two curves, in
#'   the order BSSR then fixed sample
#' @param base_size Base font size of the theme, in points
#' @param ... Further arguments, currently ignored
#'
#' @return A \code{ggplot} object
#'
#' @details
#' Setting \code{Delta.T} to 0 in \code{\link{BinaryPowerBSSR}} turns the two power columns
#' into rejection probabilities under the null hypothesis. The default reference line at
#' the target power is then out of place, and the level of significance passed through
#' \code{ref.line} together with a rescaled \code{ylim} gives a readable plot of the type I
#' error rate.
#'
#' The vertical range is imposed with \code{\link[ggplot2]{coord_cartesian}}, so points
#' outside \code{ylim} are hidden rather than removed from the data.
#'
#' @examples
#' \donttest{
#' res <- BinaryPowerBSSR(
#'   p = seq(0.19, 0.37, by = 0.03),
#'   Delta.A = 0.36, Delta.T = 0.36,
#'   N1 = 24, N2 = 24, omega = 0.5, r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
#' )
#' plot(res)
#'
#' # Type I error rate, with the reference line moved to the level of significance
#' tie <- BinaryPowerBSSR(
#'   p = seq(0.19, 0.37, by = 0.03),
#'   Delta.A = 0.36, Delta.T = 0,
#'   N1 = 24, N2 = 24, omega = 0.5, r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
#' )
#' plot(tie, ref.line = 0.025, ylim = c(0.015, 0.030), base_size = 14,
#'      show.points = FALSE,
#'      main = 'Type I error rate of the BSSR design', ylab = 'Type I error rate',
#'      legend.title = 'Design', legend.labels = c('BSSR', 'Fixed'))
#' }
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline labs coord_cartesian
#' @importFrom ggplot2 theme_bw theme element_text scale_colour_manual
plot.bbssr_powerbssr <- function(x, main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
                                 ylim = NULL, ref.line = NULL, legend.title = NULL,
                                 legend.labels = NULL, show.points = TRUE,
                                 colours = NULL, base_size = 11, ...) {
  tar.power <- attr(x, 'tar.power')
  line.y <- resolve_label(ref.line, tar.power)
  rule <- if (isTRUE(attr(x, 'restricted'))) 'Restricted' else 'Unrestricted'
  default.sub <- if (is.null(line.y)) {
    sprintf('%s rule', rule)
  } else if (is.null(ref.line)) {
    sprintf('%s rule, dashed line is the target power of %s', rule, format(tar.power))
  } else {
    sprintf('%s rule, dashed line at %s', rule,
            paste(format(line.y), collapse = ', '))
  }
  design.levels <- c('BSSR', 'Fixed sample')
  values <- if (is.null(colours)) c('steelblue', 'firebrick') else rep_len(colours, 2)
  names(values) <- design.levels
  df <- data.frame(
    p = rep(x$p, 2),
    Power = c(x$power.BSSR, x$power.TRAD),
    Design = factor(rep(design.levels, each = nrow(x)), levels = design.levels)
  )
  out <- ggplot(df, aes(x = p, y = Power, colour = Design)) +
    geom_line(linewidth = 0.8)
  if (isTRUE(show.points)) {
    out <- out + geom_point(size = 2)
  }
  if (!is.null(line.y)) {
    out <- out + geom_hline(yintercept = line.y, linetype = 'dashed', colour = 'grey40')
  }
  out +
    scale_colour_manual(
      values = values,
      labels = resolve_label(legend.labels, design.levels)
    ) +
    coord_cartesian(ylim = ylim) +
    labs(
      title = resolve_label(
        main,
        sprintf('Power of the BSSR design, %s test (%s)',
                attr(x, 'Test'), attr(x, 'alternative'))
      ),
      subtitle = resolve_label(sub, default.sub),
      x = resolve_label(xlab, 'True pooled response probability'),
      y = resolve_label(ylab, 'Power'),
      colour = resolve_label(legend.title, NULL)
    ) +
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(face = 'bold'))
}
