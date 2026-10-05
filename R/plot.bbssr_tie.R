#' Plot Method for bbssr_tie Objects
#'
#' Draws the type I error rate of a design with blinded sample size re-estimation and of
#' the corresponding fixed-sample design against the common response probability, with a
#' reference line at the nominal level.
#'
#' @param x An object of class \code{bbssr_tie}
#' @param main Title. \code{NULL} keeps the default and \code{NA} removes it
#' @param sub Subtitle. \code{NULL} keeps the default and \code{NA} removes it
#' @param xlab Label of the horizontal axis. \code{NULL} keeps the default
#' @param ylab Label of the vertical axis. \code{NULL} keeps the default
#' @param ylim Range of the vertical axis, applied through \code{coord_cartesian()}
#' @param ref.line Height of the dashed reference line. \code{NULL} draws it at the
#'   nominal level and \code{NA} removes it
#' @param legend.title Title of the legend. \code{NULL} keeps the default
#' @param legend.labels Labels of the two designs. \code{NULL} keeps the defaults
#' @param colours Two colours for the BSSR and the fixed-sample design
#' @param base_size Base font size passed to \code{theme_bw()}
#' @param ... Further arguments, ignored
#'
#' @return A \code{ggplot} object
#'
#' @examples
#' tie <- BinaryTypeIErrorBSSR(
#'   Delta.A = 0.3, N1 = 20, N2 = 20, omega = 0.5, r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Chisq', ss.method = 'standard',
#'   theta = seq(0.05, 0.95, by = 0.05)
#' )
#' plot(tie)
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_hline scale_colour_manual
#' @importFrom ggplot2 coord_cartesian labs theme_bw theme element_text
plot.bbssr_tie <- function(x, main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
                           ylim = NULL, ref.line = NULL, legend.title = NULL,
                           legend.labels = NULL, colours = NULL, base_size = 11, ...) {
  alpha <- attr(x, 'alpha')
  line.y <- resolve_label(ref.line, alpha)
  design.levels <- c('BSSR', 'Fixed sample')
  values <- if (is.null(colours)) c('steelblue', 'firebrick') else rep_len(colours, 2)
  names(values) <- design.levels
  df <- data.frame(
    theta = rep(x$theta, 2),
    TIE = c(x$TIE.BSSR, x$TIE.TRAD),
    Design = factor(rep(design.levels, each = nrow(x)), levels = design.levels)
  )
  out <- ggplot(df, aes(x = theta, y = TIE, colour = Design)) +
    geom_line(linewidth = 0.8)
  if (!is.null(line.y)) {
    out <- out + geom_hline(yintercept = line.y, linetype = 'dashed', colour = 'grey40')
  }
  out +
    scale_colour_manual(values = values,
                        labels = resolve_label(legend.labels, design.levels)) +
    coord_cartesian(ylim = ylim) +
    labs(
      title = resolve_label(main, sprintf('Type I error rate, %s test', attr(x, 'Test'))),
      subtitle = resolve_label(sub, sprintf('Dashed line at the nominal level of %s',
                                            format(alpha))),
      x = resolve_label(xlab, if (isTRUE(attr(x, 'margin') != 0)) {
        'Pooled response probability on the null boundary'
      } else {
        'Common response probability'
      }),
      y = resolve_label(ylab, 'Type I error rate'),
      colour = resolve_label(legend.title, 'Design')
    ) +
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(face = 'bold'))
}
