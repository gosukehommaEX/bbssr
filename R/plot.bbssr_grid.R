#' Plot Method for bbssr_grid Objects
#'
#' Draws the power or the expected total sample size of the designs of
#' \code{\link{BinaryGridBSSR}} against the true pooled response probability, with one
#' colour per value of a chosen column and, optionally, one panel per value of another.
#'
#' @param x An object of class \code{bbssr_grid}
#' @param what Quantity to draw. \code{'power'} (default) draws the power of the BSSR
#'   designs as solid lines and that of the fixed-sample designs as dashed lines.
#'   \code{'E.N'} draws the expected total sample size of the BSSR designs
#' @param colour.by Name of the column of \code{x} whose values are told apart by colour.
#'   The default of \code{NULL} uses the index \code{design}
#' @param facet.by Name of a column of \code{x} whose values are drawn in separate panels,
#'   or \code{NULL} (default) for a single panel
#' @param main Title. \code{NULL} keeps the default and \code{NA} removes it
#' @param sub Subtitle. \code{NULL} keeps the default and \code{NA} removes it
#' @param xlab Label of the horizontal axis. \code{NULL} keeps the default
#' @param ylab Label of the vertical axis. \code{NULL} keeps the default
#' @param ref.line Height of the dotted reference line for \code{what = 'power'}.
#'   \code{NULL} draws it at the target power when all designs share one, and \code{NA}
#'   removes it
#' @param legend.title Title of the colour legend. \code{NULL} uses \code{colour.by}
#' @param base_size Base font size passed to \code{theme_bw()}
#' @param ... Further arguments, ignored
#'
#' @return A \code{ggplot} object
#'
#' @examples
#' design <- expand.grid(Test = c('Chisq', 'Fisher'), omega = c(0.3, 0.5),
#'                       stringsAsFactors = FALSE)
#' grid <- BinaryGridBSSR(design, p = c(0.3, 0.4, 0.5), Delta.A = 0.3, N1 = 30, N2 = 30,
#'                        r = 1, alpha = 0.025, tar.power = 0.8, ss.method = 'standard')
#' plot(grid, colour.by = 'Test', facet.by = 'omega')
#' plot(grid, what = 'E.N', colour.by = 'Test', facet.by = 'omega')
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline facet_wrap labs theme_bw
#' @importFrom ggplot2 theme element_text
plot.bbssr_grid <- function(x, what = c('power', 'E.N'), colour.by = NULL,
                            facet.by = NULL, main = NULL, sub = NULL, xlab = NULL,
                            ylab = NULL, ref.line = NULL, legend.title = NULL,
                            base_size = 11, ...) {
  what <- match.arg(what)
  if (is.null(colour.by)) colour.by <- 'design'
  ok <- function(v) is.character(v) && length(v) == 1 && v %in% names(x)
  if (!ok(colour.by)) stop('colour.by must name a column of x')
  if (!is.null(facet.by) && !ok(facet.by)) stop('facet.by must name a column of x')
  col <- factor(x[[colour.by]])
  n <- nrow(x)
  if (what == 'power') {
    design.levels <- c('BSSR', 'Fixed sample')
    lev <- factor(rep(design.levels, each = n), levels = design.levels)
    df <- data.frame(p = rep(x$p, 2), value = c(x$power.BSSR, x$power.TRAD), Design = lev,
                     col = rep(col, 2))
    df$grp <- interaction(rep(x$design, 2), df$Design)
  } else {
    df <- data.frame(p = x$p, value = x$E.N, col = col, grp = factor(x$design))
  }
  if (!is.null(facet.by)) {
    # Panels in the order of the levels of a factor, or of the sorted values otherwise
    fv <- x[[facet.by]]
    lev <- if (is.factor(fv)) levels(fv) else as.character(sort(unique(fv)))
    df$facet <- factor(paste(facet.by, '=', rep_len(as.character(fv), nrow(df))),
                       levels = paste(facet.by, '=', lev))
  }
  out <- if (what == 'power') {
    ggplot(df, aes(x = p, y = value, colour = col, linetype = Design, group = grp))
  } else {
    ggplot(df, aes(x = p, y = value, colour = col, group = grp))
  }
  # Points keep a design visible when it is evaluated at a single value of p
  out <- out + geom_line(linewidth = 0.7) + geom_point(size = 1.2)
  tar.power <- attr(x, 'tar.power')
  line.y <- resolve_label(ref.line, if (length(tar.power) == 1) tar.power else NULL)
  if (what == 'power' && !is.null(line.y)) {
    out <- out + geom_hline(yintercept = line.y, linetype = 'dotted', colour = 'grey40')
  }
  if (!is.null(facet.by)) out <- out + facet_wrap(~facet)
  out +
    labs(
      title = resolve_label(main, if (what == 'power') {
        'Power of the BSSR designs'
      } else {
        'Expected total sample size of the BSSR designs'
      }),
      subtitle = resolve_label(sub, if (what == 'power') {
        'Solid lines: BSSR design, dashed lines: fixed-sample design'
      } else {
        NULL
      }),
      x = resolve_label(xlab, 'True pooled response probability'),
      y = resolve_label(ylab, if (what == 'power') {
        'Power'
      } else {
        'Expected total sample size'
      }),
      colour = resolve_label(legend.title, colour.by),
      linetype = NULL
    ) +
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(face = 'bold'))
}
