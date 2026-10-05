#' Plot Method for bbssr_crp Objects
#'
#' Draws the conditional rejection probabilities of a design with blinded sample size
#' re-estimation over the pooled responder counts of the two stages, or the type I error
#' rate given the pooled number of interim responders.
#'
#' @param x An object of class \code{bbssr_crp}
#' @param what Quantity to draw. \code{'CRP'} (default) and \code{'CRP.total'} draw the
#'   corresponding column over the pairs \code{(s, s2)}, shaded from white to the first
#'   colour, and outline the outcomes at which it exceeds the nominal level in the second
#'   colour. \code{'by.s'} draws the type I error rate given \code{s} for each value of
#'   \code{theta}, which requires an object created with \code{theta}
#' @param main Title. \code{NULL} keeps the default and \code{NA} removes it
#' @param sub Subtitle. \code{NULL} keeps the default and \code{NA} removes it
#' @param xlab Label of the horizontal axis. \code{NULL} keeps the default
#' @param ylab Label of the vertical axis. \code{NULL} keeps the default
#' @param legend.title Title of the legend. \code{NULL} keeps the default
#' @param colours Two colours, for the shading of the probabilities and for the outline of
#'   the outcomes above the nominal level. They are not used for \code{what = 'by.s'}
#' @param base_size Base font size passed to \code{theme_bw()}
#' @param ... Further arguments, ignored
#'
#' @return A \code{ggplot} object
#'
#' @examples
#' crp <- BinaryCondRejectBSSR(
#'   Delta.A = 0.3, N1 = 12, N2 = 12, n.interim = c(6, 6), r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Fisher', ss.method = 'standard',
#'   theta = c(0.3, 0.5)
#' )
#' plot(crp)
#' plot(crp, what = 'by.s')
#'
#' @export
#' @import fpCompare
#' @importFrom ggplot2 ggplot aes geom_tile geom_line geom_point geom_hline labs
#' @importFrom ggplot2 scale_fill_gradient scale_x_continuous theme_bw theme element_text
plot.bbssr_crp <- function(x, what = c('CRP', 'CRP.total', 'by.s'), main = NULL,
                           sub = NULL, xlab = NULL, ylab = NULL, legend.title = NULL,
                           colours = NULL, base_size = 11, ...) {
  what <- match.arg(what)
  alpha <- attr(x, 'alpha')
  values <- if (is.null(colours)) c('steelblue', 'firebrick') else rep_len(colours, 2)
  n1 <- attr(x, 'n1.interim') + attr(x, 'n2.interim')
  if (what == 'by.s') {
    by.s <- attr(x, 'by.s')
    if (is.null(by.s)) stop("what = 'by.s' requires an object created with theta")
    lab <- format(by.s$theta)
    df <- data.frame(s = by.s$s, TIE.s = by.s$TIE.s,
                     theta = factor(lab, levels = unique(lab)))
    out <- ggplot(df, aes(x = s, y = TIE.s, colour = theta)) +
      geom_line(linewidth = 0.6) +
      geom_point(size = 1) +
      geom_hline(yintercept = alpha, linetype = 'dashed', colour = 'grey40') +
      scale_x_continuous(breaks = integer_breaks(c(0, n1))) +
      labs(
        title = resolve_label(main, sprintf('Type I error rate given s, %s test',
                                            attr(x, 'Test'))),
        subtitle = resolve_label(sub, sprintf('Dashed line at the nominal level of %s',
                                              format(alpha))),
        x = resolve_label(xlab, 'Pooled number of interim responders s'),
        y = resolve_label(ylab, 'Type I error rate given s'),
        colour = resolve_label(legend.title, 'theta')
      )
  } else {
    df <- data.frame(s = x$s, s2 = x$s2, value = x[[what]])
    out <- ggplot(df, aes(x = s, y = s2, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = 'white', high = values[1])
    above <- df[df$value %>>% alpha, , drop = FALSE]
    if (nrow(above) > 0) {
      out <- out + geom_tile(data = above, fill = NA, colour = values[2], linewidth = 0.5)
    }
    out <- out +
      scale_x_continuous(breaks = integer_breaks(c(0, n1))) +
      labs(
        title = resolve_label(main, sprintf('Conditional rejection probability, %s test',
                                            attr(x, 'Test'))),
        subtitle = resolve_label(sub, sprintf('%s, outlined where it exceeds %s', what,
                                              format(alpha))),
        x = resolve_label(xlab, 'Pooled number of interim responders s'),
        y = resolve_label(ylab, 'Pooled number of second-stage responders s2'),
        fill = resolve_label(legend.title, what)
      )
  }
  out +
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(face = 'bold'))
}
