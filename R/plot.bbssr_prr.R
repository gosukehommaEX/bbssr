#' Plot a Rejection Region for Paired Binary Data
#'
#' Displays the possible outcomes of a matched-pairs trial and shades the outcomes for
#' which the null hypothesis is rejected.
#'
#' @param x An object of class \code{bbssr_prr} returned by \code{\link{PairedRR}}
#' @param main Title of the plot. The default of \code{NULL} builds the title from the
#'   object and a single \code{NA} drops it
#' @param sub Subtitle of the plot, following the convention of \code{main}
#' @param xlab Label of the horizontal axis, following the convention of \code{main}
#' @param ylab Label of the vertical axis, following the convention of \code{main}
#' @param legend.title Title of the legend. The default of \code{NULL} leaves the legend
#'   untitled
#' @param legend.labels Character vector of length two replacing the entries of the
#'   legend, in the order retained then rejected
#' @param colours Character vector of length two giving the fill of the tiles, in the order
#'   retained then rejected
#' @param base_size Base font size of the theme, in points
#' @param ... Further arguments, currently ignored
#'
#' @return A \code{ggplot} object
#'
#' @details
#' The vertical axis counts the pairs in which only treatment 1 responds and runs
#' downwards, and the horizontal axis counts the pairs in which only treatment 2 responds,
#' as in the map printed by \code{\link{print.bbssr_prr}}. Only the possible outcomes, with
#' at most \code{N} discordant pairs, are drawn.
#'
#' @examples
#' RR <- PairedRR(N = 20, alpha = 0.025, Test = 'McNemar-midP')
#' plot(RR)
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual labs
#' @importFrom ggplot2 scale_x_continuous scale_y_reverse
#' @importFrom ggplot2 theme_bw theme element_text coord_equal
plot.bbssr_prr <- function(x, main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
                           legend.title = NULL, legend.labels = NULL, colours = NULL,
                           base_size = 11, ...) {
  N <- attr(x, 'N')
  reject.levels <- c('retain', 'reject')
  values <- if (is.null(colours)) c('grey95', 'steelblue') else rep_len(colours, 2)
  names(values) <- reject.levels
  df <- expand.grid(a = 0:N, b = 0:N)
  df$Reject <- factor(ifelse(as.vector(x), 'reject', 'retain'), levels = reject.levels)
  df <- df[df$a + df$b <= N, , drop = FALSE]
  ggplot(df, aes(x = b, y = a, fill = Reject)) +
    geom_tile(colour = 'grey85', linewidth = 0.2) +
    scale_x_continuous(breaks = integer_breaks(c(0, N))) +
    scale_y_reverse(breaks = integer_breaks(c(0, N))) +
    scale_fill_manual(
      values = values,
      labels = resolve_label(legend.labels, reject.levels)
    ) +
    coord_equal() +
    labs(
      title = resolve_label(
        main,
        sprintf('Rejection region, %s test (%s)', attr(x, 'Test'), attr(x, 'alternative'))
      ),
      subtitle = resolve_label(
        sub,
        sprintf('N = %d pairs, alpha = %s', N, format(attr(x, 'alpha')))
      ),
      x = resolve_label(xlab, 'Pairs in which only treatment 2 responds'),
      y = resolve_label(ylab, 'Pairs in which only treatment 1 responds'),
      fill = resolve_label(legend.title, NULL)
    ) +
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(face = 'bold'))
}
