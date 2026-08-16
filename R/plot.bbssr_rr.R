#' Plot a Rejection Region
#'
#' Displays the outcome grid of a two-arm trial with a binary endpoint and shades the
#' outcomes for which the null hypothesis is rejected.
#'
#' @param x An object of class \code{bbssr_rr} returned by \code{\link{BinaryRR}}
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
#' Both axes count responders, so the breaks are restricted to integers. The vertical axis
#' runs downwards, which places the outcome with no responders in the top left corner and
#' matches the layout of the map printed by \code{\link{print.bbssr_rr}}.
#'
#' @examples
#' RR <- BinaryRR(N1 = 10, N2 = 10, alpha = 0.025, Test = 'Chisq')
#' plot(RR)
#'
#' # Enlarge the type and relabel the legend
#' plot(RR, base_size = 14, legend.title = 'Decision',
#'      legend.labels = c('do not reject', 'reject'),
#'      colours = c('white', 'grey30'))
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual labs
#' @importFrom ggplot2 scale_x_continuous scale_y_reverse
#' @importFrom ggplot2 theme_bw theme element_text coord_equal
plot.bbssr_rr <- function(x, main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
                          legend.title = NULL, legend.labels = NULL, colours = NULL,
                          base_size = 11, ...) {
  N1 <- attr(x, 'N1')
  N2 <- attr(x, 'N2')
  reject.levels <- c('retain', 'reject')
  values <- if (is.null(colours)) c('grey95', 'steelblue') else rep_len(colours, 2)
  names(values) <- reject.levels
  df <- expand.grid(x1 = 0:N1, x2 = 0:N2)
  df$Reject <- factor(
    ifelse(as.vector(x), 'reject', 'retain'),
    levels = reject.levels
  )
  ggplot(df, aes(x = x2, y = x1, fill = Reject)) +
    geom_tile(colour = 'grey85', linewidth = 0.2) +
    scale_x_continuous(breaks = integer_breaks(c(0, N2))) +
    scale_y_reverse(breaks = integer_breaks(c(0, N1))) +
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
        sprintf('N1 = %d, N2 = %d, alpha = %s', N1, N2, format(attr(x, 'alpha')))
      ),
      x = resolve_label(xlab, 'Number of responders in group 2'),
      y = resolve_label(ylab, 'Number of responders in group 1'),
      fill = resolve_label(legend.title, NULL)
    ) +
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(face = 'bold'))
}
