#' Plot a Rejection Region
#'
#' Displays the outcome grid of a two-arm trial with a binary endpoint and shades the
#' outcomes for which the null hypothesis is rejected.
#'
#' @param x An object of class \code{bbssr_rr} returned by \code{\link{BinaryRR}}
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
#' @export
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual labs
#' @importFrom ggplot2 scale_x_continuous scale_y_reverse
#' @importFrom ggplot2 theme_bw theme element_text coord_equal
plot.bbssr_rr <- function(x, ...) {
  N1 <- attr(x, 'N1')
  N2 <- attr(x, 'N2')
  df <- expand.grid(x1 = 0:N1, x2 = 0:N2)
  df$Reject <- factor(
    ifelse(as.vector(x), 'reject', 'retain'),
    levels = c('retain', 'reject')
  )
  ggplot(df, aes(x = x2, y = x1, fill = Reject)) +
    geom_tile(colour = 'grey85', linewidth = 0.2) +
    scale_x_continuous(breaks = integer_breaks(c(0, N2))) +
    scale_y_reverse(breaks = integer_breaks(c(0, N1))) +
    scale_fill_manual(values = c(retain = 'grey95', reject = 'steelblue')) +
    coord_equal() +
    labs(
      title = sprintf('Rejection region, %s test (%s)', attr(x, 'Test'), attr(x, 'alternative')),
      subtitle = sprintf('N1 = %d, N2 = %d, alpha = %s', N1, N2, format(attr(x, 'alpha'))),
      x = 'Number of responders in group 2',
      y = 'Number of responders in group 1',
      fill = NULL
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = 'bold'))
}
