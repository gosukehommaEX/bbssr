#' Plot Exact Power against the Response Probability of Group 1
#'
#' @param x An object of class \code{bbssr_power} returned by \code{\link{BinaryPower}}
#' @param ... Further arguments, currently ignored
#'
#' @return A \code{ggplot} object
#'
#' @examples
#' pw <- BinaryPower(p1 = seq(0.3, 0.8, by = 0.1), p2 = rep(0.2, 6),
#'                   N1 = 20, N2 = 20, alpha = 0.025, Test = 'Chisq')
#' plot(pw)
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs
#' @importFrom ggplot2 theme_bw theme element_text ylim
plot.bbssr_power <- function(x, ...) {
  df <- data.frame(p1 = x$p1, Power = x$Power)
  ggplot(df, aes(x = p1, y = Power)) +
    geom_line(colour = 'steelblue', linewidth = 0.8) +
    geom_point(size = 2, colour = 'steelblue') +
    ylim(0, 1) +
    labs(
      title = sprintf('Exact power, %s test (%s)', x$Test[1], x$alternative[1]),
      subtitle = sprintf('N1 = %d, N2 = %d, alpha = %s', x$N1[1], x$N2[1], format(x$alpha[1])),
      x = 'Response probability of group 1',
      y = 'Power'
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = 'bold'))
}
