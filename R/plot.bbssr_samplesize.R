#' Plot the Exact Power Curve around a Sample Size Solution
#'
#' Recomputes the exact power over a range of sample sizes of group 2 and marks the
#' selected sample size and the target power.
#'
#' @param x An object of class \code{bbssr_samplesize} returned by
#'   \code{\link{BinarySampleSize}}
#' @param N2.range Optional integer vector of sample sizes of group 2 at which the power is
#'   evaluated. By default the selected sample size plus or minus five is used
#' @param ... Further arguments, currently ignored
#'
#' @return A \code{ggplot} object
#'
#' @details
#' The power is recomputed at every point of \code{N2.range}, so a wide range combined with
#' one of the unconditional tests can take a long time to evaluate.
#'
#' @examples
#' \donttest{
#' ss <- BinarySampleSize(p1 = 0.4, p2 = 0.2, r = 1, alpha = 0.025,
#'                        tar.power = 0.8, Test = 'Chisq')
#' plot(ss)
#' }
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline geom_vline
#' @importFrom ggplot2 labs theme_bw theme element_text scale_x_continuous
plot.bbssr_samplesize <- function(x, N2.range = NULL, ...) {
  if (is.null(N2.range)) N2.range <- seq(max(1, x$N2[1] - 5), x$N2[1] + 5)
  N2.range <- sort(unique(as.integer(N2.range)))
  tsmethod <- attr(x, 'tsmethod')
  n.grid <- attr(x, 'n.grid')
  bb.gamma <- attr(x, 'bb.gamma')
  Power <- vapply(N2.range, function(n2) {
    BinaryPower(x$p1[1], x$p2[1], ceiling(x$r[1] * n2), n2, x$alpha[1], x$Test[1],
                x$alternative[1], tsmethod, n.grid, bb.gamma)$Power
  }, numeric(1))
  df <- data.frame(N2 = N2.range, Power = Power)
  ggplot(df, aes(x = N2, y = Power)) +
    geom_line(colour = 'steelblue', linewidth = 0.8) +
    geom_point(size = 2, colour = 'steelblue') +
    geom_hline(yintercept = x$tar.power[1], linetype = 'dashed', colour = 'grey40') +
    geom_vline(xintercept = x$N2[1], linetype = 'dotted', colour = 'firebrick') +
    scale_x_continuous(breaks = integer_breaks(range(N2.range))) +
    labs(
      title = sprintf('Exact power by sample size, %s test (%s)', x$Test[1], x$alternative[1]),
      subtitle = sprintf('Selected N2 = %d, dashed line is the target power of %s',
                         x$N2[1], format(x$tar.power[1])),
      x = 'Sample size of group 2',
      y = 'Power'
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = 'bold'))
}
