#' Plot the Power of a BSSR Design against the Fixed-Sample Design
#'
#' @param x An object of class \code{bbssr_powerbssr} returned by
#'   \code{\link{BinaryPowerBSSR}}
#' @param ... Further arguments, currently ignored
#'
#' @return A \code{ggplot} object
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
#' }
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline labs
#' @importFrom ggplot2 theme_bw theme element_text scale_colour_manual
plot.bbssr_powerbssr <- function(x, ...) {
  df <- data.frame(
    p = rep(x$p, 2),
    Power = c(x$power.BSSR, x$power.TRAD),
    Design = factor(rep(c('BSSR', 'Fixed sample'), each = nrow(x)),
                    levels = c('BSSR', 'Fixed sample'))
  )
  ggplot(df, aes(x = p, y = Power, colour = Design)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_hline(yintercept = attr(x, 'tar.power'), linetype = 'dashed', colour = 'grey40') +
    scale_colour_manual(values = c('BSSR' = 'steelblue', 'Fixed sample' = 'firebrick')) +
    labs(
      title = sprintf('Power of the BSSR design, %s test (%s)',
                      attr(x, 'Test'), attr(x, 'alternative')),
      subtitle = sprintf('%s rule, dashed line is the target power of %s',
                         if (isTRUE(attr(x, 'restricted'))) 'Restricted' else 'Unrestricted',
                         format(attr(x, 'tar.power'))),
      x = 'True pooled response probability',
      y = 'Power',
      colour = NULL
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = 'bold'))
}
