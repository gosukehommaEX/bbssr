#' Summary Method for bbssr_grid Objects
#'
#' Summarizes each design of \code{\link{BinaryGridBSSR}} over the true pooled response
#' probabilities at which it was evaluated.
#'
#' @param object An object of class \code{bbssr_grid}
#' @param ... Further arguments, ignored
#'
#' @return A data frame with one row per design, containing the attribute \code{designs}
#'   of \code{object} followed by the smallest and the largest power of the BSSR design
#'   (\code{power.BSSR.min}, \code{power.BSSR.max}) and of the fixed-sample design
#'   (\code{power.TRAD.min}, \code{power.TRAD.max}), and the largest expected total
#'   sample size \code{E.N.max}
#'
#' @examples
#' design <- data.frame(Test = c('Chisq', 'Fisher'))
#' grid <- BinaryGridBSSR(design, p = c(0.3, 0.4, 0.5), Delta.A = 0.3, N1 = 30, N2 = 30,
#'                        omega = 0.5, r = 1, alpha = 0.025, tar.power = 0.8,
#'                        ss.method = 'standard')
#' summary(grid)
#'
#' @export
summary.bbssr_grid <- function(object, ...) {
  d <- attr(object, 'designs')
  if (is.null(d)) stop('the object carries no table of designs')
  idx <- split(seq_len(nrow(object)), factor(object$design, levels = d$design))
  over <- function(column, fun) {
    vapply(idx, function(k) fun(object[[column]][k]), numeric(1), USE.NAMES = FALSE)
  }
  out <- cbind(d, data.frame(power.BSSR.min = over('power.BSSR', min),
                             power.BSSR.max = over('power.BSSR', max),
                             power.TRAD.min = over('power.TRAD', min),
                             power.TRAD.max = over('power.TRAD', max),
                             E.N.max = over('E.N', max)))
  rownames(out) <- NULL
  out
}
