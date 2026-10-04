#' Summary Method for bbssr_powerbssr Objects
#'
#' Summarizes the distribution of the final total sample size of a design with blinded
#' sample size re-estimation for every true pooled response probability.
#'
#' @param object An object of class \code{bbssr_powerbssr}
#' @param probs Probabilities of the reported quantiles. Default is
#'   \code{c(0.25, 0.5, 0.75)}
#' @param ... Further arguments, ignored
#'
#' @return A data frame with one row per element of \code{p} containing \code{p1},
#'   \code{p2}, \code{p}, the expected final total sample size \code{E.N}, its standard
#'   deviation \code{SD.N}, the quantiles of the final total sample size in columns named
#'   after \code{probs} (for example \code{N.q50} for the median), and the probability
#'   \code{P.N.max} that the final total reaches its largest attainable value
#'
#' @details
#' The quantile for the probability \code{q} is the smallest final total sample size
#' whose cumulative probability is at least \code{q}.
#'
#' @examples
#' res <- BinaryPowerBSSR(
#'   p = c(0.3, 0.45),
#'   Delta.A = 0.3, Delta.T = 0.3,
#'   N1 = 10, N2 = 10, omega = 0.5, r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
#' )
#' summary(res)
#'
#' @export
summary.bbssr_powerbssr <- function(object, probs = c(0.25, 0.5, 0.75), ...) {
  d <- attr(object, 'N.dist')
  if (is.null(d)) stop('the object carries no distribution of the final sample size')
  if (length(probs) < 1 || anyNA(probs) || any(probs <= 0 | probs > 1)) {
    stop('probs must be values in (0, 1]')
  }
  map <- attr(object, 'reestimation')
  N.top <- max(map$N1 + map$N2)
  rows <- lapply(seq_len(nrow(object)), function(k) {
    dk <- d[d$scenario == k, , drop = FALSE]
    agg <- tapply(dk$prob, dk$N, sum)
    N <- as.numeric(names(agg))
    prob <- as.numeric(agg)
    E <- sum(N * prob)
    cdf <- cumsum(prob)
    q <- vapply(probs, function(a) N[which(cdf >= a - 1e-12)[1]], numeric(1))
    c(E.N = E, SD.N = sqrt(sum((N - E)^2 * prob)), q, P.N.max = sum(prob[N == N.top]))
  })
  tab <- as.data.frame(do.call(rbind, rows))
  names(tab)[2 + seq_along(probs)] <- paste0('N.q', formatC(100 * probs, format = 'fg'))
  cbind(data.frame(p1 = object$p1, p2 = object$p2, p = object$p), tab)
}
