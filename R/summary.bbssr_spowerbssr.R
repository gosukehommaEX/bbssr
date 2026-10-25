#' Summary Method for bbssr_spowerbssr Objects
#'
#' Summarizes the distribution of the final total sample size of a stratified design with
#' blinded sample size re-estimation in each scenario of \code{\link{StratPowerBSSR}}.
#'
#' @param object An object of class \code{bbssr_spowerbssr}
#' @param probs Probabilities of the reported quantiles. Default is
#'   \code{c(0.25, 0.5, 0.75)}
#' @param ... Further arguments, ignored
#'
#' @return A data frame with one row per scenario containing the true pooled response
#'   probabilities \code{p.S1, p.S2, ...} of the strata, those of the two groups
#'   \code{p1.S1, ...} and \code{p2.S1, ...}, the expected value \code{E.N} and the
#'   standard deviation \code{SD.N} of the final total sample size, its quantiles
#'   (\code{N.q25}, \code{N.q50} and \code{N.q75} by default, the smallest total whose
#'   distribution function reaches the probability), and the probability \code{P.N.max}
#'   that it is the largest total the design can reach. With two strata and the fixed
#'   allocation this is the largest final total of the re-estimation over all interim
#'   outcomes. Otherwise it is \code{N.max}, and \code{P.N.max} is \code{NA} when
#'   \code{N.max} was not given.
#'
#' @details
#' The distribution is the attribute \code{N.dist} of \code{object}: exact with two strata
#' and the fixed allocation, the average over the allocation sequences with two strata and
#' a random allocation, and that of the simulated trials with three or more strata.
#'
#' @examples
#' res <- StratPowerBSSR(p = c(0.4, 0.6), Delta.A = 0.3, Delta.T = 0.3, N1 = c(10, 10),
#'                       N2 = c(10, 10), omega = 0.5, r = 1, alpha = 0.05,
#'                       tar.power = 0.8, Test = 'CMH', ss.method = 'standard')
#' summary(res)
#'
#' @export
summary.bbssr_spowerbssr <- function(object, probs = c(0.25, 0.5, 0.75), ...) {
  d <- attr(object, 'N.dist')
  if (is.null(d)) stop('the object carries no distribution of the final sample size')
  if (!is.numeric(probs) || length(probs) < 1 || anyNA(probs) ||
      any(probs <= 0 | probs > 1)) {
    stop('probs must be values in (0, 1]')
  }
  map <- attr(object, 'reestimation')
  N.top <- if (!is.null(map)) {
    max(map$N)
  } else if (!is.null(attr(object, 'N.max'))) {
    attr(object, 'N.max')
  } else {
    NA_real_
  }
  rows <- lapply(seq_len(nrow(object)), function(k) {
    dk <- d[d$scenario == k, , drop = FALSE]
    agg <- tapply(dk$prob, dk$N, sum)
    N <- as.numeric(names(agg))
    prob <- as.numeric(agg)
    E <- sum(N * prob)
    cdf <- cumsum(prob)
    q <- vapply(probs, function(a) N[which(cdf >= a - 1e-12)[1]], numeric(1))
    top <- if (is.na(N.top)) NA_real_ else sum(prob[N == N.top])
    c(E.N = E, SD.N = sqrt(sum((N - E)^2 * prob)), q, P.N.max = top)
  })
  tab <- as.data.frame(do.call(rbind, rows))
  names(tab)[2 + seq_along(probs)] <- paste0('N.q', formatC(100 * probs, format = 'fg'))
  J <- sum(grepl('^p\\.S[0-9]+$', names(object)))
  cols <- c(paste0('p.S', seq_len(J)), paste0('p1.S', seq_len(J)),
            paste0('p2.S', seq_len(J)))
  out <- cbind(as.data.frame(lapply(object[cols], as.numeric)), tab)
  rownames(out) <- NULL
  out
}
