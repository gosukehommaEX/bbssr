#' Power for Paired Binary Data
#'
#' Calculates the exact power of a test comparing two treatments on a binary endpoint in a
#' matched-pairs design, for given probabilities of the two kinds of discordant pairs.
#'
#' @param p10 Probability that only treatment 1 responds in a pair. A vector is allowed
#' @param p01 Probability that only treatment 2 responds in a pair, of the same length as
#'   \code{p10}. A single value of either is used with every value of the other
#' @inheritParams PairedRR
#'
#' @return An object of class \code{bbssr_ppower}, a data frame with one row for each
#'   element of \code{p10}, containing:
#' \describe{
#'   \item{p10, p01}{Probabilities of the two kinds of discordant pairs}
#'   \item{N}{Number of pairs}
#'   \item{alpha}{Level of significance}
#'   \item{Test}{Name of the statistical test}
#'   \item{alternative}{Direction of the alternative hypothesis}
#'   \item{Power}{Exact power}
#' }
#' The attribute \code{margin} holds the margin.
#'
#' @details
#' The numbers of the two kinds of discordant pairs and of concordant pairs follow the
#' multinomial distribution with \code{N} trials and probabilities \code{p10}, \code{p01}
#' and \code{1 - p10 - p01}, and the power is the probability of the rejection region of
#' \code{\link{PairedRR}} under this distribution, computed by complete enumeration. The
#' probabilities of the concordant pairs with both and with neither treatment responding
#' do not enter, since no test depends on how the concordant pairs split.
#'
#' @examples
#' # Average exact power of the asymptotic McNemar test with 25 pairs over the alternatives
#' # of Berger and Sidik (2003, Section 5.1)
#' alt <- do.call(rbind, lapply(seq(0.025, 0.475, by = 0.05), function(q) {
#'   data.frame(p01 = q, p10 = seq(q + 0.05, 1 - q, by = 0.05))
#' }))
#' pw <- PairedPower(alt$p10, alt$p01, N = 25, alpha = 0.05, Test = 'McNemar')
#' mean(pw$Power)
#'
#' @export
PairedPower <- function(p10, p01, N, alpha, Test,
                        alternative = c('greater', 'less', 'two.sided', 'equivalence'),
                        n.grid = 100, bb.gamma = 0, margin = 0, ref.pvalue = FALSE) {
  alternative <- match.arg(alternative)
  n <- max(length(p10), length(p01))
  if (!(length(p10) %in% c(1, n)) || !(length(p01) %in% c(1, n))) {
    stop('p10 and p01 must have the same length, or one of them a single value')
  }
  if (is.numeric(p10)) p10 <- rep_len(p10, n)
  if (is.numeric(p01)) p01 <- rep_len(p01, n)
  check_paired_probs(p10, p01)
  a <- check_paired_args(N, alpha, Test, n.grid, bb.gamma, ref.pvalue, alternative, margin)
  # The Berger-Boos parameter matters only for the unconditional test
  if (a$Test != 'McNemar-uncond') bb.gamma <- 0
  rr <- get_paired_rr(a$N, alpha, a$Test, alternative, a$n.grid, bb.gamma, a$ref.pvalue,
                      a$margin)
  Power <- paired_power(rr, pmin(1, as.double(p10)), pmin(1, as.double(p01)))
  out <- data.frame(p10 = p10, p01 = p01, N = a$N, alpha = alpha, Test = a$Test,
                    alternative = alternative, Power = Power, stringsAsFactors = FALSE)
  attr(out, 'margin') <- a$margin
  class(out) <- c('bbssr_ppower', 'data.frame')
  out
}
