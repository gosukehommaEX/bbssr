#' Exact Power of a Two-Stratum Design with Random Stratum and Group Sizes
#'
#' Internal helper returning the exact rejection probability of a two-stratum design whose
#' stratum and group sizes are random, for each row of \code{p1} and \code{p2}, as the sum
#' over the configurations of stratum and group sizes of their probability times the
#' rejection probability of the fixed design they give. With \code{alloc = 'group'} the
#' stratum sizes are \code{n1 = floor(N frac[1])} and \code{n2 = N - n1}, the rule of Jung
#' (2014, Section 3.2).
#'
#' @param N Total sample size
#' @param al List returned by \code{check_strat_alloc}
#' @param alloc \code{'stratum'}, \code{'group'} or \code{'both'}
#' @param p1 Response probabilities of group 1, a matrix with one row per scenario and two
#'   columns
#' @param p2 Response probabilities of group 2, of the same form
#' @param eps Configurations with a probability below \code{eps} are left out. Default is
#'   \code{1e-14}
#' @inheritParams StratRR
#'
#' @return A list with the vector \code{power}, the total probability \code{skipped} of
#'   the configurations left out, which bounds the error of every element of
#'   \code{power}, the number \code{configs} of configurations evaluated and the stratum
#'   sizes \code{n} of \code{alloc = 'group'} (otherwise \code{NULL})
#'
#' @keywords internal
#' @noRd
strat_alloc_power <- function(N, al, alloc, p1, p2, alpha, Test, alternative, tsmethod,
                              eps = 1e-14) {
  a <- check_strat_args(c(1, 1), c(1, 1), alpha, Test, alternative, tsmethod)
  n1 <- floor(N * al$frac[1] + 1e-9)
  nfix <- as.integer(c(n1, N - n1))
  p1 <- matrix(as.double(p1), ncol = 2)
  p2 <- matrix(as.double(p2), ncol = 2)
  res <- strat_random_power(match(alloc, c('stratum', 'group', 'both')), as.integer(N),
                            al$frac[1], nfix, al$b,
                            if (is.null(al$block)) 4L else al$block,
                            if (is.null(al$g)) c(1L, 1L) else al$g, p1, p2, a$type,
                            a$alt, a$ts, alpha, eps)
  res$n <- if (alloc == 'group') nfix else NULL
  res
}
