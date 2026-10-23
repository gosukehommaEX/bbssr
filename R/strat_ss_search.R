#' Exact Search for the Total Sample Size of a Stratified Design
#'
#' Internal helper returning, for each row of \code{p1} and \code{p2}, the nominal total
#' sample size at which the exact power of the stratified test attains the target power,
#' within the family of designs of \code{strat_sizes}. The searches are those of
#' \code{ss_exact_search}: \code{'crossing'} steps from the normal approximation of
#' \code{strat_raw_n}, rounded up, to a size that attains the target power while the size
#' one smaller does not; \code{'smallest'} scans upwards from 1; and \code{'stable'}
#' scans downwards from the limit given by \code{search.limit}, so that every size from the
#' result up to the limit attains the target power. Sizes that give the same design have
#' the same power, so the result is the smallest nominal size of its design.
#'
#' @param p1 Response probabilities of group 1, a matrix with one row per scenario and one
#'   column per stratum
#' @param p2 Response probabilities of group 2, of the same form
#' @param frac Fractions of the strata, summing to 1
#' @param r Allocation ratios of the strata
#' @param search \code{'crossing'}, \code{'smallest'} or \code{'stable'}
#' @param search.limit Factor and increment of the limit of the search \code{'stable'}
#' @inheritParams StratRR
#' @param tar.power Target power
#'
#' @return A list with the integer vectors \code{n} of nominal total sizes and
#'   \code{limit} of limits of the search \code{'stable'} (otherwise \code{NA})
#'
#' @keywords internal
#' @noRd
strat_ss_search <- function(p1, p2, frac, r, alpha, tar.power, Test, alternative, tsmethod,
                            search = 'crossing', search.limit = c(2, 50)) {
  if (search == 'stable' &&
      (!is.numeric(search.limit) || length(search.limit) != 2 || anyNA(search.limit) ||
         any(!is.finite(search.limit)) || search.limit[1] < 1 || search.limit[2] < 0)) {
    stop('search.limit must be two numbers, a factor of at least 1 and a non-negative ',
         'increment')
  }
  J <- length(frac)
  p1 <- matrix(as.double(p1), ncol = J)
  p2 <- matrix(as.double(p2), ncol = J)
  a <- check_strat_args(rep(1, J), rep(1, J), alpha, Test, alternative, tsmethod)
  init <- strat_raw_n(p1, p2, frac, r, alpha, tar.power, a$alternative, a$cc)
  if (any(!is.finite(init))) {
    stop('the response probabilities show no effect in the direction of the alternative')
  }
  power_at <- function(n, k) {
    sz <- strat_sizes(n, frac, r)
    strat_fixed_power(sz$N1[1, ], sz$N2[1, ], p1[k, , drop = FALSE],
                      p2[k, , drop = FALSE], alpha, a$Test, a$alternative, a$tsmethod)
  }
  out <- integer(nrow(p1))
  limit <- rep(NA_integer_, nrow(p1))
  for (k in seq_len(nrow(p1))) {
    n <- max(1, ceil_tol(init[k]))
    # The upward steps stop at this size, which a test with an effect in the direction of
    # the alternative reaches only if its power does not tend to 1
    cap <- max(100 * n, n + 1e4)
    up_to_target <- function(n) {
      while (power_at(n, k) %<<% tar.power) {
        n <- n + 1
        if (n > cap) {
          stop('the exact power does not attain the target power up to a total ',
               'sample size of ', cap)
        }
      }
      n
    }
    if (search == 'smallest') {
      out[k] <- as.integer(up_to_target(1))
      next
    }
    if (search == 'stable') {
      L <- max(ceiling(search.limit[1] * n), ceiling(n + search.limit[2]))
      if (power_at(L, k) %<<% tar.power) {
        stop('the total sample size ', L, ', the limit of the stable search, does not ',
             'attain the target power; raise search.limit')
      }
      n <- L
      while (n > 1 && (power_at(n - 1, k) %>=% tar.power)) n <- n - 1
      out[k] <- as.integer(n)
      limit[k] <- as.integer(L)
      next
    }
    Power <- power_at(n, k)
    if (Power %>=% tar.power) {
      while ((Power %>=% tar.power) && (n > 1)) {
        n <- n - 1
        Power <- power_at(n, k)
      }
      if (Power %<<% tar.power) n <- n + 1
    } else {
      n <- up_to_target(n + 1)
    }
    out[k] <- as.integer(n)
  }
  list(n = out, limit = limit)
}
