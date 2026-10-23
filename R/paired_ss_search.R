#' Exact Search for the Number of Pairs
#'
#' Internal helper returning, for each pair of probabilities \code{(p10[k], p01[k])}, the
#' number of pairs at which the exact power of the test of \code{\link{PairedRR}} attains
#' the target power, by the searches of \code{ss_exact_search}. The search
#' \code{'crossing'} steps from the normal approximation of \code{paired_raw_n}, rounded
#' up, to a number of pairs that attains the target power while the number one smaller
#' does not. The search \code{'smallest'} scans upwards from 1, and the search
#' \code{'stable'} scans downwards from the limit given by \code{search.limit}. The
#' rejection regions are computed once per search and shared by all pairs of
#' probabilities.
#'
#' @param p10 Probability that only treatment 1 responds, a vector
#' @param p01 Probability that only treatment 2 responds, of the same length
#' @inheritParams PairedSampleSize
#'
#' @return A list with the integer vector \code{N} of numbers of pairs and the integer
#'   vector \code{limit} of limits of the search \code{'stable'} (otherwise \code{NA})
#'
#' @keywords internal
#' @noRd
paired_ss_search <- function(p10, p01, alpha, tar.power, Test, alternative, n.grid,
                             bb.gamma, ref.pvalue, margin, search = 'crossing',
                             search.limit = c(2, 50)) {
  if (search == 'stable' &&
      (!is.numeric(search.limit) || length(search.limit) != 2 || anyNA(search.limit) ||
         any(!is.finite(search.limit)) || search.limit[1] < 1 || search.limit[2] < 0)) {
    stop('search.limit must be two numbers, a factor of at least 1 and a non-negative ',
         'increment')
  }
  store <- new.env(parent = emptyenv())
  # Rejection region at a given number of pairs, obtained once per search
  rr_at <- function(N) {
    key <- as.character(N)
    hit <- store[[key]]
    if (is.null(hit)) {
      hit <- get_paired_rr(N, alpha, Test, alternative, n.grid, bb.gamma, ref.pvalue,
                           margin)
      assign(key, hit, envir = store)
    }
    hit
  }
  power_at <- function(N, k) paired_power(rr_at(N), p10[k], p01[k])
  init <- paired_raw_n(p10, p01, alpha, tar.power, alternative, margin)
  out <- integer(length(p10))
  limit <- rep(NA_integer_, length(p10))
  for (k in seq_along(p10)) {
    N <- max(1, ceil_tol(init[k]))
    if (search == 'smallest') {
      N <- 1
      while (power_at(N, k) %<<% tar.power) N <- N + 1
      out[k] <- as.integer(N)
      next
    }
    if (search == 'stable') {
      L <- max(ceiling(search.limit[1] * N), ceiling(N + search.limit[2]))
      if (power_at(L, k) %<<% tar.power) {
        stop('the number of pairs ', L, ', the limit of the stable search, does not ',
             'attain the target power; raise search.limit')
      }
      N <- L
      while (N > 1 && (power_at(N - 1, k) %>=% tar.power)) N <- N - 1
      out[k] <- as.integer(N)
      limit[k] <- as.integer(L)
      next
    }
    Power <- power_at(N, k)
    if (Power %>=% tar.power) {
      while ((Power %>=% tar.power) && (N > 1)) {
        N <- N - 1
        Power <- power_at(N, k)
      }
      if (Power %<<% tar.power) N <- N + 1
    } else {
      while (Power %<<% tar.power) {
        N <- N + 1
        Power <- power_at(N, k)
      }
    }
    out[k] <- as.integer(N)
  }
  list(N = out, limit = limit)
}
