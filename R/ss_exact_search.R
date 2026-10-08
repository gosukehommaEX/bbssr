#' Exact Sample Size Search for Several Pairs of Response Probabilities
#'
#' Internal helper returning, for each pair of response probabilities, the size of
#' group 2 found by the search of \code{\link{BinarySampleSize}} with
#' \code{method = 'exact'}. Group 1 receives \code{ceiling(r N2)} patients.
#'
#' With \code{search = 'crossing'} the search for each pair starts from the normal
#' approximation, lowers the size of group 2 one unit at a time as long as the exact power
#' attains the target power, or otherwise raises it one unit at a time until the power
#' attains it. The returned size attains the target power while the size one unit smaller
#' does not, unless the returned size is 1, and need not be the smallest size attaining
#' the target power, since the exact power is not monotone in the sample size. With
#' \code{search = 'smallest'} the sizes are scanned upwards from 1 and the first size
#' attaining the target power is returned. With \code{search = 'stable'} the sizes are
#' scanned downwards from the limit \code{max(ceiling(a n0), ceiling(n0 + b))}, where
#' \code{n0} is the starting size of \code{'crossing'} and \code{a} and \code{b} are the
#' two elements of \code{search.limit}, and the returned size is the smallest from which
#' every size up to the limit attains the target power. The search stops with an error if
#' the limit itself does not attain it. The rejection region of each size visited is
#' obtained once and reused by every pair whose search visits that size, and the power of
#' a pair is computed by \code{power_from_rr}, so each pair receives the same sample size
#' as a separate call of \code{\link{BinarySampleSize}} would give.
#'
#' @param p1 Response probabilities of group 1
#' @param p2 Response probabilities of group 2, of the same length as \code{p1}, with no
#'   element equal to the corresponding element of \code{p1}
#' @param r Allocation ratio to group 1
#' @param alpha Level of significance
#' @param tar.power Target power
#' @param Test Type of statistical test
#' @param alternative \code{'greater'}, \code{'less'} or \code{'two.sided'}
#' @param tsmethod \code{'minlike'}, \code{'central'} or \code{'blaker'}
#' @param n.grid Number of grid points over the nuisance parameter
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure
#' @param ref.pvalue Logical. Whether the maximum over the nuisance parameter of the
#'   unconditional tests is refined between the grid points
#' @param margin Non-inferiority margin, 0 for a test of superiority on the scale of the
#'   risk difference. With a margin other than 0, or a margin on the scale of the risk
#'   ratio, the search starts from the formula of Farrington and Manning (1990)
#' @param search \code{'crossing'} (default), \code{'smallest'} or \code{'stable'}
#' @param search.limit The factor \code{a} and the increment \code{b} of the limit of
#'   \code{search = 'stable'}. Default is \code{c(2, 50)}
#' @param margin.scale \code{'RD'} or \code{'RR'}, the scale of the margin
#'
#' @return A list with the integer vectors \code{N2}, the sizes of group 2, and
#'   \code{limit}, the limits of \code{search = 'stable'} (otherwise \code{NA})
#'
#' @keywords internal
#' @noRd
#' @import fpCompare
#' @importFrom stats qnorm dbinom
ss_exact_search <- function(p1, p2, r, alpha, tar.power, Test, alternative, tsmethod,
                            n.grid, bb.gamma, ref.pvalue, margin, search = 'crossing',
                            search.limit = c(2, 50), margin.scale) {
  if (search == 'stable' &&
      (!is.numeric(search.limit) || length(search.limit) != 2 || anyNA(search.limit) ||
         any(!is.finite(search.limit)) || search.limit[1] < 1 || search.limit[2] < 0)) {
    stop('search.limit must be two numbers, a factor of at least 1 and a non-negative ',
         'increment')
  }
  store <- new.env(parent = emptyenv())
  # Rejection region at a given size of group 2, obtained once per search
  rr_at <- function(N2) {
    key <- as.character(N2)
    hit <- store[[key]]
    if (is.null(hit)) {
      hit <- get_rr(ceiling(r * N2), N2, alpha, Test, alternative, tsmethod, n.grid,
                    bb.gamma, ref.pvalue, margin, margin.scale)
      assign(key, hit, envir = store)
    }
    hit
  }
  # Exact power of pair k at a given size of group 2
  power_at <- function(N2, k) {
    N1 <- ceiling(r * N2)
    power_from_rr(rr_at(N2), dbinom(0:N1, N1, p1[k]), dbinom(0:N2, N2, p2[k]))
  }
  # Step 0 (initial sample size from the normal approximation to the chi-squared test)
  alpha.eff <- if (alternative == 'two.sided') alpha / 2 else alpha
  p <- (r * p1 + p2) / (1 + r)
  init.N2 <- '*'(
    (1 + 1 / r) / ((p1 - p2) ^ 2),
    (qnorm(alpha.eff) * sqrt(p * (1 - p)) +
       qnorm(1 - tar.power) * sqrt((p1 * (1 - p1) / r + p2 * (1 - p2)) / (1 + 1 / r))) ^ 2
  )
  if (margin.scale == 'RR' || margin != 0) {
    init.N2 <- ss_raw_n2(p1, p2, r, alpha, tar.power, alternative, 'standard', margin,
                         margin.scale)
  }
  out <- integer(length(p1))
  limit <- rep(NA_integer_, length(p1))
  for (k in seq_along(p1)) {
    N2 <- max(1, ceiling(init.N2[k]))
    if (search == 'smallest') {
      # Smallest size attaining the target power, by a scan upwards from 1
      N2 <- 1
      while (power_at(N2, k) %<<% tar.power) N2 <- N2 + 1
      out[k] <- as.integer(N2)
      next
    }
    if (search == 'stable') {
      # Smallest size from which every size up to the limit attains the target power, by a
      # scan downwards from the limit
      L <- max(ceiling(search.limit[1] * N2), ceiling(N2 + search.limit[2]))
      if (power_at(L, k) %<<% tar.power) {
        stop('the size ', L, ' of group 2, the limit of the stable search, does not ',
             'attain the target power; raise search.limit')
      }
      N2 <- L
      while (N2 > 1 && (power_at(N2 - 1, k) %>=% tar.power)) N2 <- N2 - 1
      out[k] <- as.integer(N2)
      limit[k] <- as.integer(L)
      next
    }
    # Step 1 (power calculation given the initial sample size)
    Power <- power_at(N2, k)
    # Step 2 (sample size calculation via a grid search algorithm)
    if (Power %>=% tar.power) {
      while ((Power %>=% tar.power) && (N2 > 1)) {
        N2 <- N2 - 1
        Power <- power_at(N2, k)
      }
      if (Power %<<% tar.power) N2 <- N2 + 1
    } else {
      while (Power %<<% tar.power) {
        N2 <- N2 + 1
        Power <- power_at(N2, k)
      }
    }
    out[k] <- as.integer(N2)
  }
  list(N2 = out, limit = limit)
}
