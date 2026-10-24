#' Validate the Settings of a Random Allocation
#'
#' Internal helper checking the fractions of the strata, the allocation ratios, the total
#' sample size and the block size of a two-stratum design whose stratum and group sizes are
#' random, that is \code{alloc} other than \code{'fixed'} in \code{\link{StratPower}}.
#'
#' @param alloc \code{'stratum'}, \code{'group'} or \code{'both'}
#' @param N Total sample size, or \code{NULL} when it is not yet known
#' @param frac Fractions of the strata, one positive element per stratum
#' @param r Allocation ratio of group 1 to group 2, a single value or one per stratum
#' @param block Size of the permuted blocks, used with \code{alloc = 'stratum'}
#'
#' @return A list with the fractions \code{frac}, summing to 1, the allocation ratios
#'   \code{r}, the probabilities \code{b = r / (1 + r)} of group 1 within the strata, the
#'   numbers \code{g} of patients of group 1 in a block and the integer \code{block}
#'   (both \code{NULL} unless \code{alloc = 'stratum'})
#'
#' @keywords internal
#' @noRd
check_strat_alloc <- function(alloc, N, frac, r, block) {
  if (is.null(frac)) stop('frac must be supplied with a random allocation')
  if (length(frac) != 2) {
    stop('a random allocation is available for designs with two strata')
  }
  fr <- check_strat_frac(frac, r, 2)
  if (!is.null(N) && (length(N) != 1 || !is.numeric(N) || is.na(N) || N < 1 ||
                      N != round(N))) {
    stop('N must be a single positive integer')
  }
  b <- fr$r / (1 + fr$r)
  g <- NULL
  if (alloc == 'stratum') {
    if (length(block) != 1 || !is.numeric(block) || is.na(block) || block < 2 ||
        block != round(block)) {
      stop('block must be a single integer of at least 2')
    }
    g <- block * b
    if (any(abs(g - round(g)) > 1e-9) || any(round(g) < 1) || any(round(g) >= block)) {
      stop('block must hold whole numbers of patients of the two groups in the ratio r, ',
           'such as an even block with r = 1')
    }
    g <- as.integer(round(g))
  }
  list(frac = fr$frac, r = fr$r, b = b, g = g,
       block = if (alloc == 'stratum') as.integer(block) else NULL)
}
