#' Power of a Blinded Sample Size Re-estimation Design for Binary Endpoints
#'
#' Calculates the power of a two-arm trial with a binary endpoint when blinded sample size
#' re-estimation (BSSR) is implemented, together with the power of the corresponding
#' fixed-sample design. Five exact tests are supported, each of which can be applied with a
#' one-sided or a two-sided alternative, under either a restricted or an unrestricted
#' design rule.
#'
#' @param p Vector of true pooled proportions of responders from both groups
#' @param Delta.A Assumed treatment effect (risk difference)
#' @param Delta.T True treatment effect (risk difference)
#' @param N1 Initial sample size of group 1
#' @param N2 Initial sample size of group 2
#' @param omega Fraction of the initial sample size observed at the interim analysis. The
#'   interim size of group 2 is \code{ceiling(omega N2)} and that of group 1 is
#'   \code{ceiling(r ceiling(omega N2))}, so the interim analysis keeps the allocation
#'   ratio
#' @param r Allocation ratio to group 1
#' @param alpha Level of significance for the alternative specified by \code{alternative}
#' @param tar.power Target power
#' @param Test Type of statistical test. Options: \code{'Chisq'}, \code{'Fisher'},
#'   \code{'Fisher-midP'}, \code{'Z-pool'}, or \code{'Boschloo'}
#' @param restricted Logical. If \code{TRUE}, the re-estimated sample size is not allowed
#'   to fall below the initial sample size. Default is \code{FALSE}
#' @param alternative Direction of the alternative hypothesis. Options: \code{'greater'}
#'   (default) or \code{'two.sided'}
#' @param tsmethod Convention used to construct the two-sided version of the conditional
#'   tests. Options: \code{'minlike'} (default) or \code{'central'}
#' @param n.grid Number of grid points used to search over the nuisance parameter of the
#'   unconditional tests. Default is 100
#' @param bb.gamma Confidence level parameter of the Berger-Boos procedure. The default of
#'   0 disables the procedure
#'
#' @return An object of class \code{bbssr_powerbssr}, a data frame with one row per element
#'   of \code{p} containing:
#' \describe{
#'   \item{p1}{True probability of responders for group 1}
#'   \item{p2}{True probability of responders for group 2}
#'   \item{p}{True pooled probability of responders from both groups}
#'   \item{power.BSSR}{Power of the BSSR design}
#'   \item{power.TRAD}{Power of the fixed-sample design}
#'   \item{E.N}{Expected total sample size of the BSSR design}
#' }
#' The interim sample sizes are stored as the attributes \code{n1.interim} and
#' \code{n2.interim}.
#'
#' @details
#' Both the interim and the final sample size of group 1 are obtained from the size of
#' group 2 by a single application of \code{ceiling(r ...)}, so the allocation ratio is
#' preserved as closely as whole numbers allow and is exact whenever \code{r} is a whole
#' number. The size of the second stage follows as the difference between the final and
#' the interim size, which keeps the two stages adding up to the final size. The argument
#' \code{N1} enters only through the fixed-sample comparator, and a warning is issued when
#' it is not \code{ceiling(r N2)}.
#'
#' At the interim analysis the pooled number of responders is observed without unblinding.
#' The pooled proportion is combined with the assumed treatment effect \code{Delta.A} to
#' recover group-specific proportions, from which the sample size is re-estimated. The power
#' is then averaged over the distribution of the interim outcome.
#'
#' Setting \code{Delta.T} to 0 makes the two groups identical, so \code{power.BSSR} and
#' \code{power.TRAD} become rejection probabilities under the null hypothesis. They are
#' evaluated only at the values supplied through \code{p}, since the function does not
#' search over the unit interval on its own. The largest type I error rate of the design is
#' therefore obtained by passing a grid such as \code{p = seq(0.02, 0.98, by = 0.02)} and
#' taking the maximum of the \code{power.BSSR} column. The rejection region and the
#' re-estimated sample size are shared across the elements of \code{p}, so a fine grid
#' costs much less than the same number of separate calls.
#'
#' The weighted approach available in earlier versions of the package has been removed.
#'
#' @examples
#' # Small BSSR calculation with the chi-squared test
#' BinaryPowerBSSR(
#'   p = 0.45,
#'   Delta.A = 0.3, Delta.T = 0.3,
#'   N1 = 5, N2 = 5, omega = 0.5, r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
#' )
#'
#' \donttest{
#' res <- BinaryPowerBSSR(
#'   p = seq(0.19, 0.37, by = 0.03),
#'   Delta.A = 0.36, Delta.T = 0.36,
#'   N1 = 24, N2 = 24, omega = 0.5, r = 1,
#'   alpha = 0.025, tar.power = 0.8, Test = 'Z-pool'
#' )
#' print(res)
#' plot(res)
#' }
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @export
#' @import fpCompare
#' @importFrom stats dbinom
BinaryPowerBSSR <- function(p, Delta.A, Delta.T, N1, N2, omega, r,
                            alpha, tar.power, Test, restricted = FALSE,
                            alternative = c('greater', 'two.sided'),
                            tsmethod = c('minlike', 'central'),
                            n.grid = 100, bb.gamma = 0) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  if (length(omega) != 1 || omega <= 0 || omega > 1) {
    stop('omega must be a single value in (0, 1]')
  }
  if (N1 != ceiling(r * N2)) {
    warning('N1 differs from ceiling(r * N2), the fixed-sample comparator is not ',
            'allocated in the ratio r to 1')
  }
  # Interim sample sizes. The experimental group is obtained from the control group by a
  # single rounding, so the interim allocation is as close to r to 1 as integers allow and
  # equals r times the control group exactly whenever r is a whole number
  N12 <- as.integer(ceiling(omega * N2))
  N11 <- as.integer(ceiling(r * N12))
  # Estimate the pooled proportion in a blinded fashion
  hat.p <- c(outer(0:N11, 0:N12, '+') / (N11 + N12))
  # Derive hat{p}_{1} and hat{p}_{2} using hat{p} and Delta.A
  hat.p1 <- pmin(1, hat.p + (1 / (1 + r)) * Delta.A)
  hat.p2 <- pmax(0, hat.p - (r / (1 + r)) * Delta.A)
  # Sample size re-estimation, evaluated once per distinct pair of interim proportions
  key.p <- paste(hat.p1, hat.p2, sep = '_')
  uniq.p <- !duplicated(key.p)
  ss.cache <- new.env(parent = emptyenv())
  for (u in which(uniq.p)) {
    assign(
      key.p[u],
      BinarySampleSize(hat.p1[u], hat.p2[u], r, alpha, tar.power, Test,
                       alternative, tsmethod, n.grid, bb.gamma)[['N2']],
      envir = ss.cache
    )
  }
  BSSR.N2 <- vapply(key.p, function(k) as.numeric(get(k, envir = ss.cache)),
                    numeric(1), USE.NAMES = FALSE)
  # Final sample sizes, again with a single rounding of the control group. The second
  # stage follows as the difference, so the two stages add up to the final size exactly
  hat.N2 <- if (restricted == TRUE) {
    as.integer(pmax(N2, BSSR.N2))
  } else {
    as.integer(pmax(N12, BSSR.N2))
  }
  hat.N1 <- as.integer(ceiling(r * hat.N2))
  N22 <- hat.N2 - N12
  N21 <- hat.N1 - N11
  hat.N <- hat.N1 + hat.N2
  # True proportion for each treatment group
  p1 <- p + (1 / (1 + r)) * Delta.T
  p2 <- p - (r / (1 + r)) * Delta.T
  # Omit scenarios with p_{j} outside the unit interval
  omit.ID <- union(which(p1 %>>% 1), which(p2 %<<% 0))
  if (length(omit.ID) %!=% 0) {
    p1 <- p1[-omit.ID]
    p2 <- p2[-omit.ID]
    p <- p[-omit.ID]
  }
  if (length(p) == 0) stop('no scenario has both p1 and p2 inside the unit interval')
  # Probability mass functions of the interim outcome
  dbinom1 <- outer(X = 0:N11, Y = p1, function(X, Y) dbinom(X, N11, Y))
  dbinom2 <- outer(X = 0:N12, Y = p2, function(X, Y) dbinom(X, N12, Y))
  # Interim responder counts, in the same order as hat.p
  x11 <- c(row(matrix(0L, nrow = N11 + 1L, ncol = N12 + 1L)) - 1L)
  x12 <- c(col(matrix(0L, nrow = N11 + 1L, ncol = N12 + 1L)) - 1L)
  # Rejection regions, evaluated once per distinct pair of final sample sizes
  key.N <- paste(hat.N1, hat.N2, sep = '_')
  uniq.N <- !duplicated(key.N)
  rr.cache <- new.env(parent = emptyenv())
  for (u in which(uniq.N)) {
    RR <- BinaryRR(hat.N1[u], hat.N2[u], alpha, Test, alternative, tsmethod, n.grid, bb.gamma)
    assign(key.N[u], matrix(as.vector(RR), nrow = hat.N1[u] + 1L, ncol = hat.N2[u] + 1L),
           envir = rr.cache)
  }
  # Conditional power of the second stage given each interim outcome
  power.stage2 <- do.call(
    cbind,
    lapply(seq_along(hat.p), function(i) {
      rr.full <- get(key.N[i], envir = rr.cache)
      rows <- (x11[i] + 1L):(x11[i] + 1L + N21[i])
      cols <- (x12[i] + 1L):(x12[i] + 1L + N22[i])
      rr <- matrix(rr.full[rows, cols], nrow = length(rows), ncol = length(cols))
      vapply(
        seq_along(p),
        function(j) power_from_rr(rr, dbinom(0:N21[i], N21[i], p1[j]),
                                  dbinom(0:N22[i], N22[i], p2[j])),
        numeric(1)
      )
    })
  )
  # Average over the distribution of the interim outcome
  power.BSSR <- vapply(seq_along(p), function(k) {
    sum(c(dbinom1[, k] %o% dbinom2[, k]) * power.stage2[k, ], na.rm = TRUE)
  }, numeric(1))
  E.N <- vapply(seq_along(p), function(k) {
    sum(c(dbinom1[, k] %o% dbinom2[, k]) * hat.N)
  }, numeric(1))
  # Power of the fixed-sample design
  power.TRAD <- BinaryPower(p1, p2, N1, N2, alpha, Test, alternative, tsmethod,
                            n.grid, bb.gamma)$Power
  out <- data.frame(p1, p2, p, power.BSSR, power.TRAD, E.N)
  attr(out, 'Test') <- Test
  attr(out, 'alternative') <- alternative
  attr(out, 'alpha') <- alpha
  attr(out, 'tar.power') <- tar.power
  attr(out, 'restricted') <- restricted
  attr(out, 'N1') <- N1
  attr(out, 'N2') <- N2
  attr(out, 'omega') <- omega
  attr(out, 'n1.interim') <- N11
  attr(out, 'n2.interim') <- N12
  attr(out, 'Delta.A') <- Delta.A
  attr(out, 'Delta.T') <- Delta.T
  class(out) <- c('bbssr_powerbssr', 'data.frame')
  out
}
