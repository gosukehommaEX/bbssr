#' Rejection Region of a Stratified Test
#'
#' Returns the critical values of a test for a two-arm trial with a binary endpoint
#' stratified into two or more strata, for every combination of the numbers of responders
#' of the strata. The tests are based on \code{S}, the number of responders of group 1
#' summed over the strata. Given the numbers of responders \code{z = (z_1, ..., z_J)} of
#' the strata, each test rejects the null hypothesis when \code{S} is at most a lower or
#' at least an upper critical value.
#'
#' @param N1 Sizes of group 1 in the strata, one element per stratum
#' @param N2 Sizes of group 2 in the strata, one element per stratum
#' @param alpha Level of significance, for the alternative specified
#' @param Test Statistical test. Options: \code{'CMH'} (Cochran-Mantel-Haenszel test
#'   with the hypergeometric variance, without continuity correction), \code{'CMH-cc'}
#'   (the same with the continuity correction of \code{\link[stats]{mantelhaen.test}}),
#'   \code{'Cochran'} (Cochran's test with the binomial variance, the test of Woolson et
#'   al. (1986)), \code{'Cochran-cc'} (the same with the continuity correction, which
#'   rejects the same outcomes as that of Nam (1992) at one-sided levels up to 0.078) or
#'   \code{'Strat-Fisher'} (the stratified exact conditional test of Jung (2014)). See
#'   Details
#' @param alternative Direction of the alternative hypothesis. Options:
#'   \code{'greater'} (default), under which group 1 responds more often than group 2,
#'   \code{'less'} or \code{'two.sided'}
#' @param tsmethod Convention of the two-sided p-value of \code{'Strat-Fisher'}, as for
#'   the Fisher exact test of \code{\link{BinaryRR}}. Options: \code{'minlike'}
#'   (default), \code{'central'} or \code{'blaker'}. Ignored by the other tests and for
#'   one-sided alternatives
#'
#' @return An object of class \code{bbssr_srr}, a data frame with one row for each
#'   combination of the numbers of responders of the strata, the first stratum varying
#'   fastest, containing:
#' \describe{
#'   \item{z.S1, z.S2, ...}{Numbers of responders of the strata}
#'   \item{S.min, S.max}{Smallest and largest possible value of \code{S}}
#'   \item{lower}{Lower critical value: the null hypothesis is rejected when
#'     \code{S <= lower}. \code{NA} when no value of \code{S} is rejected in the lower
#'     tail}
#'   \item{upper}{Upper critical value: the null hypothesis is rejected when
#'     \code{S >= upper}. \code{NA} when no value of \code{S} is rejected in the upper
#'     tail}
#' }
#'
#' @details
#' In stratum \code{j}, with \code{N_j = N1[j] + N2[j]} patients and \code{z_j}
#' responders, the null expectation of the number of responders of group 1 is
#' \code{N1[j] z_j / N_j}, and \code{E} is the sum of these over the strata. The
#' Cochran-Mantel-Haenszel statistic is \code{(S - E) / sqrt(V)} with the hypergeometric
#' variance
#' \code{V = sum_j N1[j] N2[j] z_j (N_j - z_j) / (N_j^2 (N_j - 1))}, and Cochran's
#' statistic uses the binomial variance \code{sum_j N1[j] N2[j] z_j (N_j - z_j) / N_j^3}
#' instead. The continuity correction reduces \code{|S - E|} by 0.5 when it is at least 0.5,
#' as in \code{\link[stats]{mantelhaen.test}}. Nam (1992) reduces \code{S - E} by 0.5 in
#' every outcome, which makes a difference only where \code{|S - E| < 0.5}. There \code{S}
#' is the integer next to \code{E}, and since the hypergeometric variance is at least
#' \code{f (1 - f)}, with \code{f} the fractional part of \code{E}, and the binomial
#' variance at least half of it, the one-sided p-value of \code{'Cochran-cc'} exceeds
#' 0.078. The two corrections therefore reject the same outcomes at one-sided levels up to
#' 0.078. The statistics are referred to the standard
#' normal distribution, with the absolute value for a two-sided alternative, and an
#' outcome whose variance is 0 is not rejected. The stratified exact conditional test uses
#' the distribution of \code{S} given \code{z}, the convolution of the hypergeometric
#' distributions of the strata. The two-sided p-value \code{'central'} is twice the smaller
#' tail probability (Jung, 2014), and \code{'minlike'} and \code{'blaker'} order the values
#' of \code{S} by their probability and by their smaller tail probability, with the ties
#' of \code{\link{BinaryRR}}. As in the other functions of the package, a p-value is
#' rejected when it is below \code{alpha} by more than about 1.49e-8.
#'
#' The function stores the critical values for reuse during the session, like the
#' rejection regions of the other designs.
#'
#' @references
#' Woolson RF, Bean JA, Rojas PB (1986). Sample size for case-control studies using
#' Cochran's statistic. \emph{Biometrics}, 42(4), 927-932.
#'
#' Nam J (1992). Sample size determination for case-control studies and the comparison of
#' stratified and unstratified analyses. \emph{Biometrics}, 48(2), 389-395.
#'
#' Jung SH (2014). Stratified Fisher's exact test and its sample size calculation.
#' \emph{Biometrical Journal}, 56(1), 129-140.
#'
#' @examples
#' rr <- StratRR(N1 = c(10, 12), N2 = c(10, 12), alpha = 0.025, Test = 'CMH')
#' print(rr)
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{StratPower}}, \code{\link{StratSampleSize}}
#' @export
StratRR <- function(N1, N2, alpha, Test, alternative = c('greater', 'less', 'two.sided'),
                    tsmethod = c('minlike', 'central', 'blaker')) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  a <- check_strat_args(N1, N2, alpha, Test, alternative, tsmethod)
  J <- length(a$N1)
  if (prod(as.double(a$N1 + a$N2 + 1)) > 1e6) {
    stop('the design has more than 1e6 combinations of the numbers of responders of the ',
         'strata')
  }
  crit <- get_strat_crit(a$N1, a$N2, alpha, a$Test, a$alternative, a$tsmethod)
  z <- as.matrix(expand.grid(lapply(a$N1 + a$N2, function(n) 0:n)))
  colnames(z) <- paste0('z.S', seq_len(J))
  S.min <- as.integer(rowSums(pmax(sweep(z, 2, a$N2), 0)))
  S.max <- as.integer(rowSums(pmin(z, matrix(a$N1, nrow(z), J, byrow = TRUE))))
  lower <- ifelse(crit$lo >= S.min, crit$lo, NA_integer_)
  upper <- ifelse(crit$hi <= S.max, crit$hi, NA_integer_)
  out <- data.frame(z, S.min = S.min, S.max = S.max, lower = as.integer(lower),
                    upper = as.integer(upper))
  attr(out, 'N1') <- a$N1
  attr(out, 'N2') <- a$N2
  attr(out, 'alpha') <- alpha
  attr(out, 'Test') <- a$Test
  attr(out, 'alternative') <- a$alternative
  attr(out, 'tsmethod') <- a$tsmethod
  class(out) <- c('bbssr_srr', 'data.frame')
  out
}
