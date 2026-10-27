#' Exact Power of a Stratified Test
#'
#' Calculates the exact power of a test of \code{\link{StratRR}} for a two-arm trial with a
#' binary endpoint, with fixed group sizes in each stratum or, with two strata, with
#' stratum and group sizes that are random.
#'
#' @param p1 Response probabilities of group 1 in the strata: a vector with one element
#'   per stratum, or a matrix with one row per scenario and one column per stratum
#' @param p2 Response probabilities of group 2 in the strata, of the same form as
#'   \code{p1}
#' @param N1 Sizes of group 1 in the strata, one element per stratum, with
#'   \code{alloc = 'fixed'}
#' @param N2 Sizes of group 2 in the strata, one element per stratum, with
#'   \code{alloc = 'fixed'}
#' @param alloc Allocation of the patients to the strata and the groups. \code{'fixed'}
#'   (default) uses the group sizes \code{N1} and \code{N2}. The others need two strata
#'   and make the sizes random, see Details: \code{'stratum'} (stratum sizes binomial,
#'   permuted blocks within the strata), \code{'group'} (stratum sizes fixed, group sizes
#'   binomial) and \code{'both'} (stratum and group sizes binomial)
#' @param N Total sample size, with a random allocation
#' @param frac Fractions of the two strata, with a random allocation. They are rescaled to
#'   sum to 1
#' @param r Allocation ratio of group 1 to group 2 within the strata, a single value or
#'   one value per stratum, with a random allocation. Default is 1
#' @param block Size of the permuted blocks of \code{alloc = 'stratum'}. Each block holds
#'   \code{block r / (1 + r)} patients of group 1, which must be a whole number. Default
#'   is 4
#' @inheritParams StratRR
#'
#' @return An object of class \code{bbssr_spower}, a data frame with one row per scenario
#'   containing the response probabilities \code{p1.S1, p1.S2, ...} and
#'   \code{p2.S1, p2.S2, ...} and the exact power \code{Power}. With a random allocation
#'   the attribute \code{skipped} is the total probability of the configurations left out
#'   of the sum, see Details, and \code{n.strata} the stratum sizes of
#'   \code{alloc = 'group'}.
#'
#' @details
#' The numbers of responders of the two groups in each stratum are independent binomial
#' variables. The power is the probability of the rejection region of
#' \code{\link{StratRR}}, computed exactly by enumerating the numbers of responders of the
#' strata and carrying the distribution of the partial sum of \code{S} over the strata.
#' The response probabilities may differ between the strata in any way, so the power is
#' available also when the odds ratio or the risk difference is not common to the strata.
#'
#' With \code{alloc} other than \code{'fixed'} the design has two strata and a total of
#' \code{N} patients, and the sizes of the strata and of the groups within them are
#' random. With \code{'group'} the stratum sizes are fixed at \code{n1 = floor(N frac[1])}
#' and \code{n2 = N - n1}, and each patient of stratum \code{j} is in group 1 with the
#' probability \code{r[j] / (1 + r[j])}, so that the size of group 1 is binomial (Jung,
#' 2014, Section 3.2, where the size of stratum 1 is the round-off of \code{N frac[1]};
#' the integer part used here reproduces all but one of the sample sizes of the column
#' 'n fixed' of Table 4 of the article). With \code{'both'} the size of stratum 1 is
#' binomial as well, with
#' the probability \code{frac[1]} (Jung, 2014, Section 3.1). With \code{'stratum'} the
#' stratum sizes are binomial and the groups are allocated by permuted blocks of size
#' \code{block} within each stratum: every complete block holds
#' \code{block r[j] / (1 + r[j])} patients of group 1, and the number of patients of
#' group 1 among the patients of the last, incomplete block is hypergeometric. The power is
#' the sum over the stratum and group sizes of their probability times the power of the
#' fixed design they give. Configurations with a probability below 1e-14 are left out, and
#' their total probability, which bounds the error of the power, is returned in the
#' attribute \code{skipped}. A group without patients in a stratum leaves that stratum
#' without variation of \code{S}, and a stratum without patients leaves a design with one
#' stratum. The number of configurations grows with \code{N} most quickly for
#' \code{'both'}: at \code{N = 100} the power of one scenario takes a fraction of a second
#' with \code{'stratum'}, about a second with \code{'group'} and tens of seconds with
#' \code{'both'}.
#'
#' @references
#' Jung SH (2014). Stratified Fisher's exact test and its sample size calculation.
#' \emph{Biometrical Journal}, 56(1), 129-140.
#'
#' @examples
#' StratPower(p1 = c(0.5, 0.7), p2 = c(0.3, 0.5), N1 = c(40, 40), N2 = c(40, 40),
#'            alpha = 0.025, Test = 'CMH')
#'
#' # Several scenarios, and the stratified exact test
#' StratPower(p1 = rbind(c(0.5, 0.7), c(0.45, 0.65)), p2 = rbind(c(0.3, 0.5), c(0.3, 0.5)),
#'            N1 = c(40, 40), N2 = c(40, 40), alpha = 0.025, Test = 'Strat-Fisher')
#'
#' # Random stratum sizes, with permuted blocks of 4 within the strata
#' StratPower(p1 = c(0.5, 0.7), p2 = c(0.3, 0.5), alpha = 0.025, Test = 'CMH',
#'            alloc = 'stratum', N = 160, frac = c(0.5, 0.5))
#'
#' @author Gosuke Homma (\email{my.name.is.gosuke@@gmail.com})
#' @seealso \code{\link{StratRR}}, \code{\link{StratSampleSize}}
#' @export
StratPower <- function(p1, p2, N1 = NULL, N2 = NULL, alpha, Test,
                       alternative = c('greater', 'less', 'two.sided'),
                       tsmethod = c('minlike', 'central', 'blaker'),
                       alloc = c('fixed', 'stratum', 'group', 'both'), N = NULL,
                       frac = NULL, r = 1, block = 4) {
  alternative <- match.arg(alternative)
  tsmethod <- match.arg(tsmethod)
  alloc <- match.arg(alloc)
  if (alloc == 'fixed') {
    if (is.null(N1) || is.null(N2)) stop("N1 and N2 must be supplied with alloc = 'fixed'")
    if (!is.null(N) || !is.null(frac)) {
      stop("N and frac describe a random allocation; with alloc = 'fixed' the design is ",
           'given by N1 and N2')
    }
    a <- check_strat_args(N1, N2, alpha, Test, alternative, tsmethod)
  } else {
    if (!is.null(N1) || !is.null(N2)) {
      stop('with a random allocation the design is given by N, frac and r, not by N1 ',
           'and N2')
    }
    if (is.null(N)) stop('N must be supplied with a random allocation')
    al <- check_strat_alloc(alloc, N, frac, r, block)
    a <- check_strat_args(c(1, 1), c(1, 1), alpha, Test, alternative, tsmethod)
  }
  J <- if (alloc == 'fixed') length(a$N1) else 2L
  p1 <- check_strat_probs(p1, J, 'p1')
  p2 <- check_strat_probs(p2, J, 'p2')
  if (nrow(p1) != nrow(p2)) {
    if (nrow(p1) == 1) p1 <- p1[rep(1, nrow(p2)), , drop = FALSE]
    else if (nrow(p2) == 1) p2 <- p2[rep(1, nrow(p1)), , drop = FALSE]
    else stop('p1 and p2 must have the same number of scenarios')
  }
  if (alloc == 'fixed') {
    Power <- strat_fixed_power(a$N1, a$N2, p1, p2, alpha, a$Test, a$alternative,
                               a$tsmethod)
  } else {
    res <- strat_alloc_power(N, al, alloc, p1, p2, alpha, a$Test, a$alternative,
                             a$tsmethod)
    Power <- res$power
  }
  out <- data.frame(p1, p2, Power = Power)
  names(out) <- c(paste0('p1.S', seq_len(J)), paste0('p2.S', seq_len(J)), 'Power')
  attr(out, 'alloc') <- alloc
  if (alloc == 'fixed') {
    attr(out, 'N1') <- a$N1
    attr(out, 'N2') <- a$N2
  } else {
    attr(out, 'N') <- as.integer(N)
    attr(out, 'frac') <- al$frac
    attr(out, 'r') <- al$r
    attr(out, 'block') <- al$block
    attr(out, 'n.strata') <- res$n
    attr(out, 'skipped') <- res$skipped
  }
  attr(out, 'alpha') <- alpha
  attr(out, 'Test') <- a$Test
  attr(out, 'alternative') <- a$alternative
  attr(out, 'tsmethod') <- a$tsmethod
  class(out) <- c('bbssr_spower', 'data.frame')
  out
}
