#' Unrounded Sample Size of Group 2 for a Test of Equivalence
#'
#' Internal helper returning the unrounded size of group 2 at which the normal
#' approximation to the power of the two one-sided tests reaches the target power. Group
#' 1 receives \code{r} times as many patients. Let \code{d_L}, \code{v0_L} and
#' \code{v1_L} be the ingredients of \code{ss_raw_parts} for the lower boundary
#' (alternative \code{'greater'} with the margin \code{lower}) and \code{d_U},
#' \code{v0_U} and \code{v1_U} those for the upper boundary (alternative \code{'less'}
#' with the margin \code{upper}). Each one-sided test is at the level \code{alpha}, and
#' the approximate power at the size \code{n} of group 2 is
#' \code{Phi((d_L sqrt(n) - z.a sqrt(v0_L)) / sqrt(v1_L)) +
#' Phi((d_U sqrt(n) - z.a sqrt(v0_U)) / sqrt(v1_U)) - 1}, with \code{z.a} the upper
#' \code{alpha} quantile of the standard normal distribution. Each term is the normal
#' approximation to the power of one of the two tests that underlies the one-sided size
#' of \code{ss_raw_n2}, and the sum minus 1 is a lower bound for the probability that
#' both reject. The approximation increases with \code{n}, and the size is the root of
#' the approximation minus the target power. A term whose variance \code{v1} is 0 is 1
#' when its numerator is positive and 0 otherwise. A pair of probabilities outside the
#' equivalence region, or on its boundary, gives \code{NaN}, and non-finite ingredients
#' give \code{NA}. When the approximation reaches the target at no patients, which
#' happens only with \code{alpha} above 0.5, the size is 0. The function is vectorized
#' over \code{p1} and \code{p2}.
#'
#' @param p1 Response probability of group 1
#' @param p2 Response probability of group 2
#' @param r Allocation ratio to group 1
#' @param alpha Level of significance of each one-sided test
#' @param tar.power Target power
#' @param ss.method \code{'standard'}, \code{'null.variance'} or
#'   \code{'alternative.variance'}
#' @param margin Margins of the test of equivalence, see \code{norm_margin}
#' @param margin.scale \code{'RD'} or \code{'RR'}, the scale of the margins
#'
#' @return A numeric vector
#'
#' @keywords internal
#' @noRd
#' @importFrom stats pnorm qnorm uniroot
ss_raw_n2_equiv <- function(p1, p2, r, alpha, tar.power, ss.method, margin,
                            margin.scale) {
  m <- norm_margin('equivalence', margin, margin.scale)
  z.a <- qnorm(1 - alpha)
  n <- max(length(p1), length(p2))
  lo <- lapply(ss_raw_parts(p1, p2, r, 'greater', ss.method, m[1], margin.scale),
               rep_len, length.out = n)
  up <- lapply(ss_raw_parts(p1, p2, r, 'less', ss.method, m[2], margin.scale),
               rep_len, length.out = n)
  term <- function(d, v0, v1, x) {
    num <- d * sqrt(x) - z.a * sqrt(v0)
    if (v1 > 0) pnorm(num / sqrt(v1)) else as.numeric(num > 0)
  }
  vapply(seq_len(n), function(i) {
    v <- c(lo$d[i], lo$v0[i], lo$v1[i], up$d[i], up$v0[i], up$v1[i])
    if (anyNA(v) || any(!is.finite(v))) return(NA_real_)
    if (!(lo$d[i] > 0 && up$d[i] > 0)) return(NaN)
    gap <- function(x) {
      term(lo$d[i], lo$v0[i], lo$v1[i], x) + term(up$d[i], up$v0[i], up$v1[i], x) - 1 -
        tar.power
    }
    # With alpha above 0.5 the approximation can exceed the target without patients
    if (gap(0) >= 0) return(0)
    hi <- 1
    while (gap(hi) < 0) hi <- 2 * hi
    uniroot(gap, c(0, hi), tol = 1e-12)$root
  }, numeric(1))
}
