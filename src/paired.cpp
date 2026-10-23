#include <Rcpp.h>
#include <vector>
#include "paired_core.h"
using namespace Rcpp;

// Rcpp interfaces of the computations for the paired design, see paired_core.h.

// Regions of a list, which must be square logical matrices, as pointers to their cells.
// The matrices are not copied, so the pointers stay valid while the list exists. Ns
// receives the numbers of pairs.
static std::vector<const int*> region_pointers(const List& rr_list, int n1,
                                               std::vector<int>& Ns) {
  std::vector<const int*> regions;
  Ns.clear();
  for (R_xlen_t k = 0; k < rr_list.size(); ++k) {
    SEXP s = rr_list[k];
    if (TYPEOF(s) != LGLSXP || !Rf_isMatrix(s)) {
      stop("the rejection regions must be logical matrices");
    }
    LogicalMatrix R(s);
    if (R.ncol() != R.nrow()) stop("the rejection regions must be square");
    if (R.nrow() - 1 < n1) stop("a final number of pairs is below the interim number");
    regions.push_back(R.begin());
    Ns.push_back(R.nrow() - 1);
  }
  return regions;
}

// Largest element of a vector of numbers of pairs, at least 0
static int max_pairs(const std::vector<int>& Ns) {
  int Nmax = 0;
  for (int N : Ns) if (N > Nmax) Nmax = N;
  return Nmax;
}

// Null tail probabilities of sorted cells maximized over the discordance probability.
// dbinomD (N + 1) by G matrix of b(D; N, theta[g]); D and w the number of discordant pairs
// and the null weight b(a; D, 1/2) of each sorted cell; idx_last, g_lo and g_hi zero
// based, as in max_tail_prob().
// [[Rcpp::export]]
NumericVector paired_tail_max(NumericMatrix dbinomD, IntegerVector D, NumericVector w,
                              IntegerVector idx_last, IntegerVector g_lo, IntegerVector g_hi,
                              NumericVector theta, bool refine) {
  const int N = dbinomD.nrow() - 1;
  const int G = dbinomD.ncol();
  const int n = D.size();
  if (w.size() != n || idx_last.size() != n || g_lo.size() != n || g_hi.size() != n) {
    stop("D, w, idx_last, g_lo and g_hi must have the same length");
  }
  if (theta.size() != G) stop("theta must have one element per column of dbinomD");
  // A tie group ends at or after its first member and the groups follow each other, so
  // that the loop over the groups in paired::tail_max() terminates
  for (int j = 0; j < n; ++j) {
    if (idx_last[j] < j || idx_last[j] >= n || (j > 0 && idx_last[j] < idx_last[j - 1])) {
      stop("idx_last does not describe consecutive tie groups");
    }
    if (D[j] < 0 || D[j] > N) stop("D must lie between 0 and the number of pairs");
    if (g_lo[j] < 0 || g_hi[j] >= G || g_lo[j] > g_hi[j]) stop("invalid grid range");
  }
  std::vector<double> dd(dbinomD.begin(), dbinomD.end());
  std::vector<int> Dv(D.begin(), D.end()), lv(idx_last.begin(), idx_last.end()),
    lo(g_lo.begin(), g_lo.end()), hi(g_hi.begin(), g_hi.end());
  std::vector<double> wv(w.begin(), w.end()), th(theta.begin(), theta.end());
  std::vector<double> out = paired::tail_max(dd, N, G, Dv, wv, lv, lo, hi, th, refine);
  return wrap(out);
}

// Probability of the rejection region R of N pairs for each pair (p10[i], p01[i])
// [[Rcpp::export]]
NumericVector paired_power(LogicalMatrix R, NumericVector p10, NumericVector p01) {
  const int N = R.nrow() - 1;
  if (N < 0 || R.ncol() != N + 1) stop("the rejection region must be square");
  if (p10.size() != p01.size()) stop("p10 and p01 must have the same length");
  const std::vector<double> lf = paired::log_factorials(N);
  NumericVector out(p10.size());
  for (R_xlen_t i = 0; i < p10.size(); ++i) {
    out[i] = paired::power_region(R.begin(), N, p10[i], p01[i], lf);
  }
  return out;
}

// Rejection probability of a re-estimation design for each pair (p10[i], p01[i]).
// rr_list holds the rejection regions of the distinct final numbers of pairs, and
// rr_id[D1] is the zero based position in rr_list of the region used after D1
// discordant pairs among the n1 interim pairs, D1 = 0, ..., n1.
// [[Rcpp::export]]
NumericVector paired_bssr_power(List rr_list, IntegerVector rr_id, int n1, NumericVector p10,
                                NumericVector p01) {
  if (rr_id.size() != n1 + 1) stop("rr_id must have n1 + 1 elements");
  if (p10.size() != p01.size()) stop("p10 and p01 must have the same length");
  std::vector<int> Ns;
  const std::vector<const int*> regions = region_pointers(rr_list, n1, Ns);
  for (int D1 = 0; D1 <= n1; ++D1) {
    if (rr_id[D1] < 0 || rr_id[D1] >= static_cast<int>(Ns.size())) {
      stop("rr_id must give positions in rr_list");
    }
  }
  const std::vector<double> lf = paired::log_factorials(max_pairs(Ns));
  NumericVector out(p10.size());
  for (R_xlen_t i = 0; i < p10.size(); ++i) {
    out[i] = paired::bssr_power(regions, Ns, rr_id.begin(), n1, p10[i], p01[i], lf);
  }
  return out;
}

// Conditional power of the final analysis with the region R, given D1[i] discordant
// pairs among n1 interim pairs: sum_D2 b(D2; N - n1, psi[i]) g_N(D1[i] + D2) with the
// directions binomial with probability pi_dir[i]. Elements with the same pi_dir should be
// adjacent, since g_N is computed once for each run of equal values
// [[Rcpp::export]]
NumericVector paired_cp(LogicalMatrix R, int n1, IntegerVector D1, NumericVector psi,
                        NumericVector pi_dir) {
  const int N = R.nrow() - 1;
  if (N < 0 || R.ncol() != N + 1) stop("the rejection region must be square");
  if (n1 < 0 || n1 > N) stop("n1 must lie between 0 and the number of pairs");
  const int n = D1.size();
  if (psi.size() != n || pi_dir.size() != n) {
    stop("D1, psi and pi_dir must have the same length");
  }
  for (int i = 0; i < n; ++i) {
    if (D1[i] < 0 || D1[i] > n1) stop("D1 must lie between 0 and n1");
    if (!(psi[i] >= 0.0 && psi[i] <= 1.0) || !(pi_dir[i] >= 0.0 && pi_dir[i] <= 1.0)) {
      stop("psi and pi_dir must lie in [0, 1]");
    }
  }
  const std::vector<int> d(D1.begin(), D1.end());
  const std::vector<double> p(psi.begin(), psi.end()), q(pi_dir.begin(), pi_dir.end());
  const std::vector<double> lf = paired::log_factorials(N);
  return wrap(paired::conditional_power(R.begin(), N, n1, d, p, q, lf));
}

// Bernstein coefficients of degree N of sum R(a, b) W[a + b] P(a, b) on a boundary of
// the null hypothesis, see paired::bern_coef()
// [[Rcpp::export]]
NumericVector paired_bern_coef(LogicalMatrix R, NumericVector W, int type, double m) {
  const int N = R.nrow() - 1;
  if (N < 0 || R.ncol() != N + 1) stop("the rejection region must be square");
  if (W.size() != N + 1) stop("W must have N + 1 elements");
  if (type < 0 || type > 2) stop("type must be 0, 1 or 2");
  if (type > 0 && !(m >= 0.0 && m < 1.0)) stop("m must lie in [0, 1)");
  const std::vector<double> lf = paired::log_factorials(N);
  return wrap(paired::bern_coef(R.begin(), W.begin(), N, type, m, lf));
}
