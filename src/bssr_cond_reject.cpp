#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// Conditional rejection probabilities of a design with blinded sample size re-estimation
// under the null hypothesis of equal response probabilities.
//
// rr_list  rejection regions of the distinct final sample sizes, as logical matrices
// rr_id    zero based index into rr_list of the final sample size reached from each
//          pooled number of interim responders s = 0, ..., n11 + n12
// n21, n22 second-stage sample sizes of group 1 and group 2 for each element of rr_list
// n11, n12 interim sample sizes of group 1 and group 2
//
// Given the pooled number s of interim responders and the pooled number s2 of
// second-stage responders, the responder counts of group 1 in the two stages are
// independent hypergeometric counts, so the rejection probability is a double sum over
// them that does not involve the common response probability. For comparison, the
// rejection probability given only the total s + s2 is also returned. It treats the
// responder count of group 1 as one hypergeometric count among all final patients.
//
// The result lists every pair (s, s2), with s2 running over 0, ..., n21 + n22 of the
// final sample size reached from s, ordered by s and then by s2.
// [[Rcpp::export]]
List bssr_cond_reject(List rr_list, IntegerVector rr_id, IntegerVector n21,
                      IntegerVector n22, int n11, int n12) {
  const int K = rr_list.size();
  const int n1 = n11 + n12;
  if (n21.size() != K || n22.size() != K) {
    stop("n21 and n22 must have one element per rejection region");
  }
  if (rr_id.size() != n1 + 1) {
    stop("rr_id must have one element per pooled number of interim responders");
  }
  for (int s = 0; s <= n1; ++s) {
    if (rr_id[s] < 0 || rr_id[s] >= K) stop("rr_id lies outside rr_list");
  }
  // Rejection probability given the total number of responders, for each final size
  std::vector<std::vector<double> > given_total(K);
  for (int k = 0; k < K; ++k) {
    LogicalMatrix rr = rr_list[k];
    const int N1 = n11 + n21[k];
    const int N2 = n12 + n22[k];
    if (n21[k] < 0 || n22[k] < 0 || rr.nrow() != N1 + 1 || rr.ncol() != N2 + 1) {
      stop("the dimension of a rejection region does not match the final sample size");
    }
    given_total[k].assign(N1 + N2 + 1, 0.0);
    for (int t = 0; t <= N1 + N2; ++t) {
      double v = 0.0;
      for (int x1 = std::max(0, t - N2); x1 <= std::min(N1, t); ++x1) {
        if (rr(x1, t - x1) == TRUE) v += R::dhyper(x1, N1, N2, t, 0);
      }
      given_total[k][t] = v;
    }
  }
  int n_out = 0;
  for (int s = 0; s <= n1; ++s) n_out += n21[rr_id[s]] + n22[rr_id[s]] + 1;
  IntegerVector s_out(n_out), s2_out(n_out);
  NumericVector crp(n_out), crp_total(n_out);
  std::vector<double> h1(n11 + 1), h2;
  int idx = 0;
  for (int s = 0; s <= n1; ++s) {
    const int k = rr_id[s];
    LogicalMatrix rr = rr_list[k];
    const int m1 = n21[k];
    const int m2 = n22[k];
    // Interim responders of group 1 given s
    const int lo1 = std::max(0, s - n12);
    const int hi1 = std::min(n11, s);
    for (int a = lo1; a <= hi1; ++a) h1[a - lo1] = R::dhyper(a, n11, n12, s, 0);
    for (int s2 = 0; s2 <= m1 + m2; ++s2) {
      // Second-stage responders of group 1 given s2
      const int lo2 = std::max(0, s2 - m2);
      const int hi2 = std::min(m1, s2);
      h2.assign(hi2 - lo2 + 1, 0.0);
      for (int b = lo2; b <= hi2; ++b) h2[b - lo2] = R::dhyper(b, m1, m2, s2, 0);
      double c = 0.0;
      for (int a = lo1; a <= hi1; ++a) {
        double inner = 0.0;
        for (int b = lo2; b <= hi2; ++b) {
          if (rr(a + b, (s - a) + (s2 - b)) == TRUE) inner += h2[b - lo2];
        }
        c += h1[a - lo1] * inner;
      }
      s_out[idx] = s;
      s2_out[idx] = s2;
      crp[idx] = c;
      crp_total[idx] = given_total[k][s + s2];
      ++idx;
    }
  }
  return List::create(_["s"] = s_out, _["s2"] = s2_out, _["crp"] = crp,
                      _["crp.total"] = crp_total);
}
