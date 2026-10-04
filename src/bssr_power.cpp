#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// Rejected cells of a rejection region, stored column by column as runs of consecutive
// group 1 responder counts. The runs of column j are lo[q] to hi[q] for q running from
// ptr[j] to ptr[j + 1] - 1. A one-sided region has at most one run per column and a
// two-sided region at most two, so a conditional probability over a column reduces to a
// few differences of a binomial distribution function.
struct ColumnRuns {
  std::vector<int> ptr;
  std::vector<int> lo;
  std::vector<int> hi;
};

static ColumnRuns column_runs(LogicalMatrix rr) {
  const int n1 = rr.nrow();
  const int n2 = rr.ncol();
  ColumnRuns out;
  out.ptr.assign(n2 + 1, 0);
  for (int j = 0; j < n2; ++j) {
    out.ptr[j] = static_cast<int>(out.lo.size());
    int i = 0;
    while (i < n1) {
      if (rr(i, j) == TRUE) {
        const int start = i;
        while (i + 1 < n1 && rr(i + 1, j) == TRUE) ++i;
        out.lo.push_back(start);
        out.hi.push_back(i);
      }
      ++i;
    }
  }
  out.ptr[n2] = static_cast<int>(out.lo.size());
  return out;
}

// Rejection probability of a design with blinded sample size re-estimation.
//
// rr_list  rejection regions of the distinct final sample sizes, as logical matrices
// rr_id    zero based index into rr_list of the final sample size of each interim cell
// x11, x12 interim responder counts of group 1 and group 2 of each interim cell
// n21, n22 second-stage sample sizes of group 1 and group 2 for each element of rr_list
// p1, p2   response probabilities of group 1 and group 2, one element per scenario
// n11, n12 interim sample sizes of group 1 and group 2
//
// For each scenario the conditional rejection probability of every interim cell is the
// sum, over the second-stage responder count of group 2, of the probability that the
// final responder count of group 1 falls in a rejected run of the corresponding column.
// The result is averaged over the binomial distribution of the interim cells.
// [[Rcpp::export]]
NumericVector bssr_power(List rr_list, IntegerVector rr_id, IntegerVector x11,
                         IntegerVector x12, IntegerVector n21, IntegerVector n22,
                         NumericVector p1, NumericVector p2, int n11, int n12) {
  const int K = rr_list.size();
  if (n21.size() != K || n22.size() != K) {
    stop("n21 and n22 must have one element per rejection region");
  }
  std::vector<ColumnRuns> runs(K);
  for (int k = 0; k < K; ++k) {
    LogicalMatrix rr = rr_list[k];
    if (rr.nrow() != n11 + n21[k] + 1 || rr.ncol() != n12 + n22[k] + 1) {
      stop("the dimension of a rejection region does not match the final sample size");
    }
    runs[k] = column_runs(rr);
  }
  const int n_cell = x11.size();
  const int n_scen = p1.size();
  NumericVector out(n_scen);
  std::vector<std::vector<double> > pmf2(K), lower1(K), upper1(K);
  for (int s = 0; s < n_scen; ++s) {
    // Second-stage distributions for each final sample size
    for (int k = 0; k < K; ++k) {
      const int m21 = n21[k];
      const int m22 = n22[k];
      pmf2[k].assign(m22 + 1, 0.0);
      for (int b = 0; b <= m22; ++b) pmf2[k][b] = R::dbinom(b, m22, p2[s], 0);
      std::vector<double> pmf1(m21 + 1);
      for (int a = 0; a <= m21; ++a) pmf1[a] = R::dbinom(a, m21, p1[s], 0);
      // Lower and upper cumulative sums are both kept, so that a tail probability is
      // never obtained as one minus a probability close to one
      lower1[k].assign(m21 + 1, 0.0);
      upper1[k].assign(m21 + 1, 0.0);
      double acc = 0.0;
      for (int a = 0; a <= m21; ++a) {
        acc += pmf1[a];
        lower1[k][a] = acc;
      }
      acc = 0.0;
      for (int a = m21; a >= 0; --a) {
        acc += pmf1[a];
        upper1[k][a] = acc;
      }
    }
    double total = 0.0;
    for (int c = 0; c < n_cell; ++c) {
      const int k = rr_id[c];
      const ColumnRuns& cr = runs[k];
      const int a0 = x11[c];
      const int b0 = x12[c];
      const int m21 = n21[k];
      const int m22 = n22[k];
      const std::vector<double>& lo1 = lower1[k];
      const std::vector<double>& up1 = upper1[k];
      double cp = 0.0;
      for (int b = 0; b <= m22; ++b) {
        const int col = b0 + b;
        double col_prob = 0.0;
        for (int q = cr.ptr[col]; q < cr.ptr[col + 1]; ++q) {
          // Second-stage counts of group 1 that bring the final count into the run
          int lo = cr.lo[q] - a0;
          int hi = cr.hi[q] - a0;
          if (lo < 0) lo = 0;
          if (hi > m21) hi = m21;
          if (lo > hi) continue;
          if (lo == 0) {
            col_prob += lo1[hi];
          } else if (hi == m21) {
            col_prob += up1[lo];
          } else {
            col_prob += lo1[hi] - lo1[lo - 1];
          }
        }
        cp += pmf2[k][b] * col_prob;
      }
      const double v = R::dbinom(a0, n11, p1[s], 0) * R::dbinom(b0, n12, p2[s], 0) * cp;
      if (!ISNAN(v)) total += v;
    }
    out[s] = total;
  }
  return out;
}
