#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// Null tail probability of an ordering statistic, maximized over a grid of the nuisance
// parameter. The cells are supplied already sorted from the most extreme to the least
// extreme. For each grid point the running sum of the null cell probabilities is formed
// along that ordering, and the value taken by each cell is the running sum evaluated at
// the last member of its tie group.
//
// dbinom1  (N1 + 1) by G matrix of null probabilities of the group 1 responder count
// dbinom2  (N2 + 1) by G matrix of null probabilities of the group 2 responder count
// x1, x2   zero based responder counts of the sorted cells
// idx_last zero based position of the last member of the tie group of each sorted cell
// g_lo     zero based first grid index over which each cell is maximized
// g_hi     zero based last grid index over which each cell is maximized
// [[Rcpp::export]]
NumericVector max_tail_prob(NumericMatrix dbinom1, NumericMatrix dbinom2,
                            IntegerVector x1, IntegerVector x2,
                            IntegerVector idx_last,
                            IntegerVector g_lo, IntegerVector g_hi) {
  const int n = x1.size();
  const int G = dbinom1.ncol();
  NumericVector out(n);
  std::vector<double> cum(n);
  for (int g = 0; g < G; ++g) {
    double run = 0.0;
    for (int m = 0; m < n; ++m) {
      run += dbinom1(x1[m], g) * dbinom2(x2[m], g);
      cum[m] = run;
    }
    for (int k = 0; k < n; ++k) {
      if (g >= g_lo[k] && g <= g_hi[k]) {
        const double v = cum[idx_last[k]];
        if (v > out[k]) out[k] = v;
      }
    }
  }
  return out;
}

// Power of a test defined by a rejection region, that is the sum over the rejection
// region of the product of the two independent binomial probability mass functions.
// [[Rcpp::export]]
double power_from_rr(LogicalMatrix RR, NumericVector dbinom1, NumericVector dbinom2) {
  const int n1 = RR.nrow();
  const int n2 = RR.ncol();
  double total = 0.0;
  for (int i = 0; i < n1; ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < n2; ++j) {
      if (RR(i, j)) row_sum += dbinom2[j];
    }
    total += dbinom1[i] * row_sum;
  }
  return total;
}
