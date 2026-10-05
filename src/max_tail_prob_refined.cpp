#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

namespace {

// Value and first two derivatives of a polynomial in Bernstein form.
struct BernsteinValue {
  double f;
  double d1;
  double d2;
};

// Polynomial f(t) = sum_s h[s] b(s; N, t) at an interior point 0 < t < 1, where b is the
// binomial probability mass function. The terms are generated from the mode outwards by
// the ratio of successive binomial probabilities, which keeps each term positive and
// avoids evaluating the mass function N + 1 times. The derivatives follow from
// d b(s; N, t) / dt = b(s; N, t) (s - N t) / (t (1 - t)).
BernsteinValue bernstein_eval(const std::vector<double>& h, const std::vector<double>& up,
                              const std::vector<double>& down, int N, double t) {
  const double w = t * (1.0 - t);
  const double odds = t / (1.0 - t);
  int mode = static_cast<int>(std::floor((N + 1) * t));
  if (mode > N) mode = N;
  double s0 = 0.0, s1 = 0.0, s2 = 0.0;
  const double bm = R::dbinom(mode, N, t, 0);
  double b = bm;
  {
    const double u = mode - N * t;
    const double hb = h[mode] * b;
    s0 += hb;
    s1 += hb * u;
    s2 += hb * u * u;
  }
  for (int s = mode; s < N && b > 0.0; ++s) {
    b *= up[s] * odds;
    const double u = (s + 1) - N * t;
    const double hb = h[s + 1] * b;
    s0 += hb;
    s1 += hb * u;
    s2 += hb * u * u;
  }
  b = bm;
  for (int s = mode; s > 0 && b > 0.0; --s) {
    b *= down[s] / odds;
    const double u = (s - 1) - N * t;
    const double hb = h[s - 1] * b;
    s0 += hb;
    s1 += hb * u;
    s2 += hb * u * u;
  }
  BernsteinValue out;
  out.f = s0;
  out.d1 = s1 / w;
  out.d2 = (s2 - N * w * s0 - (1.0 - 2.0 * t) * s1) / (w * w);
  return out;
}

// Largest value of the polynomial found by a safeguarded Newton iteration for a zero of
// the derivative inside the bracket (lo, hi), started at t0. A Newton step is taken when
// the polynomial is concave at the current point and the step stays inside the bracket,
// and a bisection step otherwise. The bracket is shrunk by the sign of the derivative.
double local_max(const std::vector<double>& h, const std::vector<double>& up,
                 const std::vector<double>& down, int N, double lo, double hi, double t0) {
  double best = 0.0;
  double t = (t0 > lo && t0 < hi) ? t0 : 0.5 * (lo + hi);
  for (int it = 0; it < 200; ++it) {
    const BernsteinValue e = bernstein_eval(h, up, down, N, t);
    if (e.f > best) best = e.f;
    if (e.d1 > 0.0) {
      lo = t;
    } else {
      hi = t;
    }
    double next = 0.5 * (lo + hi);
    if (e.d2 < 0.0) {
      const double newton = t - e.d1 / e.d2;
      if (newton > lo && newton < hi) next = newton;
    }
    if (std::fabs(next - t) <= 1e-13 || hi - lo <= 1e-13) break;
    t = next;
  }
  return best;
}

} // namespace

// Tail probability maximized over the nuisance parameter with local refinement of the
// grid maximum. The arguments are those of max_tail_prob together with the grid itself.
// The grid values are accumulated in the same order as in max_tail_prob, so the grid
// maximum of every cell is reproduced exactly. In addition, the tail probability of each
// tie group is written as a polynomial in Bernstein form, sum_s h[s] b(s; N, theta), where
// h[s] is the conditional probability of the tail set given s responders in total. Every
// local maximum of the grid values within the range of a cell is then refined over the
// interval between its two neighbouring grid points, grid points closer than 1e-10 being
// passed over, and the result is the larger of the grid maximum and the refined values.
// Refinement is skipped when the grid maximum is within 1e-12 of one, since the tail
// probability cannot exceed one.
// [[Rcpp::export]]
NumericVector max_tail_prob_refined(NumericMatrix dbinom1, NumericMatrix dbinom2,
                                    IntegerVector x1, IntegerVector x2,
                                    IntegerVector idx_last,
                                    IntegerVector g_lo, IntegerVector g_hi,
                                    NumericVector theta) {
  const int n = x1.size();
  const int G = dbinom1.ncol();
  const int N1 = dbinom1.nrow() - 1;
  const int N2 = dbinom2.nrow() - 1;
  const int N = N1 + N2;
  if (dbinom2.ncol() != G || theta.size() != G) {
    stop("the grid does not match the probability matrices");
  }
  if (x2.size() != n || idx_last.size() != n || g_lo.size() != n || g_hi.size() != n) {
    stop("the cell vectors differ in length");
  }
  std::vector<double> up(N + 1, 0.0), down(N + 1, 0.0);
  for (int s = 0; s <= N; ++s) {
    up[s] = static_cast<double>(N - s) / (s + 1);
    down[s] = static_cast<double>(s) / (N - s + 1);
  }
  std::vector<double> v(G, 0.0);
  std::vector<double> h(N + 1, 0.0);
  NumericVector out(n);
  int start = 0;
  for (int m = 0; m < n; ++m) {
    if (x1[m] < 0 || x1[m] > N1 || x2[m] < 0 || x2[m] > N2 || idx_last[m] < m) {
      stop("invalid cell index");
    }
    for (int g = 0; g < G; ++g) {
      v[g] += dbinom1(x1[m], g) * dbinom2(x2[m], g);
    }
    const int s = x1[m] + x2[m];
    h[s] += R::dhyper(x1[m], N1, N2, s, 0);
    if (idx_last[m] != m) continue;
    // End of a tie group, which covers the cells start, ..., m
    int prev_lo = -1, prev_hi = -1;
    double prev_val = 0.0;
    for (int k = start; k <= m; ++k) {
      const int lo = g_lo[k], hi = g_hi[k];
      if (lo < 0 || hi >= G || lo > hi) stop("invalid grid range");
      if (lo == prev_lo && hi == prev_hi) {
        out[k] = prev_val;
        continue;
      }
      double best = 0.0;
      for (int g = lo; g <= hi; ++g) {
        if (v[g] > best) best = v[g];
      }
      if (best < 1.0 - 1e-12) {
        double refined = best;
        for (int g = lo; g <= hi; ++g) {
          const double left = (g > lo) ? v[g - 1] : -1.0;
          const double right = (g < hi) ? v[g + 1] : -1.0;
          if (v[g] > 0.0 && v[g] >= left && v[g] > right) {
            // The bracket ends at the nearest grid points that are not duplicates of
            // theta[g] up to rounding, so that it never collapses on one side
            int ga = g, gb = g;
            while (ga > lo && theta[g] - theta[ga] <= 1e-10) --ga;
            while (gb < hi && theta[gb] - theta[g] <= 1e-10) ++gb;
            const double a = theta[ga];
            const double b = theta[gb];
            if (b > a) {
              const double val = local_max(h, up, down, N, a, b, theta[g]);
              if (val > refined) refined = val;
            }
          }
        }
        best = refined;
      }
      out[k] = best;
      prev_lo = lo;
      prev_hi = hi;
      prev_val = best;
    }
    start = m + 1;
  }
  return out;
}
