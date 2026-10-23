// Computations for the paired (matched-pairs) binary design. Plain C++ without Rcpp, so
// that the routines can also be compiled and checked outside R; the Rcpp interfaces are
// in paired.cpp.
//
// A trial of N pairs is summarized by the numbers a = x10 of pairs in which only
// treatment 1 responds and b = x01 of pairs in which only treatment 2 responds, with
// a + b <= N. A rejection region is an (N + 1) by (N + 1) matrix of 0 and 1 in column
// major order, indexed by (a, b), with 0 in the cells a + b > N. Given the number
// D = a + b of discordant pairs, a is binomial with D trials and success probability
// pi = p10 / (p10 + p01), and D is binomial with N trials and success probability
// psi = p10 + p01.
#ifndef BBSSR_PAIRED_CORE_H
#define BBSSR_PAIRED_CORE_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace paired {

// lf[k] = log(k!) for k = 0, ..., n
inline std::vector<double> log_factorials(int n) {
  std::vector<double> lf(n + 1);
  for (int k = 0; k <= n; ++k) lf[k] = std::lgamma(k + 1.0);
  return lf;
}

// k log(p), with 0 log(0) = 0
inline double xlogy(int k, double p) {
  if (k == 0) return 0.0;
  if (p <= 0.0) return -std::numeric_limits<double>::infinity();
  return k * std::log(p);
}

// Binomial probability b(x; n, p)
inline double binom_pmf(int x, int n, double p, const std::vector<double>& lf) {
  if (x < 0 || x > n) return 0.0;
  if (p <= 0.0) return x == 0 ? 1.0 : 0.0;
  if (p >= 1.0) return x == n ? 1.0 : 0.0;
  return std::exp(lf[n] - lf[x] - lf[n - x] + x * std::log(p) + (n - x) * std::log1p(-p));
}

// b(x; n, p) for x = 0, ..., n, written to out. The probability at the mode is computed
// from the log factorials and the others by the ratios of consecutive probabilities,
// moving away from the mode, so that each costs one multiplication. Probabilities far
// from the mode underflow to 0
inline void binom_row(int n, double p, const std::vector<double>& lf,
                      std::vector<double>& out) {
  out.assign(n + 1, 0.0);
  if (p <= 0.0) {
    out[0] = 1.0;
    return;
  }
  if (p >= 1.0) {
    out[n] = 1.0;
    return;
  }
  int m = static_cast<int>(std::floor((n + 1) * p));
  if (m > n) m = n;
  if (m < 0) m = 0;
  out[m] = std::exp(lf[n] - lf[m] - lf[n - m] + m * std::log(p) + (n - m) * std::log1p(-p));
  const double odds = p / (1.0 - p);
  for (int x = m + 1; x <= n; ++x) out[x] = out[x - 1] * (n - x + 1) / x * odds;
  for (int x = m - 1; x >= 0; --x) out[x] = out[x + 1] * (x + 1) / (n - x) / odds;
}

// b(x; n, p) for x = 0, ..., n
inline std::vector<double> binom_vector(int n, double p, const std::vector<double>& lf) {
  std::vector<double> out;
  binom_row(n, p, lf, out);
  return out;
}

// g[D] = sum over a of b(a; D, pi) R(a, D - a), the rejection probability given D
// discordant pairs, for D = d_lo, ..., d_hi (the other entries are 0)
inline std::vector<double> g_values(const int* R, int N, double pi, int d_lo, int d_hi,
                                    const std::vector<double>& lf) {
  std::vector<double> g(N + 1, 0.0), row;
  if (d_lo < 0) d_lo = 0;
  if (d_hi > N) d_hi = N;
  for (int D = d_lo; D <= d_hi; ++D) {
    binom_row(D, pi, lf, row);
    double s = 0.0;
    for (int a = 0; a <= D; ++a) {
      if (R[a + (D - a) * (N + 1)]) s += row[a];
    }
    g[D] = s;
  }
  return g;
}

// Probability of the rejection region R of N pairs, a sum of trinomial probabilities
inline double power_region(const int* R, int N, double p10, double p01,
                           const std::vector<double>& lf) {
  double pc = 1.0 - p10 - p01;
  if (pc < 0.0) pc = 0.0;
  double total = 0.0;
  for (int b = 0; b <= N; ++b) {
    for (int a = 0; a + b <= N; ++a) {
      if (!R[a + b * (N + 1)]) continue;
      const int c = N - a - b;
      const double l = lf[N] - lf[a] - lf[b] - lf[c] + xlogy(a, p10) + xlogy(b, p01) +
        xlogy(c, pc);
      if (std::isfinite(l)) total += std::exp(l);
    }
  }
  return total;
}

// Rejection probability of a design with blinded sample size re-estimation. The interim
// analysis has n1 pairs, the final number of pairs after D1 discordant interim pairs is
// Ns[id[D1]], and the final analysis uses the region regions[id[D1]]. Since the
// re-estimation uses D1 only, and given the numbers of discordant pairs the directions
// are binomial with the same pi in both stages,
// P = sum_D1 b(D1; n1, psi) sum_D2 b(D2; N - n1, psi) g_N(D1 + D2).
inline double bssr_power(const std::vector<const int*>& regions, const std::vector<int>& Ns,
                         const int* id, int n1, double p10, double p01,
                         const std::vector<double>& lf) {
  const double psi = p10 + p01;
  if (psi <= 0.0) {
    // Every pair is concordant: D1 = D2 = 0
    return regions[id[0]][0] ? 1.0 : 0.0;
  }
  const double pi = p10 / psi;
  const std::vector<double> b1 = binom_vector(n1, psi, lf);
  double total = 0.0;
  std::vector<int> done(regions.size(), 0);
  std::vector<std::vector<double> > g(regions.size()), b2(regions.size());
  for (int D1 = 0; D1 <= n1; ++D1) {
    const int k = id[D1];
    const int N = Ns[k];
    const int n2 = N - n1;
    if (!done[k]) {
      g[k] = g_values(regions[k], N, pi, 0, N, lf);
      b2[k] = binom_vector(n2, psi, lf);
      done[k] = 1;
    }
    double inner = 0.0;
    for (int D2 = 0; D2 <= n2; ++D2) inner += b2[k][D2] * g[k][D1 + D2];
    total += b1[D1] * inner;
  }
  return total;
}

// Conditional power of the final analysis with the region R of N pairs, given D1[i]
// discordant pairs among n1 interim pairs, under the discordance probability psi[i] and
// the direction probability pi[i]: sum_D2 b(D2; N - n1, psi[i]) g_N(D1[i] + D2; pi[i]).
// The values g_N(D; pi) are computed once for each distinct pi, over all D.
inline std::vector<double> conditional_power(const int* R, int N, int n1,
                                             const std::vector<int>& D1,
                                             const std::vector<double>& psi,
                                             const std::vector<double>& pi,
                                             const std::vector<double>& lf) {
  const int n2 = N - n1;
  std::vector<double> out(D1.size(), 0.0), b2, g;
  double pi_done = std::numeric_limits<double>::quiet_NaN();
  for (std::size_t i = 0; i < D1.size(); ++i) {
    if (!(pi[i] == pi_done)) {
      g = g_values(R, N, pi[i], 0, N, lf);
      pi_done = pi[i];
    }
    binom_row(n2, psi[i], lf, b2);
    double s = 0.0;
    for (int D2 = 0; D2 <= n2; ++D2) s += b2[D2] * g[D1[i] + D2];
    out[i] = s;
  }
  return out;
}

// Value at psi of sum_D h[D] b(D; N, psi)
inline double tail_value(const std::vector<double>& h, int N, double psi,
                         const std::vector<double>& lf) {
  double s = 0.0;
  for (int D = 0; D <= N; ++D) if (h[D] != 0.0) s += h[D] * binom_pmf(D, N, psi, lf);
  return s;
}

// Largest value over [lo, hi] found by a golden-section search for a local maximum
inline double golden_max(const std::vector<double>& h, int N, double lo, double hi,
                         const std::vector<double>& lf) {
  const double r = 0.5 * (std::sqrt(5.0) - 1.0);
  double best = std::max(tail_value(h, N, lo, lf), tail_value(h, N, hi, lf));
  double x1 = hi - r * (hi - lo), x2 = lo + r * (hi - lo);
  double f1 = tail_value(h, N, x1, lf), f2 = tail_value(h, N, x2, lf);
  for (int it = 0; it < 200 && hi - lo > 1e-13; ++it) {
    if (f1 >= f2) {
      hi = x2; x2 = x1; f2 = f1;
      x1 = hi - r * (hi - lo); f1 = tail_value(h, N, x1, lf);
    } else {
      lo = x1; x1 = x2; f1 = f2;
      x2 = lo + r * (hi - lo); f2 = tail_value(h, N, x2, lf);
    }
    if (f1 > best) best = f1;
    if (f2 > best) best = f2;
  }
  return best;
}

// Null tail probabilities of an ordering of the cells, maximized over the discordance
// probability psi. The n cells are sorted from the most extreme to the least extreme;
// cell j has D[j] discordant pairs and the null weight w[j] = b(a; D, 1/2), so that its
// null probability is b(D; N, psi) w[j]. last[j] is the position of the last member of
// the tie group of cell j, and the tail probability of cell j is the running sum up to
// last[j]. dbinomD[D + g * (N + 1)] = b(D; N, theta[g]). The maximum over the grid points
// g_lo[j], ..., g_hi[j] is taken, and with refine every local maximum on these points is
// refined by a golden-section search between its neighbours.
inline std::vector<double> tail_max(const std::vector<double>& dbinomD, int N, int G,
                                    const std::vector<int>& D, const std::vector<double>& w,
                                    const std::vector<int>& last,
                                    const std::vector<int>& g_lo,
                                    const std::vector<int>& g_hi,
                                    const std::vector<double>& theta, bool refine) {
  const int n = static_cast<int>(D.size());
  const std::vector<double> lf = log_factorials(N);
  std::vector<double> out(n, 0.0), prof(G, 0.0), h(N + 1, 0.0);
  int m = 0;
  while (m < n) {
    const int end = last[m];
    for (int j = m; j <= end; ++j) {
      for (int g = 0; g < G; ++g) prof[g] += dbinomD[D[j] + g * (N + 1)] * w[j];
      h[D[j]] += w[j];
    }
    int prev_lo = -1, prev_hi = -1;
    double prev = 0.0;
    for (int k = m; k <= end; ++k) {
      if (g_lo[k] == prev_lo && g_hi[k] == prev_hi) {
        out[k] = prev;
        continue;
      }
      double best = 0.0;
      for (int g = g_lo[k]; g <= g_hi[k]; ++g) if (prof[g] > best) best = prof[g];
      if (refine) {
        for (int g = g_lo[k]; g <= g_hi[k]; ++g) {
          const bool left = (g == g_lo[k]) || prof[g] >= prof[g - 1];
          const bool right = (g == g_hi[k]) || prof[g] >= prof[g + 1];
          if (!(left && right)) continue;
          const double lo = theta[g > g_lo[k] ? g - 1 : g];
          const double hi = theta[g < g_hi[k] ? g + 1 : g];
          if (hi > lo) {
            const double v = golden_max(h, N, lo, hi, lf);
            if (v > best) best = v;
          }
        }
      }
      out[k] = best;
      prev_lo = g_lo[k];
      prev_hi = g_hi[k];
      prev = best;
    }
    m = end + 1;
  }
  return out;
}

// Bernstein coefficients, of degree N in t on [0, 1], of
// sum over (a, b) of R(a, b) W[a + b] P(a, b) on a boundary of the null hypothesis:
//   type 0: p10 = p01 = t / 2 (psi = t)
//   type 1: p10 - p01 = -m, psi = m + (1 - m) t, so that p10 = (1 - m) t / 2,
//           p01 = m (1 - t) + (1 + m) t / 2 and 1 - psi = (1 - m) (1 - t)
//   type 2: p10 - p01 = m, the same with p10 and p01 exchanged
inline std::vector<double> bern_coef(const int* R, const double* W, int N, int type, double m,
                                     const std::vector<double>& lf) {
  std::vector<double> beta(N + 1, 0.0);
  if (type == 0 || m <= 0.0) {
    for (int D = 0; D <= N; ++D) {
      if (W[D] == 0.0) continue;
      double s = 0.0;
      for (int a = 0; a <= D; ++a) {
        if (R[a + (D - a) * (N + 1)]) {
          s += std::exp(lf[D] - lf[a] - lf[D - a] - D * std::log(2.0));
        }
      }
      beta[D] = W[D] * s;
    }
    return beta;
  }
  const double lh = std::log((1.0 - m) / 2.0), lq = std::log1p(-m), lm = std::log(m),
    lfp = std::log((1.0 + m) / 2.0);
  for (int b = 0; b <= N; ++b) {
    for (int a = 0; a + b <= N; ++a) {
      if (!R[a + b * (N + 1)] || W[a + b] == 0.0) continue;
      const int c = N - a - b;
      // u is the count whose probability is proportional to t, v the other discordant
      // count, whose probability is m (1 - t) + (1 + m) t / 2
      const int u = (type == 1) ? a : b;
      const int v = (type == 1) ? b : a;
      const double base = std::log(W[a + b]) + u * lh + c * lq - lf[u] - lf[c];
      for (int k = 0; k <= v; ++k) {
        const int j = u + k;
        beta[j] += std::exp(base + lf[j] + lf[N - j] - lf[k] - lf[v - k] + (v - k) * lm +
                            k * lfp);
      }
    }
  }
  return beta;
}

}  // namespace paired

#endif
