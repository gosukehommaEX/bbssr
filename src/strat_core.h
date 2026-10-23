// Computations for stratified designs with a binary endpoint. Plain C++ without Rcpp, so
// that the routines can also be compiled and checked outside R; the Rcpp interfaces are
// in strat.cpp.
//
// A design has J strata, and stratum j has N1[j] patients in group 1 and N2[j] in group
// 2. The outcome of stratum j is summarized by the number x of responders of group 1 and
// the number z of responders of the stratum. All the tests use the statistic
// S = sum_j x_j, and given z = (z_1, ..., z_J) they reject when S is at most a lower or at
// least an upper critical value. A vector of stratum totals z is stored at the position
// z_1 + (N_1 + 1) (z_2 + (N_2 + 1) (z_3 + ...)), where N_j = N1[j] + N2[j], the first
// stratum varying fastest.
#ifndef BBSSR_STRAT_CORE_H
#define BBSSR_STRAT_CORE_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace strat {

// Rejection tolerance of the package: a p-value is rejected when it is below the level
// by more than this amount, the default tolerance of fpCompare
const double REJECT_TOL = 1.4901161193847656e-08;
// Relative tolerance of tie_groups(), which groups equal keys of the two-sided p-values
const double TIE_TOL = 1e-10;

enum TestType { CMH = 0, CMH_CC = 1, COCHRAN = 2, COCHRAN_CC = 3, FISHER = 4 };
enum Alternative { GREATER = 0, LESS = 1, TWO_SIDED = 2 };
enum TsMethod { MINLIKE = 0, CENTRAL = 1, BLAKER = 2 };

struct Test {
  int type;
  int alternative;
  int tsmethod;
  double alpha;
};

struct Design {
  std::vector<int> N1;
  std::vector<int> N2;
  int J() const { return static_cast<int>(N1.size()); }
  int total(int j) const { return N1[j] + N2[j]; }
  // Number of vectors of stratum totals
  std::size_t cells() const {
    std::size_t n = 1;
    for (int j = 0; j < J(); ++j) n *= static_cast<std::size_t>(total(j) + 1);
    return n;
  }
};

// Critical values for one vector of stratum totals: reject when S <= lo or S >= hi.
// lo = Smin - 1 and hi = Smax + 1 when the corresponding tail is empty
struct Crit {
  int lo;
  int hi;
};

// lf[k] = log(k!) for k = 0, ..., n
inline std::vector<double> log_factorials(int n) {
  std::vector<double> lf(n + 1);
  for (int k = 0; k <= n; ++k) lf[k] = std::lgamma(k + 1.0);
  return lf;
}

// A distribution on the integers lo, ..., lo + p.size() - 1
struct Dist {
  int lo;
  std::vector<double> p;
  int hi() const { return lo + static_cast<int>(p.size()) - 1; }
};

// b(x; n, p) for x = 0, ..., n. The probability at the mode is computed from the log
// factorials and the others by the ratios of consecutive probabilities, moving away from
// the mode. Probabilities far from the mode underflow to 0
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

// Hypergeometric distribution of the number x of group 1 responders given z responders
// among m1 + m2 patients: choose(m1, x) choose(m2, z - x) / choose(m1 + m2, z) for
// x = max(0, z - m2), ..., min(z, m1). Computed from the mode by the ratios of
// consecutive probabilities
inline Dist hyper(int m1, int m2, int z, const std::vector<double>& lf) {
  Dist d;
  const int lo = std::max(0, z - m2);
  const int hi = std::min(z, m1);
  d.lo = lo;
  d.p.assign(hi - lo + 1, 0.0);
  // Mode of the hypergeometric distribution
  int mode = static_cast<int>(std::floor((z + 1.0) * (m1 + 1.0) / (m1 + m2 + 2.0)));
  if (mode < lo) mode = lo;
  if (mode > hi) mode = hi;
  const int n = m1 + m2;
  d.p[mode - lo] = std::exp(lf[m1] - lf[mode] - lf[m1 - mode] + lf[m2] - lf[z - mode] -
                            lf[m2 - z + mode] - (lf[n] - lf[z] - lf[n - z]));
  for (int x = mode + 1; x <= hi; ++x) {
    d.p[x - lo] = d.p[x - 1 - lo] * (m1 - x + 1.0) * (z - x + 1.0) /
      (x * (m2 - z + x + 0.0));
  }
  for (int x = mode - 1; x >= lo; --x) {
    d.p[x - lo] = d.p[x + 1 - lo] * (x + 1.0) * (m2 - z + x + 1.0) /
      ((m1 - x + 0.0) * (z - x + 0.0));
  }
  return d;
}

inline Dist convolve(const Dist& a, const Dist& b) {
  Dist c;
  c.lo = a.lo + b.lo;
  c.p.assign(a.p.size() + b.p.size() - 1, 0.0);
  for (std::size_t i = 0; i < a.p.size(); ++i) {
    const double ai = a.p[i];
    if (ai == 0.0) continue;
    double* out = &c.p[i];
    for (std::size_t k = 0; k < b.p.size(); ++k) out[k] += ai * b.p[k];
  }
  return c;
}

// Standard normal upper tail probability
inline double norm_upper(double x) { return 0.5 * std::erfc(x / std::sqrt(2.0)); }

inline bool rejects(double p, double alpha) { return alpha - p > REJECT_TOL; }

// Mean E and variance V of S given the stratum totals under the null hypothesis, with the
// hypergeometric variance (CMH) or the binomial variance (Cochran)
inline void null_moments(const Design& d, const int* z, bool cochran, double& E, double& V) {
  E = 0.0;
  V = 0.0;
  for (int j = 0; j < d.J(); ++j) {
    const double m1 = d.N1[j], m2 = d.N2[j], n = m1 + m2, zj = z[j];
    E += m1 * zj / n;
    if (n > 1.0) {
      V += cochran ? m1 * m2 * zj * (n - zj) / (n * n * n)
                   : m1 * m2 * zj * (n - zj) / (n * n * (n - 1.0));
    }
  }
}

// p-value of an asymptotic test for the statistic S, given E and V. With the continuity
// correction the deviation |S - E| is reduced by 0.5 when it is at least 0.5, as in
// stats::mantelhaen.test(). A zero variance gives the p-value 1
inline double asymptotic_pvalue(double S, double E, double V, bool cc, int alternative) {
  if (!(V > 0.0)) return 1.0;
  const double delta = S - E;
  const double yates = (cc && std::fabs(delta) >= 0.5) ? 0.5 : 0.0;
  const double sd = std::sqrt(V);
  if (alternative == TWO_SIDED) {
    return std::min(1.0, 2.0 * norm_upper((std::fabs(delta) - yates) / sd));
  }
  const double sgn = delta > 0.0 ? 1.0 : (delta < 0.0 ? -1.0 : 0.0);
  const double Z = sgn * (std::fabs(delta) - yates) / sd;
  return alternative == GREATER ? norm_upper(Z) : norm_upper(-Z);
}

// Null distribution of S given the stratum totals: the convolution of the hypergeometric
// distributions of the strata
inline Dist null_dist(const Design& d, const int* z, const std::vector<double>& lf) {
  Dist f = hyper(d.N1[0], d.N2[0], z[0], lf);
  for (int j = 1; j < d.J(); ++j) f = convolve(f, hyper(d.N1[j], d.N2[j], z[j], lf));
  return f;
}

// p-values of the stratified exact conditional test for every S of the null distribution
// f. The two-sided p-values of minlike and blaker sum the probabilities of the values whose
// key is at most that of S, with the ties of tie_groups()
inline std::vector<double> fisher_pvalues(const Dist& f, int alternative, int tsmethod) {
  const std::size_t n = f.p.size();
  std::vector<double> up(n), lo(n), pv(n);
  double acc = 0.0;
  for (std::size_t k = n; k-- > 0;) {
    acc += f.p[k];
    up[k] = acc;
  }
  acc = 0.0;
  for (std::size_t k = 0; k < n; ++k) {
    acc += f.p[k];
    lo[k] = acc;
  }
  if (alternative == GREATER) return up;
  if (alternative == LESS) return lo;
  if (tsmethod == CENTRAL) {
    for (std::size_t k = 0; k < n; ++k) pv[k] = std::min(1.0, 2.0 * std::min(lo[k], up[k]));
    return pv;
  }
  std::vector<double> key(n);
  for (std::size_t k = 0; k < n; ++k) {
    key[k] = tsmethod == BLAKER ? std::min(lo[k], up[k]) : f.p[k];
  }
  std::vector<std::size_t> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(),
                   [&key](std::size_t a, std::size_t b) { return key[a] < key[b]; });
  std::vector<double> cum(n);
  acc = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    acc += f.p[ord[i]];
    cum[i] = acc;
  }
  // Groups of adjacent sorted keys within the relative tolerance; every member takes the
  // cumulative probability at the end of its group
  std::size_t start = 0;
  while (start < n) {
    std::size_t end = start;
    while (end + 1 < n) {
      const double a = key[ord[end]], b = key[ord[end + 1]];
      const double ref = std::max(std::fabs(a), std::fabs(b));
      if (std::fabs(b - a) > TIE_TOL * ref) break;
      ++end;
    }
    for (std::size_t i = start; i <= end; ++i) pv[ord[i]] = std::min(1.0, cum[end]);
    start = end + 1;
  }
  return pv;
}

// Critical values from the rejection indicators of S = Smin, ..., Smax. The rejected
// values must form a lower and an upper tail; otherwise an error is raised
inline Crit tails_of(const std::vector<char>& rej, int Smin) {
  const int n = static_cast<int>(rej.size());
  int a = 0;
  while (a < n && rej[a]) ++a;
  int b = n - 1;
  while (b >= a && rej[b]) --b;
  for (int k = a; k <= b; ++k) {
    if (rej[k]) throw std::runtime_error("the rejection region is not formed by two tails");
  }
  Crit c;
  c.lo = Smin + a - 1;
  c.hi = Smin + b + 1;
  if (a == n) c.hi = Smin + n;  // every value rejected: the lower tail holds them all
  return c;
}

// Critical values of the test for the stratum totals z
inline Crit critical(const Design& d, const Test& t, const int* z,
                     const std::vector<double>& lf) {
  int Smin = 0, Smax = 0;
  for (int j = 0; j < d.J(); ++j) {
    Smin += std::max(0, z[j] - d.N2[j]);
    Smax += std::min(z[j], d.N1[j]);
  }
  Crit c;
  c.lo = Smin - 1;
  c.hi = Smax + 1;
  if (t.type == FISHER) {
    const Dist f = null_dist(d, z, lf);
    const std::vector<double> pv = fisher_pvalues(f, t.alternative, t.tsmethod);
    std::vector<char> rej(pv.size());
    for (std::size_t k = 0; k < pv.size(); ++k) rej[k] = rejects(pv[k], t.alpha);
    return tails_of(rej, Smin);
  }
  const bool cochran = t.type == COCHRAN || t.type == COCHRAN_CC;
  const bool cc = t.type == CMH_CC || t.type == COCHRAN_CC;
  double E, V;
  null_moments(d, z, cochran, E, V);
  if (!(V > 0.0)) return c;
  auto rej = [&](int S) {
    return rejects(asymptotic_pvalue(S, E, V, cc, t.alternative), t.alpha);
  };
  // Values above E can only be rejected in favour of 'greater' and values below E in
  // favour of 'less'. Above E the rejected values are closed upwards and below E
  // downwards, also with the continuity correction of mantelhaen.test(), so each boundary
  // is found by bisection
  const int Su = std::max(Smin, static_cast<int>(std::floor(E)) + 1);
  const int Sl = std::min(Smax, static_cast<int>(std::ceil(E)) - 1);
  if (t.alternative != LESS && Su <= Smax && rej(Smax)) {
    int b = Smax;  // rejected
    if (rej(Su)) {
      b = Su;
    } else {
      int a = Su;  // not rejected
      while (b - a > 1) {
        const int mid = a + (b - a) / 2;
        if (rej(mid)) b = mid; else a = mid;
      }
    }
    c.hi = b;
  }
  if (t.alternative != GREATER && Sl >= Smin && rej(Smin)) {
    int a = Smin;  // rejected
    if (rej(Sl)) {
      a = Sl;
    } else {
      int b = Sl;  // not rejected
      while (b - a > 1) {
        const int mid = a + (b - a) / 2;
        if (rej(mid)) a = mid; else b = mid;
      }
    }
    c.lo = a;
  }
  return c;
}

// Critical values for every vector of stratum totals of the design, in the storage order
inline void critical_grid(const Design& d, const Test& t, const std::vector<double>& lf,
                          std::vector<int>& lo, std::vector<int>& hi) {
  const std::size_t n = d.cells();
  lo.assign(n, 0);
  hi.assign(n, 0);
  std::vector<int> z(d.J(), 0);
  for (std::size_t k = 0; k < n; ++k) {
    const Crit c = critical(d, t, z.data(), lf);
    lo[k] = c.lo;
    hi[k] = c.hi;
    for (int j = 0; j < d.J(); ++j) {
      if (++z[j] <= d.total(j)) break;
      z[j] = 0;
    }
  }
}

// Largest p-value rejected at the level t.alpha (pr) and smallest p-value not rejected
// (pn) over all outcomes of the design, -1 and 2 when there is none, from which
// report_level() chooses a level with the same rejection regions
inline void level_bracket(const Design& d, const Test& t, const std::vector<double>& lf,
                          double& pr, double& pn) {
  pr = -1.0;
  pn = 2.0;
  const bool cochran = t.type == COCHRAN || t.type == COCHRAN_CC;
  const bool cc = t.type == CMH_CC || t.type == COCHRAN_CC;
  std::vector<int> z(d.J(), 0);
  const std::size_t n = d.cells();
  for (std::size_t k = 0; k < n; ++k) {
    int Smin = 0, Smax = 0;
    for (int j = 0; j < d.J(); ++j) {
      Smin += std::max(0, z[j] - d.N2[j]);
      Smax += std::min(z[j], d.N1[j]);
    }
    std::vector<double> pv;
    if (t.type == FISHER) {
      pv = fisher_pvalues(null_dist(d, z.data(), lf), t.alternative, t.tsmethod);
    } else {
      double E, V;
      null_moments(d, z.data(), cochran, E, V);
      pv.resize(Smax - Smin + 1);
      for (int S = Smin; S <= Smax; ++S) {
        pv[S - Smin] = asymptotic_pvalue(S, E, V, cc, t.alternative);
      }
    }
    for (double p : pv) {
      if (rejects(p, t.alpha)) {
        pr = std::max(pr, p);
      } else {
        pn = std::min(pn, p);
      }
    }
    for (int j = 0; j < d.J(); ++j) {
      if (++z[j] <= d.total(j)) break;
      z[j] = 0;
    }
  }
}

// p-value of the test for the stratum totals z and the statistic S
inline double pvalue(const Design& d, const Test& t, const int* z, int S,
                     const std::vector<double>& lf) {
  if (t.type == FISHER) {
    const Dist f = null_dist(d, z, lf);
    if (S < f.lo || S > f.hi()) throw std::runtime_error("S is outside its range");
    return fisher_pvalues(f, t.alternative, t.tsmethod)[S - f.lo];
  }
  const bool cochran = t.type == COCHRAN || t.type == COCHRAN_CC;
  const bool cc = t.type == CMH_CC || t.type == COCHRAN_CC;
  double E, V;
  null_moments(d, z, cochran, E, V);
  return asymptotic_pvalue(S, E, V, cc, t.alternative);
}

// Distribution of the number x of group 1 responders of a stratum jointly with its number
// z of responders, unnormalized in z: g(z, x) = b(x; m1, p1) b(z - x; m2, p2) for
// x = max(0, z - m2), ..., min(z, m1)
inline Dist joint_slice(const std::vector<double>& b1, const std::vector<double>& b2, int z) {
  const int m1 = static_cast<int>(b1.size()) - 1, m2 = static_cast<int>(b2.size()) - 1;
  Dist g;
  g.lo = std::max(0, z - m2);
  const int hi = std::min(z, m1);
  g.p.resize(hi - g.lo + 1);
  for (int x = g.lo; x <= hi; ++x) g.p[x - g.lo] = b1[x] * b2[z - x];
  return g;
}

// Probability of rejection over the strata j, ..., J - 1 given the totals of the earlier
// strata in z and the unnormalized distribution A of the partial sum of S over them
inline void power_rec(int j, const Dist& A, const Design& d, const Test& t,
                      const std::vector<std::vector<double>>& b1,
                      const std::vector<std::vector<double>>& b2, std::vector<int>& z,
                      const std::vector<double>& lf, double& total) {
  if (j < d.J() - 1) {
    for (int zj = 0; zj <= d.total(j); ++zj) {
      z[j] = zj;
      power_rec(j + 1, convolve(A, joint_slice(b1[j], b2[j], zj)), d, t, b1, b2, z, lf,
                total);
    }
    return;
  }
  // Last stratum: tail sums of the partial distribution, up[i] = A(S' >= A.lo + i) and
  // down[i] = A(S' <= A.lo + i)
  const int n = static_cast<int>(A.p.size());
  std::vector<double> up(n + 1, 0.0), down(n, 0.0);
  for (int i = n - 1; i >= 0; --i) up[i] = up[i + 1] + A.p[i];
  double acc = 0.0;
  for (int i = 0; i < n; ++i) {
    acc += A.p[i];
    down[i] = acc;
  }
  for (int zl = 0; zl <= d.total(j); ++zl) {
    z[j] = zl;
    const Dist g = joint_slice(b1[j], b2[j], zl);
    const Crit c = critical(d, t, z.data(), lf);
    double s = 0.0;
    for (int i = 0; i < static_cast<int>(g.p.size()); ++i) {
      const double gx = g.p[i];
      if (gx == 0.0) continue;
      const int x = g.lo + i;
      // A(S' >= c.hi - x) + A(S' <= c.lo - x)
      const int u = c.hi - x - A.lo;
      const int l = c.lo - x - A.lo;
      double r = 0.0;
      if (u <= 0) r += up[0]; else if (u < n) r += up[u];
      if (l >= n - 1) r += down[n - 1]; else if (l >= 0) r += down[l];
      s += gx * r;
    }
    total += s;
  }
}

// Rejection probability of a fixed design with probabilities p1[j] and p2[j], for any
// number of strata, by enumeration of the stratum totals. The distribution of the partial
// sum of S over the first strata is carried down the enumeration, and at the last stratum
// the probability of rejection is read from its tail sums
inline double power_fixed(const Design& d, const Test& t, const double* p1, const double* p2,
                          const std::vector<double>& lf) {
  const int J = d.J();
  std::vector<std::vector<double>> b1(J), b2(J);
  for (int j = 0; j < J; ++j) {
    binom_row(d.N1[j], p1[j], lf, b1[j]);
    binom_row(d.N2[j], p2[j], lf, b2[j]);
  }
  std::vector<int> z(J, 0);
  double total = 0.0;
  Dist A0;
  A0.lo = 0;
  A0.p.assign(1, 1.0);
  power_rec(0, A0, d, t, b1, b2, z, lf, total);
  return total;
}

// ---------------------------------------------------------------------------------------
// Two-stratum designs with blinded sample size re-estimation. The interim analysis has
// n1[j] and n2[j] patients in stratum j, and the final design used after the interim
// pooled numbers of responders (s_1, s_2) is designs[id[s_1 + (n1[0] + n2[0] + 1) s_2]].
// Each final design carries the critical values of every vector of its stratum totals.

struct FinalDesign {
  Design d;
  std::vector<int> lo;
  std::vector<int> hi;
};

// Joint distribution of the final number x of group 1 responders and the second-stage
// number w of responders of a stratum, given s interim responders, unnormalized by the
// probability of s: sum_a b(a; n1, p1) b(s - a; n2, p2) b(x - a; m1, p1) b(w - x + a; m2, p2)
// for x = 0, ..., n1 + m1 and w = 0, ..., m1 + m2, stored as J[x * (m1 + m2 + 1) + w]
inline std::vector<double> stage_joint(int n1, int n2, int m1, int m2, int s,
                                       const std::vector<double>& bi1,
                                       const std::vector<double>& bi2,
                                       const std::vector<double>& bs1,
                                       const std::vector<double>& bs2) {
  const int W = m1 + m2 + 1;
  std::vector<double> J((n1 + m1 + 1) * W, 0.0);
  for (int a = std::max(0, s - n2); a <= std::min(s, n1); ++a) {
    const double wa = bi1[a] * bi2[s - a];
    if (wa == 0.0) continue;
    for (int u = 0; u <= m1; ++u) {
      const double wu = wa * bs1[u];
      if (wu == 0.0) continue;
      double* row = &J[(a + u) * W + u];
      for (int v = 0; v <= m2; ++v) row[v] += wu * bs2[v];
    }
  }
  return J;
}

// Rejection probability of the re-estimation design for the probabilities p1[j], p2[j]
inline double bssr_power2(const std::vector<int>& n1, const std::vector<int>& n2,
                          const std::vector<FinalDesign>& designs, const std::vector<int>& id,
                          const double* p1, const double* p2, const std::vector<double>& lf) {
  const int S0 = n1[0] + n2[0], S1 = n1[1] + n2[1];
  std::vector<double> bi[2][2];
  for (int j = 0; j < 2; ++j) {
    binom_row(n1[j], p1[j], lf, bi[j][0]);
    binom_row(n2[j], p2[j], lf, bi[j][1]);
  }
  // Interim cells of each design, ordered by s_2
  std::vector<std::vector<std::pair<int, int>>> cells(designs.size());
  for (int s1 = 0; s1 <= S1; ++s1) {
    for (int s0 = 0; s0 <= S0; ++s0) cells[id[s0 + (S0 + 1) * s1]].push_back({s0, s1});
  }
  double total = 0.0;
  for (std::size_t k = 0; k < designs.size(); ++k) {
    if (cells[k].empty()) continue;
    const Design& d = designs[k].d;
    const int m1[2] = {d.N1[0] - n1[0], d.N1[1] - n1[1]};
    const int m2[2] = {d.N2[0] - n2[0], d.N2[1] - n2[1]};
    std::vector<double> bs[2][2];
    for (int j = 0; j < 2; ++j) {
      binom_row(m1[j], p1[j], lf, bs[j][0]);
      binom_row(m2[j], p2[j], lf, bs[j][1]);
    }
    const int W0 = m1[0] + m2[0] + 1, W1 = m1[1] + m2[1] + 1;
    const int X1 = d.N1[1] + 1;
    const int Z0 = d.total(0) + 1;
    std::map<int, std::vector<double>> J0;
    int cur = -1;
    std::vector<double> up, down;  // tails of stratum 2 in x, for each w
    for (const auto& cell : cells[k]) {
      const int s0 = cell.first, s1 = cell.second;
      if (s1 != cur) {
        cur = s1;
        const std::vector<double> J1 = stage_joint(n1[1], n2[1], m1[1], m2[1], s1, bi[1][0],
                                                   bi[1][1], bs[1][0], bs[1][1]);
        // up[y * W1 + w] = sum_{x >= y} J1(x, w), down[y * W1 + w] = sum_{x <= y} J1(x, w)
        up.assign((X1 + 1) * W1, 0.0);
        down.assign(X1 * W1, 0.0);
        for (int w = 0; w < W1; ++w) {
          double acc = 0.0;
          for (int x = X1 - 1; x >= 0; --x) {
            acc += J1[x * W1 + w];
            up[x * W1 + w] = acc;
          }
          acc = 0.0;
          for (int x = 0; x < X1; ++x) {
            acc += J1[x * W1 + w];
            down[x * W1 + w] = acc;
          }
        }
      }
      auto it = J0.find(s0);
      if (it == J0.end()) {
        it = J0.emplace(s0, stage_joint(n1[0], n2[0], m1[0], m2[0], s0, bi[0][0], bi[0][1],
                                        bs[0][0], bs[0][1])).first;
      }
      const std::vector<double>& J0s = it->second;
      const int a_lo = std::max(0, s0 - n2[0]), a_hi = std::min(s0, n1[0]);
      double crp = 0.0;
      for (int w0 = 0; w0 < W0; ++w0) {
        const int z0 = s0 + w0;
        const int x_lo = a_lo + std::max(0, w0 - m2[0]), x_hi = a_hi + std::min(w0, m1[0]);
        for (int x = x_lo; x <= x_hi; ++x) {
          const double h = J0s[x * W0 + w0];
          if (h == 0.0) continue;
          double acc = 0.0;
          for (int w1 = 0; w1 < W1; ++w1) {
            const std::size_t cell_z = z0 + static_cast<std::size_t>(Z0) * (s1 + w1);
            const int u = designs[k].hi[cell_z] - x;
            const int l = designs[k].lo[cell_z] - x;
            if (u <= 0) acc += up[w1]; else if (u < X1) acc += up[u * W1 + w1];
            if (l >= X1 - 1) acc += down[(X1 - 1) * W1 + w1];
            else if (l >= 0) acc += down[l * W1 + w1];
          }
          crp += h * acc;
        }
      }
      total += crp;
    }
  }
  return total;
}

// ---------------------------------------------------------------------------------------
// Type I error rate of a two-stratum re-estimation design as a polynomial in the common
// response probabilities theta_1 and theta_2 of the strata. Under the null hypothesis the
// numbers of group 1 responders given the interim and second-stage numbers of responders
// of a stratum are hypergeometric whatever theta, so the rate is
// sum_k sum_z Q_k(z) B_{z_1, T_1}(theta_1) B_{z_2, T_2}(theta_2), where T_j is the total of
// stratum j in design k, B the Bernstein basis polynomials, and Q_k(z) the sum over the
// interim cells s that lead to design k of the hypergeometric probability of s given z
// (the split of the stratum responders between the two stages) times the probability of
// rejection given s and z. The coefficients of all designs, raised to the largest degree
// in each variable, are returned as a (D_1 + 1) by (D_2 + 1) matrix in row major order.

// Coefficients c of degree n raised to degree D in one variable:
// c'_i = sum_k c_k choose(n, k) choose(D - n, i - k) / choose(D, i)
inline void elevate_rows(const std::vector<double>& c, int n, int D, int stride, int count,
                         const std::vector<double>& lf, std::vector<double>& out,
                         int out_stride) {
  // c[k * stride + r] for k = 0..n, r = 0..count-1; out[i * out_stride + r]
  for (int i = 0; i <= D; ++i) {
    const int k_lo = std::max(0, i - (D - n)), k_hi = std::min(i, n);
    for (int k = k_lo; k <= k_hi; ++k) {
      const double w = std::exp(lf[n] - lf[k] - lf[n - k] + lf[D - n] - lf[i - k] -
                                lf[D - n - i + k] - (lf[D] - lf[i] - lf[D - i]));
      const double* src = &c[k * stride];
      double* dst = &out[i * out_stride];
      for (int r = 0; r < count; ++r) dst[r] += w * src[r];
    }
  }
}

inline std::vector<double> tie_coef2(const std::vector<int>& n1, const std::vector<int>& n2,
                                     const std::vector<FinalDesign>& designs,
                                     const std::vector<int>& id, int& D0, int& D1,
                                     const std::vector<double>& lf) {
  const int S0 = n1[0] + n2[0], S1 = n1[1] + n2[1];
  std::vector<std::vector<std::pair<int, int>>> cells(designs.size());
  for (int s1 = 0; s1 <= S1; ++s1) {
    for (int s0 = 0; s0 <= S0; ++s0) cells[id[s0 + (S0 + 1) * s1]].push_back({s0, s1});
  }
  D0 = 0;
  D1 = 0;
  for (std::size_t k = 0; k < designs.size(); ++k) {
    if (cells[k].empty()) continue;
    D0 = std::max(D0, designs[k].d.total(0));
    D1 = std::max(D1, designs[k].d.total(1));
  }
  std::vector<double> coef((D0 + 1) * (D1 + 1), 0.0);
  for (std::size_t k = 0; k < designs.size(); ++k) {
    if (cells[k].empty()) continue;
    const Design& d = designs[k].d;
    const int T0 = d.total(0), T1 = d.total(1);
    const int m1[2] = {d.N1[0] - n1[0], d.N1[1] - n1[1]};
    const int m2[2] = {d.N2[0] - n2[0], d.N2[1] - n2[1]};
    const int W0 = m1[0] + m2[0] + 1, W1 = m1[1] + m2[1] + 1;
    std::vector<double> Q((T0 + 1) * (T1 + 1), 0.0);  // Q[z0 * (T1 + 1) + z1]
    // Null distributions of x for stratum 1 given (s0, w0), cached by s0
    std::map<int, std::vector<Dist>> phi0;
    int cur = -1;
    std::vector<std::vector<double>> up1, down1;  // tails for stratum 2 by w1
    std::vector<int> lo1(W1);
    for (const auto& cell : cells[k]) {
      const int s0 = cell.first, s1 = cell.second;
      if (s1 != cur) {
        cur = s1;
        const Dist a = hyper(n1[1], n2[1], s1, lf);
        up1.assign(W1, std::vector<double>());
        down1.assign(W1, std::vector<double>());
        for (int w = 0; w < W1; ++w) {
          const Dist f = convolve(a, hyper(m1[1], m2[1], w, lf));
          const int n = static_cast<int>(f.p.size());
          lo1[w] = f.lo;
          up1[w].assign(n + 1, 0.0);
          down1[w].assign(n, 0.0);
          for (int i = n - 1; i >= 0; --i) up1[w][i] = up1[w][i + 1] + f.p[i];
          double acc = 0.0;
          for (int i = 0; i < n; ++i) {
            acc += f.p[i];
            down1[w][i] = acc;
          }
        }
      }
      auto it = phi0.find(s0);
      if (it == phi0.end()) {
        std::vector<Dist> v(W0);
        const Dist a = hyper(n1[0], n2[0], s0, lf);
        for (int w = 0; w < W0; ++w) v[w] = convolve(a, hyper(m1[0], m2[0], w, lf));
        it = phi0.emplace(s0, v).first;
      }
      const std::vector<Dist>& ph = it->second;
      // Stage split probabilities choose(n, s) choose(m, w) / choose(n + m, s + w)
      const int nn0 = S0, mm0 = W0 - 1, nn1 = S1, mm1 = W1 - 1;
      for (int w0 = 0; w0 < W0; ++w0) {
        const int z0 = s0 + w0;
        const double h0 = std::exp(lf[nn0] - lf[s0] - lf[nn0 - s0] + lf[mm0] - lf[w0] -
                                   lf[mm0 - w0] - (lf[T0] - lf[z0] - lf[T0 - z0]));
        const Dist& f0 = ph[w0];
        for (int w1 = 0; w1 < W1; ++w1) {
          const int z1 = s1 + w1;
          const double h1 = std::exp(lf[nn1] - lf[s1] - lf[nn1 - s1] + lf[mm1] - lf[w1] -
                                     lf[mm1 - w1] - (lf[T1] - lf[z1] - lf[T1 - z1]));
          const std::size_t cz = z0 + static_cast<std::size_t>(T0 + 1) * z1;
          const int chi = designs[k].hi[cz], clo = designs[k].lo[cz];
          const std::vector<double>& up = up1[w1];
          const std::vector<double>& down = down1[w1];
          const int n = static_cast<int>(down.size());
          double T = 0.0;
          for (int i = 0; i < static_cast<int>(f0.p.size()); ++i) {
            const double fx = f0.p[i];
            if (fx == 0.0) continue;
            const int x = f0.lo + i;
            const int u = chi - x - lo1[w1];
            const int l = clo - x - lo1[w1];
            double r = 0.0;
            if (u <= 0) r += up[0]; else if (u < n) r += up[u];
            if (l >= n - 1) r += down[n - 1]; else if (l >= 0) r += down[l];
            T += fx * r;
          }
          Q[z0 * (T1 + 1) + z1] += h0 * h1 * T;
        }
      }
    }
    // Raise the degrees to (D0, D1) and add
    std::vector<double> tmp((D0 + 1) * (T1 + 1), 0.0);
    elevate_rows(Q, T0, D0, T1 + 1, T1 + 1, lf, tmp, T1 + 1);
    // Along the second variable: each row of tmp is a polynomial of degree T1
    for (int i = 0; i <= D0; ++i) {
      std::vector<double> row(&tmp[i * (T1 + 1)], &tmp[i * (T1 + 1)] + (T1 + 1));
      std::vector<double> out(D1 + 1, 0.0);
      elevate_rows(row, T1, D1, 1, 1, lf, out, 1);
      for (int r = 0; r <= D1; ++r) coef[i * (D1 + 1) + r] += out[r];
    }
  }
  return coef;
}

// ---------------------------------------------------------------------------------------
// Bivariate Bernstein polynomials on a rectangle: coefficients c[i * (d2 + 1) + k] of
// B_{i,d1}(x) B_{k,d2}(y)

struct Poly2 {
  int d1;
  int d2;
  std::vector<double> c;
};

// de Casteljau at t along one axis of a vector of coefficients of degree n: left holds the
// coefficients on [0, t] and right those on [t, 1] of the rescaled variable
inline void casteljau(const std::vector<double>& c, double t, std::vector<double>& left,
                      std::vector<double>& right) {
  const std::size_t n = c.size();
  std::vector<double> w(c);
  left.resize(n);
  right.resize(n);
  for (std::size_t j = 0; j < n; ++j) {
    left[j] = w[0];
    right[n - 1 - j] = w[n - 1 - j];
    for (std::size_t i = 0; i + 1 < n - j; ++i) w[i] = (1.0 - t) * w[i] + t * w[i + 1];
  }
}

// Split along axis 1 (x) or 2 (y) at the parameter t of the current rectangle
inline void split2(const Poly2& P, int axis, double t, Poly2& L, Poly2& R) {
  L = P;
  R = P;
  const int m = P.d2 + 1;
  std::vector<double> col, l, r;
  if (axis == 1) {
    col.resize(P.d1 + 1);
    for (int k = 0; k < m; ++k) {
      for (int i = 0; i <= P.d1; ++i) col[i] = P.c[i * m + k];
      casteljau(col, t, l, r);
      for (int i = 0; i <= P.d1; ++i) {
        L.c[i * m + k] = l[i];
        R.c[i * m + k] = r[i];
      }
    }
  } else {
    for (int i = 0; i <= P.d1; ++i) {
      col.assign(P.c.begin() + i * m, P.c.begin() + (i + 1) * m);
      casteljau(col, t, l, r);
      std::copy(l.begin(), l.end(), L.c.begin() + i * m);
      std::copy(r.begin(), r.end(), R.c.begin() + i * m);
    }
  }
}

// Coefficients of the polynomial restricted to [a1, b1] x [a2, b2] within [0, 1]^2
inline Poly2 restrict2(const Poly2& P, double a1, double b1, double a2, double b2) {
  Poly2 Q = P, L, R;
  if (b1 < 1.0) {
    split2(Q, 1, b1, L, R);
    Q = L;
  }
  if (a1 > 0.0) {
    split2(Q, 1, a1 / b1, L, R);
    Q = R;
  }
  if (b2 < 1.0) {
    split2(Q, 2, b2, L, R);
    Q = L;
  }
  if (a2 > 0.0) {
    split2(Q, 2, a2 / b2, L, R);
    Q = R;
  }
  return Q;
}

// Value at (x, y) in [0, 1]^2
inline double eval2(const Poly2& P, double x, double y, const std::vector<double>& lf) {
  std::vector<double> bx, by;
  binom_row(P.d1, x, lf, bx);
  binom_row(P.d2, y, lf, by);
  double s = 0.0;
  for (int i = 0; i <= P.d1; ++i) {
    if (bx[i] == 0.0) continue;
    double r = 0.0;
    const double* row = &P.c[i * (P.d2 + 1)];
    for (int k = 0; k <= P.d2; ++k) r += row[k] * by[k];
    s += bx[i] * r;
  }
  return s;
}

struct Max2 {
  double x;
  double y;
  double value;
  double bound;
  bool complete;
};

// Maximum over [0, 1]^2. On any rectangle the polynomial lies between the smallest and the
// largest of its coefficients there, and the corner coefficients are its values at the
// corners. Rectangles are halved along their longer side until the largest coefficient
// does not exceed the best value found by more than tol; the largest coefficient of the
// rectangles set aside bounds the maximum. At most max_pieces rectangles are split
inline Max2 bernstein_max2(const Poly2& P, double tol, int max_pieces) {
  struct Piece {
    Poly2 P;
    double a1, b1, a2, b2;
  };
  Max2 out;
  out.value = -std::numeric_limits<double>::infinity();
  out.x = 0.0;
  out.y = 0.0;
  out.complete = true;
  auto corners = [&out](const Piece& p) {
    const int m = p.P.d2 + 1, n = p.P.d1;
    const double v[4] = {p.P.c[0], p.P.c[p.P.d2], p.P.c[n * m], p.P.c[n * m + p.P.d2]};
    const double xs[4] = {p.a1, p.a1, p.b1, p.b1}, ys[4] = {p.a2, p.b2, p.a2, p.b2};
    for (int i = 0; i < 4; ++i) {
      if (v[i] > out.value) {
        out.value = v[i];
        out.x = xs[i];
        out.y = ys[i];
      }
    }
  };
  std::vector<Piece> stack;
  Piece p0{P, 0.0, 1.0, 0.0, 1.0};
  corners(p0);
  double bound = out.value;
  stack.push_back(p0);
  int split = 0;
  while (!stack.empty()) {
    Piece p = std::move(stack.back());
    stack.pop_back();
    const double cmax = *std::max_element(p.P.c.begin(), p.P.c.end());
    if (cmax <= out.value + tol) {
      bound = std::max(bound, cmax);
      continue;
    }
    if (split >= max_pieces) {
      out.complete = false;
      bound = std::max(bound, cmax);
      continue;
    }
    ++split;
    const int axis = (p.b1 - p.a1 >= p.b2 - p.a2) ? 1 : 2;
    Piece L = p, R = p;
    split2(p.P, axis, 0.5, L.P, R.P);
    if (axis == 1) {
      const double mid = 0.5 * (p.a1 + p.b1);
      L.b1 = mid;
      R.a1 = mid;
    } else {
      const double mid = 0.5 * (p.a2 + p.b2);
      L.b2 = mid;
      R.a2 = mid;
    }
    corners(L);
    corners(R);
    stack.push_back(std::move(R));
    stack.push_back(std::move(L));
  }
  out.bound = std::max(bound, out.value);
  return out;
}

}  // namespace strat

#endif
