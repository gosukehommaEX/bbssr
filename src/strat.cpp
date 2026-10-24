#include <Rcpp.h>
#include <stdexcept>
#include <vector>
#include "strat_core.h"
using namespace Rcpp;

// Rcpp interfaces of the computations for stratified designs, see strat_core.h.

// A design from the per-stratum group sizes, which must be positive
static strat::Design make_design(const IntegerVector& N1, const IntegerVector& N2) {
  if (N1.size() != N2.size() || N1.size() < 1) {
    stop("N1 and N2 must have the same positive length");
  }
  strat::Design d;
  for (R_xlen_t j = 0; j < N1.size(); ++j) {
    if (N1[j] == NA_INTEGER || N2[j] == NA_INTEGER || N1[j] < 1 || N2[j] < 1) {
      stop("the group sizes of every stratum must be positive");
    }
    d.N1.push_back(N1[j]);
    d.N2.push_back(N2[j]);
  }
  return d;
}

static strat::Test make_test(int type, int alternative, int tsmethod, double alpha) {
  if (type < 0 || type > 7) stop("type must be 0, ..., 7");
  if (alternative < 0 || alternative > 2) stop("alternative must be 0, 1 or 2");
  if (tsmethod < 0 || tsmethod > 2) stop("tsmethod must be 0, 1 or 2");
  if (!(alpha > 0.0 && alpha < 1.0)) stop("alpha must lie in (0, 1)");
  strat::Test t{type, alternative, tsmethod, alpha};
  return t;
}

static int largest_total(const strat::Design& d) {
  int m = 0;
  for (int j = 0; j < d.J(); ++j) m = std::max(m, d.total(j));
  return m;
}

// Total size of the design, the largest argument of the log factorials of the
// unstratified tests, which pool the strata
static int total_size(const strat::Design& d) {
  int m = 0;
  for (int j = 0; j < d.J(); ++j) m += d.total(j);
  return m;
}

// Critical values for every vector of stratum totals, in the storage order of
// strat_core.h: reject when S <= lo or S >= hi
// [[Rcpp::export]]
List strat_crit(IntegerVector N1, IntegerVector N2, int type, int alternative, int tsmethod,
                double alpha) {
  const strat::Design d = make_design(N1, N2);
  const strat::Test t = make_test(type, alternative, tsmethod, alpha);
  const std::vector<double> lf = strat::log_factorials(total_size(d));
  std::vector<int> lo, hi;
  try {
    strat::critical_grid(d, t, lf, lo, hi);
  } catch (const std::exception& e) {
    stop(e.what());
  }
  return List::create(_["lo"] = wrap(lo), _["hi"] = wrap(hi));
}

// Exact power of a fixed design for each row of p1 and p2 (one column per stratum)
// [[Rcpp::export]]
NumericVector strat_power(IntegerVector N1, IntegerVector N2, NumericMatrix p1,
                          NumericMatrix p2, int type, int alternative, int tsmethod,
                          double alpha) {
  const strat::Design d = make_design(N1, N2);
  const strat::Test t = make_test(type, alternative, tsmethod, alpha);
  const int J = d.J();
  if (p1.ncol() != J || p2.ncol() != J || p1.nrow() != p2.nrow()) {
    stop("p1 and p2 must have one row per scenario and one column per stratum");
  }
  const std::vector<double> lf = strat::log_factorials(total_size(d));
  NumericVector out(p1.nrow());
  std::vector<double> q1(J), q2(J);
  for (int i = 0; i < p1.nrow(); ++i) {
    for (int j = 0; j < J; ++j) {
      q1[j] = p1(i, j);
      q2[j] = p2(i, j);
      if (!(q1[j] >= 0.0 && q1[j] <= 1.0 && q2[j] >= 0.0 && q2[j] <= 1.0)) {
        stop("the response probabilities must lie in [0, 1]");
      }
    }
    try {
      out[i] = strat::power_fixed(d, t, q1.data(), q2.data(), lf);
    } catch (const std::exception& e) {
      stop(e.what());
    }
  }
  return out;
}

// Final designs of a two-stratum re-estimation design from a list whose elements hold
// N1, N2 (length 2) and the critical values lo and hi of strat_crit()
static std::vector<strat::FinalDesign> final_designs(const List& designs,
                                                     const IntegerVector& n1,
                                                     const IntegerVector& n2, int& Tmax) {
  if (n1.size() != 2 || n2.size() != 2) stop("n1 and n2 must have length 2");
  for (int j = 0; j < 2; ++j) {
    if (n1[j] == NA_INTEGER || n2[j] == NA_INTEGER || n1[j] < 0 || n2[j] < 0) {
      stop("the interim sizes must be non-negative");
    }
  }
  std::vector<strat::FinalDesign> out;
  Tmax = n1[0] + n2[0] + n1[1] + n2[1];
  for (R_xlen_t k = 0; k < designs.size(); ++k) {
    List e = designs[k];
    strat::FinalDesign f;
    f.d = make_design(as<IntegerVector>(e["N1"]), as<IntegerVector>(e["N2"]));
    if (f.d.J() != 2) stop("the final designs must have two strata");
    for (int j = 0; j < 2; ++j) {
      if (f.d.N1[j] < n1[j] || f.d.N2[j] < n2[j]) {
        stop("a final group size is below the interim size");
      }
    }
    IntegerVector lo = e["lo"], hi = e["hi"];
    if (static_cast<std::size_t>(lo.size()) != f.d.cells() ||
        static_cast<std::size_t>(hi.size()) != f.d.cells()) {
      stop("the critical values do not match the design");
    }
    f.lo.assign(lo.begin(), lo.end());
    f.hi.assign(hi.begin(), hi.end());
    Tmax = std::max(Tmax, largest_total(f.d));
    out.push_back(f);
  }
  return out;
}

static void check_id(const IntegerVector& id, const IntegerVector& n1, const IntegerVector& n2,
                     std::size_t K) {
  const R_xlen_t cells = static_cast<R_xlen_t>(n1[0] + n2[0] + 1) * (n1[1] + n2[1] + 1);
  if (id.size() != cells) stop("id must have one element per interim cell");
  for (R_xlen_t i = 0; i < cells; ++i) {
    if (id[i] == NA_INTEGER || id[i] < 0 || id[i] >= static_cast<int>(K)) {
      stop("id must give zero based positions in designs");
    }
  }
}

// Exact power of a two-stratum re-estimation design for each row of p1 and p2. id[s1 +
// (n1[0] + n2[0] + 1) * s2] is the zero based position in designs of the final design used
// after s1 and s2 interim responders in the two strata
// [[Rcpp::export]]
NumericVector strat_bssr_power(IntegerVector n1, IntegerVector n2, List designs,
                               IntegerVector id, NumericMatrix p1, NumericMatrix p2) {
  int Tmax = 0;
  const std::vector<strat::FinalDesign> ds = final_designs(designs, n1, n2, Tmax);
  check_id(id, n1, n2, ds.size());
  if (p1.ncol() != 2 || p2.ncol() != 2 || p1.nrow() != p2.nrow()) {
    stop("p1 and p2 must have one row per scenario and two columns");
  }
  const std::vector<int> a1(n1.begin(), n1.end()), a2(n2.begin(), n2.end()),
    ids(id.begin(), id.end());
  const std::vector<double> lf = strat::log_factorials(Tmax);
  NumericVector out(p1.nrow());
  for (int i = 0; i < p1.nrow(); ++i) {
    const double q1[2] = {p1(i, 0), p1(i, 1)}, q2[2] = {p2(i, 0), p2(i, 1)};
    for (int j = 0; j < 2; ++j) {
      if (!(q1[j] >= 0.0 && q1[j] <= 1.0 && q2[j] >= 0.0 && q2[j] <= 1.0)) {
        stop("the response probabilities must lie in [0, 1]");
      }
    }
    out[i] = strat::bssr_power2(a1, a2, ds, ids, q1, q2, lf);
  }
  return out;
}

// Coefficients of the type I error rate of a two-stratum re-estimation design in the
// Bernstein basis B_{i,D1}(theta_1) B_{k,D2}(theta_2), as a (D1 + 1) by (D2 + 1) matrix
// [[Rcpp::export]]
NumericMatrix strat_tie_coef(IntegerVector n1, IntegerVector n2, List designs,
                             IntegerVector id) {
  int Tmax = 0;
  const std::vector<strat::FinalDesign> ds = final_designs(designs, n1, n2, Tmax);
  check_id(id, n1, n2, ds.size());
  const std::vector<int> a1(n1.begin(), n1.end()), a2(n2.begin(), n2.end()),
    ids(id.begin(), id.end());
  const std::vector<double> lf = strat::log_factorials(Tmax);
  int D0 = 0, D1 = 0;
  const std::vector<double> c = strat::tie_coef2(a1, a2, ds, ids, D0, D1, lf);
  NumericMatrix out(D0 + 1, D1 + 1);
  for (int i = 0; i <= D0; ++i) {
    for (int k = 0; k <= D1; ++k) out(i, k) = c[i * (D1 + 1) + k];
  }
  return out;
}

// Conditional rejection probabilities of a two-stratum re-estimation design under the null
// hypothesis, see strat::cond_reject2(). The points (theta1[p], theta2[p]), possibly none,
// are those of the decomposition of the type I error rate. The vectors per interim cell
// are in the storage order of id; w1 and w2 are the second-stage numbers of responders of
// the two strata
// [[Rcpp::export]]
List strat_cond_reject(IntegerVector n1, IntegerVector n2, List designs, IntegerVector id,
                       double alpha, bool full, NumericVector theta1, NumericVector theta2) {
  int Tmax = 0;
  const std::vector<strat::FinalDesign> ds = final_designs(designs, n1, n2, Tmax);
  check_id(id, n1, n2, ds.size());
  if (!(alpha > 0.0 && alpha < 1.0)) stop("alpha must lie in (0, 1)");
  if (theta1.size() != theta2.size()) stop("theta1 and theta2 must have the same length");
  for (R_xlen_t p = 0; p < theta1.size(); ++p) {
    if (!(theta1[p] >= 0.0 && theta1[p] <= 1.0 && theta2[p] >= 0.0 && theta2[p] <= 1.0)) {
      stop("the points must lie in the unit square");
    }
  }
  const std::vector<int> a1(n1.begin(), n1.end()), a2(n2.begin(), n2.end()),
    ids(id.begin(), id.end());
  const std::vector<double> t1(theta1.begin(), theta1.end()), t2(theta2.begin(), theta2.end());
  const std::vector<double> lf = strat::log_factorials(Tmax);
  const strat::CondReject cr = strat::cond_reject2(a1, a2, ds, ids, alpha, full, t1, t2, lf);
  const int ncell = static_cast<int>(cr.crp_max.size()), P = static_cast<int>(t1.size());
  NumericMatrix tie(ncell, P);
  for (int i = 0; i < ncell; ++i) {
    for (int p = 0; p < P; ++p) tie(i, p) = cr.tie_s[static_cast<std::size_t>(i) * P + p];
  }
  return List::create(_["crp.max"] = wrap(cr.crp_max), _["crp.min"] = wrap(cr.crp_min),
                      _["total.max"] = wrap(cr.total_max), _["exceed"] = wrap(cr.exceed),
                      _["exceed.total"] = wrap(cr.exceed_total),
                      _["w1.max"] = wrap(cr.w0_max), _["w2.max"] = wrap(cr.w1_max),
                      _["w1.total"] = wrap(cr.w0_total), _["w2.total"] = wrap(cr.w1_total),
                      _["tie.s"] = tie, _["s1"] = wrap(cr.s0), _["s2"] = wrap(cr.s1),
                      _["w1"] = wrap(cr.w0), _["w2"] = wrap(cr.w1), _["crp"] = wrap(cr.crp),
                      _["total"] = wrap(cr.total));
}

static strat::Poly2 make_poly(const NumericMatrix& coef) {
  if (coef.nrow() < 1 || coef.ncol() < 1) stop("coef must not be empty");
  strat::Poly2 P;
  P.d1 = coef.nrow() - 1;
  P.d2 = coef.ncol() - 1;
  P.c.resize(coef.nrow() * coef.ncol());
  for (int i = 0; i <= P.d1; ++i) {
    for (int k = 0; k <= P.d2; ++k) P.c[i * (P.d2 + 1) + k] = coef(i, k);
  }
  return P;
}

// Maximum over [a1, b1] x [a2, b2] of the polynomial with the coefficients coef over
// [0, 1]^2, see strat::bernstein_max2(). x and y are the location in the original variables
// [[Rcpp::export]]
List strat_bern_max(NumericMatrix coef, double a1, double b1, double a2, double b2,
                    double tol, int max_pieces) {
  if (!(0.0 <= a1 && a1 < b1 && b1 <= 1.0 && 0.0 <= a2 && a2 < b2 && b2 <= 1.0)) {
    stop("the rectangle must lie in the unit square and have a positive width");
  }
  const strat::Poly2 P = strat::restrict2(make_poly(coef), a1, b1, a2, b2);
  const strat::Max2 m = strat::bernstein_max2(P, tol, max_pieces);
  return List::create(_["x"] = a1 + (b1 - a1) * m.x, _["y"] = a2 + (b2 - a2) * m.y,
                      _["value"] = m.value, _["bound"] = m.bound,
                      _["complete"] = m.complete);
}

// Values of the polynomial with the coefficients coef at the points (x[i], y[i])
// [[Rcpp::export]]
NumericVector strat_bern_eval(NumericMatrix coef, NumericVector x, NumericVector y) {
  if (x.size() != y.size()) stop("x and y must have the same length");
  const strat::Poly2 P = make_poly(coef);
  const std::vector<double> lf = strat::log_factorials(std::max(P.d1, P.d2));
  NumericVector out(x.size());
  for (R_xlen_t i = 0; i < x.size(); ++i) {
    if (!(x[i] >= 0.0 && x[i] <= 1.0 && y[i] >= 0.0 && y[i] <= 1.0)) {
      stop("the points must lie in the unit square");
    }
    out[i] = strat::eval2(P, x[i], y[i], lf);
  }
  return out;
}

// p-values of the test for observed data: row i gives, for each stratum, the group sizes
// N1 and N2 and the numbers of responders x1 and x2 of the two groups
// [[Rcpp::export]]
NumericVector strat_pvalue_cases(IntegerMatrix N1, IntegerMatrix N2, IntegerMatrix x1,
                                 IntegerMatrix x2, int type, int alternative, int tsmethod) {
  const int n = N1.nrow(), J = N1.ncol();
  if (N2.nrow() != n || x1.nrow() != n || x2.nrow() != n || N2.ncol() != J ||
      x1.ncol() != J || x2.ncol() != J) {
    stop("N1, N2, x1 and x2 must have the same dimensions");
  }
  const strat::Test t = make_test(type, alternative, tsmethod, 0.5);
  int Tmax = 1;
  for (int i = 0; i < n; ++i) {
    int T = 0;
    for (int j = 0; j < J; ++j) T += N1(i, j) + N2(i, j);
    Tmax = std::max(Tmax, T);
  }
  const std::vector<double> lf = strat::log_factorials(Tmax);
  NumericVector out(n);
  strat::Design d;
  d.N1.resize(J);
  d.N2.resize(J);
  std::vector<int> z(J);
  for (int i = 0; i < n; ++i) {
    int S = 0;
    for (int j = 0; j < J; ++j) {
      d.N1[j] = N1(i, j);
      d.N2[j] = N2(i, j);
      if (d.N1[j] < 1 || d.N2[j] < 1) stop("the group sizes must be positive");
      if (x1(i, j) < 0 || x1(i, j) > d.N1[j] || x2(i, j) < 0 || x2(i, j) > d.N2[j]) {
        stop("the numbers of responders must lie between 0 and the group sizes");
      }
      z[j] = x1(i, j) + x2(i, j);
      S += x1(i, j);
    }
    out[i] = strat::pvalue(d, t, z.data(), S, lf);
  }
  return out;
}

// Largest p-value rejected at the level alpha and smallest p-value not rejected over all
// outcomes of the design, -1 and 2 when there is none
// [[Rcpp::export]]
NumericVector strat_level_bracket(IntegerVector N1, IntegerVector N2, int type,
                                  int alternative, int tsmethod, double alpha) {
  const strat::Design d = make_design(N1, N2);
  const strat::Test t = make_test(type, alternative, tsmethod, alpha);
  const std::vector<double> lf = strat::log_factorials(total_size(d));
  double pr = -1.0, pn = 2.0;
  strat::level_bracket(d, t, lf, pr, pn);
  return NumericVector::create(pr, pn);
}
