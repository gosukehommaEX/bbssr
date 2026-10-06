#include <Rcpp.h>
#include <algorithm>
#include <utility>
#include <vector>
using namespace Rcpp;

// Bernstein coefficients of the same polynomial on the two halves of its interval, by the
// de Casteljau algorithm at the midpoint. The last coefficient of the left half and the
// first of the right half are the value of the polynomial at the midpoint.
static void split_half(const std::vector<double>& c, std::vector<double>& left,
                       std::vector<double>& right) {
  const std::size_t n = c.size();
  std::vector<double> w(c);
  left.resize(n);
  right.resize(n);
  for (std::size_t j = 0; j < n; ++j) {
    left[j] = w[0];
    right[n - 1 - j] = w[n - 1 - j];
    for (std::size_t i = 0; i + 1 < n - j; ++i) w[i] = 0.5 * (w[i] + w[i + 1]);
  }
}

// Maximum over [0, 1] of the polynomial sum_k coef[k] B_{k,n}(t), where B_{k,n} are the
// Bernstein basis polynomials of degree n = length(coef) - 1.
//
// On any interval the polynomial lies between the smallest and the largest of its
// Bernstein coefficients on that interval, and the first and last coefficients are its
// values at the two ends. The interval is halved repeatedly; a piece is set aside once
// its largest coefficient does not exceed the best value found by more than tol, and the
// largest coefficient of the pieces set aside is an upper bound of the maximum. At most
// max_pieces pieces are split; if this limit is reached the remaining pieces are set aside
// as they are, the bound remains valid but may exceed the best value by more than tol,
// and complete is FALSE.
//
// Returns the location x and the value of the best point found, the upper bound, and
// whether the search ended with the bound within tol of the value.
// [[Rcpp::export]]
List bernstein_max(NumericVector coef, double tol, int max_pieces) {
  struct Piece {
    std::vector<double> c;
    double a;
    double b;
  };
  const std::size_t n = coef.size();
  if (n == 0) stop("coef must not be empty");
  Piece p0;
  p0.c.assign(coef.begin(), coef.end());
  p0.a = 0.0;
  p0.b = 1.0;
  double best = p0.c[0];
  double x = 0.0;
  if (p0.c[n - 1] > best) {
    best = p0.c[n - 1];
    x = 1.0;
  }
  double bound = best;
  bool complete = true;
  int split = 0;
  std::vector<Piece> stack;
  stack.push_back(p0);
  while (!stack.empty()) {
    Piece p = std::move(stack.back());
    stack.pop_back();
    const double cmax = *std::max_element(p.c.begin(), p.c.end());
    if (cmax <= best + tol) {
      bound = std::max(bound, cmax);
      continue;
    }
    if (split >= max_pieces) {
      complete = false;
      bound = std::max(bound, cmax);
      continue;
    }
    ++split;
    Piece l;
    Piece r;
    const double mid = 0.5 * (p.a + p.b);
    split_half(p.c, l.c, r.c);
    l.a = p.a;
    l.b = mid;
    r.a = mid;
    r.b = p.b;
    if (l.c[n - 1] > best) {
      best = l.c[n - 1];
      x = mid;
    }
    stack.push_back(std::move(r));
    stack.push_back(std::move(l));
  }
  bound = std::max(bound, best);
  return List::create(_["x"] = x, _["value"] = best, _["bound"] = bound,
                      _["complete"] = complete);
}

// Matrix M of dimension (N + 1) by (N + 1) with M[X, j] = P(Bin(N - j, a) + Bin(j, c) = X),
// the two binomial variables being independent. If p(t) = a (1 - t) + c t, then
// dbinom(X, N, p(t)) = sum_j M[X, j] B_{j,N}(t), so M converts a binomial probability in
// a response probability that is linear in t into Bernstein coefficients in t.
// [[Rcpp::export]]
NumericMatrix binom_conv_matrix(int N, double a, double c) {
  if (N < 0) stop("N must be non-negative");
  NumericMatrix M(N + 1, N + 1);
  std::vector<double> u;
  std::vector<double> v;
  for (int j = 0; j <= N; ++j) {
    u.resize(N - j + 1);
    v.resize(j + 1);
    for (int k = 0; k <= N - j; ++k) u[k] = R::dbinom(k, N - j, a, 0);
    for (int k = 0; k <= j; ++k) v[k] = R::dbinom(k, j, c, 0);
    for (int k = 0; k <= N - j; ++k) {
      if (u[k] == 0.0) continue;
      for (int l = 0; l <= j; ++l) M(k + l, j) += u[k] * v[l];
    }
  }
  return M;
}
