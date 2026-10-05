"""Independent reference values for the unit tests of bbssr.

Run by the developer in an environment with numpy and scipy:
    python3 tools/reference/reference_values.py > tools/reference/reference-values.txt
then check that the test files carry these values:
    python3 tools/reference/check_test_constants.py

The computations follow the definitions of the designs directly (binomial sums over the
interim and second-stage outcomes) and share no code with the package. Each block names
the test file whose expected values it produces.
"""
import math
import sys
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import beta as beta_dist
from scipy.stats import binom, hypergeom, norm

TOL = 1.4901161193847656e-08  # tolerance of fpCompare
OUT = []


def emit(test_file, key, values, rtol=1e-10):
    """Print one key. rtol is the relative tolerance used by check_test_constants.py."""
    values = [float(v) for v in np.atleast_1d(values)]
    OUT.append((test_file, key, values))
    print(f"{test_file}\t{key}\t{rtol:g}\t" + " ".join("%.15g" % v for v in values))
    sys.stdout.flush()


def ceil_tol(x):
    return math.ceil(x - 1e-9)


# ---------------------------------------------------------------------------------------
# Test statistics and rejection regions
def zstat(N1, N2):
    x1 = np.arange(N1 + 1)[:, None]
    x2 = np.arange(N2 + 1)[None, :]
    d = x1 / N1 - x2 / N2
    p = (x1 + x2) / (N1 + N2)
    with np.errstate(all="ignore"):
        z = d / np.sqrt(p * (1 - p) * (1 / N1 + 1 / N2))
    z[~np.isfinite(z)] = 0
    return z


def fisher_upper(N1, N2):
    P = np.full((N1 + 1, N2 + 1), np.nan)
    for s in range(N1 + N2 + 1):
        k = np.arange(max(0, s - N2), min(N1, s) + 1)
        d = hypergeom.pmf(k, N1 + N2, N1, s)
        P[k, s - k] = np.minimum(1, np.cumsum(d[::-1])[::-1])
    return P


def unconditional(stat, N1, N2, decreasing, G=100):
    s = stat.ravel()
    x1 = np.repeat(np.arange(N1 + 1), N2 + 1)
    x2 = np.tile(np.arange(N2 + 1), N1 + 1)
    order = np.argsort(-s if decreasing else s, kind="stable")
    so = s[order]
    ref = np.maximum(np.abs(so[:-1]), np.abs(so[1:]))
    new = np.abs(np.diff(so)) > 1e-10 * ref
    ends = np.r_[np.where(new)[0], len(so) - 1]
    last = np.repeat(ends, np.diff(np.r_[-1, ends]))
    th = np.linspace(0, 1, G)
    b1 = binom.pmf(np.arange(N1 + 1)[:, None], N1, th[None, :])
    b2 = binom.pmf(np.arange(N2 + 1)[:, None], N2, th[None, :])
    cum = np.cumsum(b1[x1[order]] * b2[x2[order]], axis=0)
    out = np.empty_like(s)
    out[order] = np.minimum(1, cum[last].max(axis=1))
    return out.reshape(stat.shape)


@lru_cache(maxsize=None)
def pvalues(N1, N2, test, alt):
    if alt == "two.sided":
        assert test == "Chisq"
        return np.minimum(2 * norm.sf(np.abs(zstat(N1, N2))), 1)
    if test == "Chisq":
        return norm.sf(zstat(N1, N2))
    if test == "Fisher":
        return fisher_upper(N1, N2)
    if test == "Z-pool":
        return unconditional(zstat(N1, N2), N1, N2, True)
    if test == "Boschloo":
        return unconditional(fisher_upper(N1, N2), N1, N2, False)
    raise ValueError(test)


def reject(N1, N2, test, alt, alpha):
    return pvalues(N1, N2, test, alt) < alpha - TOL


def power(p1, p2, N1, N2, test, alt, alpha):
    R = reject(N1, N2, test, alt, alpha).astype(float)
    return binom.pmf(np.arange(N1 + 1), N1, p1) @ R @ binom.pmf(np.arange(N2 + 1), N2, p2)


# ---------------------------------------------------------------------------------------
# Sample size rules
def exact_n2(p1, p2, r, alpha, tp, test, alt):
    pa = lambda n2: power(p1, p2, math.ceil(r * n2), n2, test, alt, alpha)
    ae = alpha / 2 if alt == "two.sided" else alpha
    p = (r * p1 + p2) / (1 + r)
    init = (1 + 1 / r) / (p1 - p2) ** 2 * (
        norm.ppf(ae) * math.sqrt(p * (1 - p))
        + norm.ppf(1 - tp) * math.sqrt((p1 * (1 - p1) / r + p2 * (1 - p2)) / (1 + 1 / r))
    ) ** 2
    ge = lambda a, b: a - b > -TOL
    lt = lambda a, b: b - a > TOL
    n2 = max(1, math.ceil(init))
    P = pa(n2)
    if ge(P, tp):
        while ge(P, tp) and n2 > 1:
            n2 -= 1
            P = pa(n2)
        if lt(P, tp):
            n2 += 1
    else:
        while lt(P, tp):
            n2 += 1
            P = pa(n2)
    return n2


def normal_n2(p1, p2, r, alpha, tp, alt, method):
    ae = alpha / 2 if alt == "two.sided" else alpha
    za, zb = norm.ppf(1 - ae), norm.ppf(tp)
    p = (r * p1 + p2) / (1 + r)
    v0 = max(p * (1 - p), 0) * (1 + 1 / r)
    v1 = (max(p1 * (1 - p1), 0) / r + max(p2 * (1 - p2), 0)) if method == "standard" else v0
    return (za * math.sqrt(v0) + zb * math.sqrt(v1)) ** 2 / (p1 - p2) ** 2


# ---------------------------------------------------------------------------------------
# Re-estimation designs: final group sizes for each pooled interim count s
def final_sizes_rd(DA, r, n11, n12, alpha, tp, test, alt, method):
    n1 = n11 + n12
    out = {}
    for s in range(n1 + 1):
        ph = s / n1
        a, b = ph + DA / (1 + r), ph - r * DA / (1 + r)
        if method == "exact":
            n2 = exact_n2(min(1, a), max(0, b), r, alpha, tp, test, alt)
        else:
            n2 = ceil_tol(normal_n2(a, b, r, alpha, tp, alt, method))
        N2 = max(n12, n2)
        out[s] = (max(n11, math.ceil(r * N2)), N2)
    return out


def bssr_reject_prob(sizes, n11, n12, p1, p2, test, alt, alpha):
    """Sum over every interim cell and every second-stage outcome."""
    w1 = binom.pmf(np.arange(n11 + 1), n11, p1)
    w2 = binom.pmf(np.arange(n12 + 1), n12, p2)
    tot = 0.0
    for x11 in range(n11 + 1):
        for x12 in range(n12 + 1):
            N1, N2 = sizes[x11 + x12]
            m1, m2 = N1 - n11, N2 - n12
            R = reject(N1, N2, test, alt, alpha)[x11:x11 + m1 + 1, x12:x12 + m2 + 1]
            cp = binom.pmf(np.arange(m1 + 1), m1, p1) @ R.astype(float) @ binom.pmf(
                np.arange(m2 + 1), m2, p2)
            tot += w1[x11] * w2[x12] * cp
    return tot


def interim_total_pmf(n11, n12, p1, p2):
    return np.convolve(binom.pmf(np.arange(n11 + 1), n11, p1),
                       binom.pmf(np.arange(n12 + 1), n12, p2))


def expected_n(sizes, n11, n12, p1, p2):
    ps = interim_total_pmf(n11, n12, p1, p2)
    return float(sum(ps[s] * (N1 + N2) for s, (N1, N2) in sizes.items()))


def refined_max(f, grid):
    v = np.array([f(t) for t in grid])
    i = int(v.argmax())
    lo, hi = grid[max(0, i - 1)], grid[min(len(grid) - 1, i + 1)]
    o = minimize_scalar(lambda t: -f(t), bounds=(lo, hi), method="bounded",
                        options={"xatol": 1e-10})
    return (-o.fun, o.x) if -o.fun > v[i] else (v[i], grid[i])


# ---------------------------------------------------------------------------------------
def block_split_pooled():
    f = "test-split_pooled.R"

    def split_or(p, psi, r):
        g = lambda x: (r * (psi * x / (1 - x + psi * x)) + x) / (1 + r) - p
        x = brentq(g, 0, 1, xtol=1e-15, rtol=1e-15)
        return psi * x / (1 - x + psi * x), x

    for (p, psi, r) in [(0.3, 2.0, 1), (0.42, 0.5, 2), (0.8, 3.0, 2)]:
        emit(f, f"OR p={p} psi={psi} r={r}", split_or(p, psi, r))


def block_ss_raw_n2():
    f = "test-ss_raw_n2.R"
    for args in [(0.42, 0.27, 1, 0.025, 0.8, "greater", "standard"),
                 (0.25, 0.05, 1, 0.05, 0.8, "two.sided", "null.variance"),
                 (0.95, 0.75, 3, 0.05, 0.8, "two.sided", "null.variance"),
                 (0.6, 0.3, 2, 0.025, 0.9, "greater", "standard")]:
        emit(f, "ss_raw_n2 " + " ".join(map(str, args)), normal_n2(*args))


def block_bssr_exact():
    f = "test-binary-power-bssr.R"
    cases = [("Chisq", "greater", 1), ("Chisq", "two.sided", 2), ("Boschloo", "greater", 2),
             ("Z-pool", "greater", 1), ("Fisher", "greater", 1)]
    for test, alt, r in cases:
        n12 = math.ceil(0.5 * 10)
        n11 = math.ceil(r * n12)
        sizes = final_sizes_rd(0.3, r, n11, n12, 0.025, 0.8, test, alt, "exact")
        for DT in [0.3, 0.0]:
            pw, en = [], []
            for p in [0.3, 0.45]:
                p1, p2 = p + DT / (1 + r), p - r * DT / (1 + r)
                pw.append(bssr_reject_prob(sizes, n11, n12, p1, p2, test, alt, 0.025))
                en.append(expected_n(sizes, n11, n12, p1, p2))
            emit(f, f"{test} {alt} r={r} Delta.T={DT} power", pw)
            emit(f, f"{test} {alt} r={r} Delta.T={DT} E.N", en)


def block_type1():
    # Delta.A = 0.3, N1 = N2 = 39, interim 20 + 20, chi-squared, standard formula
    n11 = n12 = 20
    sizes = final_sizes_rd(0.3, 1, n11, n12, 0.025, 0.8, "Chisq", "greater", "standard")

    @lru_cache(maxsize=None)
    def boundary(N1, N2, a):
        # Smallest rejected count of group 1 in each column (one-sided regions are monotone)
        R = reject(N1, N2, "Chisq", "greater", a)
        return np.array([np.argmax(R[:, j]) if R[:, j].any() else N1 + 1
                         for j in range(N2 + 1)])

    def tie(t, a=0.025):
        tot = 0.0
        w = binom.pmf(np.arange(n11 + 1), n11, t)
        for s in range(n11 + n12 + 1):
            N1, N2 = sizes[s]
            m1, m2 = N1 - n11, N2 - n12
            B = boundary(N1, N2, a)
            x11 = np.arange(max(0, s - n12), min(n11, s) + 1)
            x12 = s - x11
            b = np.arange(m2 + 1)
            need = B[x12[:, None] + b[None, :]] - x11[:, None]
            tot += np.sum(w[x11] * w[x12] * (binom.sf(need - 1, m1, t) @ binom.pmf(b, m2, t)))
        return tot

    def tfix(t, a=0.025):
        B = boundary(39, 39, a)
        return float(np.sum(binom.pmf(np.arange(40), 39, t) * binom.sf(B - 1, 39, t)))

    g = np.round(np.arange(0.05, 0.951, 0.05), 2)
    f = "test-BinaryTypeIErrorBSSR.R"
    emit(f, "grid BSSR", [tie(t) for t in g])
    emit(f, "grid TRAD", [tfix(t) for t in g])
    G = np.round(np.arange(0.005, 0.996, 0.005), 3)
    mb, mf = refined_max(tie, G), refined_max(tfix, G)
    emit(f, "refined max", [mb[0], mf[0]])
    # The location is determined only to the tolerance of the optimizer
    emit(f, "refined location folded", [min(mb[1], 1 - mb[1]), min(mf[1], 1 - mf[1])],
         rtol=1e-6)

    def adjusted(fun):
        lo, hi = 0.0, 0.025
        while hi - lo > 1e-9 * 0.025:
            mid = (lo + hi) / 2
            if refined_max(lambda t: fun(t, mid), G)[0] <= 0.025:
                lo = mid
            else:
                hi = mid
        return lo

    ab, af = adjusted(tie), adjusted(tfix)
    f = "test-BinaryAlphaAdjBSSR.R"
    emit(f, "max at nominal", [mb[0], mf[0]])
    emit(f, "adjusted level", [ab, af])
    emit(f, "max at adjusted level",
         [refined_max(lambda t: tie(t, ab), G)[0], refined_max(lambda t: tfix(t, af), G)[0]])

    # Distribution of the final total for p = 0.4 and Delta.T = 0.3
    ps = interim_total_pmf(n11, n12, 0.55, 0.25)
    N = np.array([sum(sizes[s]) for s in range(n11 + n12 + 1)])
    E = float(ps @ N)
    order = np.argsort(N, kind="stable")
    cdf = np.cumsum(ps[order])
    q = lambda a: N[order][np.searchsorted(cdf, a - 1e-12)]
    f = "test-summary.bbssr_powerbssr.R"
    emit(f, "E.N", E)
    emit(f, "SD.N", math.sqrt(float(ps @ (N - E) ** 2)))
    emit(f, "quartiles", [q(0.25), q(0.5), q(0.75)])


# ---------------------------------------------------------------------------------------
# Refined maximization over the nuisance parameter (ref.pvalue = TRUE)
def tail_coefficients(stat, N1, N2, decreasing):
    """Bernstein coefficients of the tail probability of every cell, in row-major order.

    The tail probability of a cell is sum_s h[s] b(s; N, theta), where h[s] is the
    hypergeometric probability of the tail set given s responders in total."""
    N = N1 + N2
    s = stat.ravel()
    x1 = np.repeat(np.arange(N1 + 1), N2 + 1)
    x2 = np.tile(np.arange(N2 + 1), N1 + 1)
    order = np.argsort(-s if decreasing else s, kind="stable")
    so = s[order]
    ref = np.maximum(np.abs(so[:-1]), np.abs(so[1:]))
    new = np.abs(np.diff(so)) > 1e-10 * ref
    ends = np.r_[np.where(new)[0], len(so) - 1]
    last = np.repeat(ends, np.diff(np.r_[-1, ends]))
    tot = x1[order] + x2[order]
    H = np.zeros((len(s), N + 1))
    H[np.arange(len(s)), tot] = hypergeom.pmf(x1[order], N, N1, tot)
    H = np.cumsum(H, axis=0)[last]
    out = np.empty_like(H)
    out[order] = H
    return out, x1 + x2


def de_casteljau(c, t):
    """Bernstein coefficients of the same polynomial on [0, t] and on [t, 1]."""
    n = len(c)
    left, right, w = np.empty(n), np.empty(n), c.copy()
    for j in range(n):
        left[j], right[n - 1 - j] = w[0], w[-1]
        w = (1 - t) * w[:-1] + t * w[1:]
    return left, right


def certified_max(h, a=0.0, b=1.0, tol=1e-13):
    """Maximum of sum_s h[s] b(s; N, theta) over [a, b] by branch and bound.

    The largest Bernstein coefficient on an interval bounds the polynomial there, and the
    end coefficients are values of the polynomial, so the interval is split until every
    bound is within tol of the best value found. Returns the best value and the largest
    bound left."""
    c = h
    if b < 1:
        c, _ = de_casteljau(c, b)
    if a > 0:
        _, c = de_casteljau(c, a / b)
    best = max(c[0], c[-1])
    stack, bounds = [c], [best]
    while stack:
        c = stack.pop()
        if c.max() <= best + tol or len(bounds) > 10 ** 6:
            bounds.append(c.max())
            continue
        left, right = de_casteljau(c, 0.5)
        best = max(best, left[-1])
        stack += [left, right]
    return best, max(bounds)


def certified_pvalues(stat, N1, N2, decreasing, gamma=0.0):
    H, tot = tail_coefficients(stat, N1, N2, decreasing)
    N = N1 + N2
    if gamma > 0:
        s = np.arange(N + 1)
        lo = beta_dist.ppf(gamma / 2, np.maximum(s, 1), N - s + 1)
        up = beta_dist.ppf(1 - gamma / 2, s + 1, np.maximum(N - s, 1))
        lo[s == 0], up[s == N] = 0.0, 1.0
    p, gap, memo = np.empty(len(H)), 0.0, {}
    for k in range(len(H)):
        a, b = (lo[tot[k]], up[tot[k]]) if gamma > 0 else (0.0, 1.0)
        key = (H[k].tobytes(), a, b)
        if key not in memo:
            memo[key] = certified_max(H[k], a, b)
        best, bound = memo[key]
        p[k], gap = best, max(gap, bound - best)
    assert gap < 1e-12, gap
    return np.minimum(1, p + gamma).reshape(N1 + 1, N2 + 1)


def certified_size(R, N1, N2):
    N = N1 + N2
    h = np.zeros(N + 1)
    for i, j in zip(*np.nonzero(R)):
        h[i + j] += hypergeom.pmf(i, N, N1, i + j)
    best, bound = certified_max(h)
    assert bound - best < 1e-12
    return best


def block_refined():
    f = "test-max_tail_prob_refined.R"
    N1 = N2 = 32
    for test, stat, dec in [("Z-pool", zstat(N1, N2), True), ("Boschloo", fisher_upper(N1, N2), False)]:
        cert = certified_pvalues(stat, N1, N2, dec)
        grid = pvalues(N1, N2, test, "greater")
        emit(f, f"{test} 32 x 32 sum", cert.sum())
        # Cells whose decision at the level 0.025 changes with the refinement
        flip = np.argwhere((grid < 0.025 - TOL) & ~(cert < 0.025 - TOL))
        emit(f, f"{test} 32 x 32 flipped cells", flip.ravel(), rtol=0)
        emit(f, f"{test} 32 x 32 flipped p-values", cert[tuple(flip.T)])
        f2 = "test-binary-rr.R"
        emit(f2, f"{test} 32 x 32 rejected", [(grid < 0.025 - TOL).sum(), (cert < 0.025 - TOL).sum()], rtol=0)
        size = [certified_size(grid < 0.025 - TOL, N1, N2), certified_size(cert < 0.025 - TOL, N1, N2)]
        emit(f2, f"{test} 32 x 32 size", size)
        if test == "Z-pool":
            emit("test-BinaryTypeIErrorBSSR.R", "Z-pool 32 x 32 fixed-design size", size)
    # Berger-Boos interval
    cert = certified_pvalues(fisher_upper(20, 15), 20, 15, False, gamma=0.001)
    emit(f, "Boschloo 20 x 15 Berger-Boos sum", cert.sum())
    # A peak between two points of the uniform grid that the arcsine grid resolves
    H, _ = tail_coefficients(zstat(150, 60), 150, 60, True)
    v = certified_max(H[103 * 61 + 32])[0]
    emit("test-binary-rr.R", "Z-pool 150 x 60 cell (103, 32)", v)
    emit(f, "Z-pool 150 x 60 cell (103, 32)", v)
    # Grids on which a point of the arcsine grid falls within rounding of a point of the
    # uniform grid (n.grid = 101 contains 0.25, 0.5 and 0.75). The certified maximum does
    # not depend on the grid
    emit("test-binary-rr.R", "Z-pool 6 x 12 sum", certified_pvalues(zstat(6, 12), 6, 12, True).sum())
    emit("test-binary-rr.R", "Boschloo 5 x 89 sum",
         certified_pvalues(fisher_upper(5, 89), 5, 89, False).sum())


if __name__ == "__main__":
    print("# test file\tkey\trelative tolerance\tvalues (15 significant digits)")
    block_split_pooled()
    block_ss_raw_n2()
    block_type1()
    block_bssr_exact()
    block_refined()
