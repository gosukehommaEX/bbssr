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


def fisher_two_sided(N1, N2, ts, midp=False):
    """Two-sided Fisher p-value from the hypergeometric tails of every total s.

    minlike orders the tables by their probability, blaker by the smaller of their two
    tail probabilities (formula (2) of Mehrotra, Chan and Berger, 2003), and central
    doubles the smaller tail. With midp the tables tied with the observed one in the
    ordering (the observed one included) contribute half of their probability."""
    P = np.full((N1 + 1, N2 + 1), np.nan)
    for s in range(N1 + N2 + 1):
        k = np.arange(max(0, s - N2), min(N1, s) + 1)
        f = hypergeom.pmf(k, N1 + N2, N1, s)
        lo = hypergeom.cdf(k, N1 + N2, N1, s)
        up = hypergeom.sf(k - 1, N1 + N2, N1, s)
        if ts == "central":
            if midp:
                lo, up = lo - f / 2, up - f / 2
            p = 2 * np.minimum(lo, up)
        else:
            key = f if ts == "minlike" else np.minimum(lo, up)
            p = np.empty(len(k))
            for i in range(len(k)):
                tied = np.abs(key - key[i]) <= 1e-10 * np.maximum(key, key[i])
                more = (key < key[i]) & ~tied
                p[i] = f[more].sum() + (0.5 if midp else 1.0) * f[tied].sum()
        P[k, s - k] = np.minimum(1, p)
    return P


@lru_cache(maxsize=None)
def pvalues(N1, N2, test, alt, ts=None):
    if alt == "two.sided":
        if test == "Chisq":
            return np.minimum(2 * norm.sf(np.abs(zstat(N1, N2))), 1)
        if test in ("Fisher", "Fisher-midP"):
            return fisher_two_sided(N1, N2, ts, midp=(test == "Fisher-midP"))
        if test == "Boschloo":
            return unconditional(fisher_two_sided(N1, N2, ts), N1, N2, False)
        raise ValueError(test)
    if test == "Chisq":
        return norm.sf(zstat(N1, N2))
    if test == "Fisher":
        return fisher_upper(N1, N2)
    if test == "Z-pool":
        return unconditional(zstat(N1, N2), N1, N2, True)
    if test == "Boschloo":
        return unconditional(fisher_upper(N1, N2), N1, N2, False)
    raise ValueError(test)


def reject(N1, N2, test, alt, alpha, ts=None):
    return pvalues(N1, N2, test, alt, ts) < alpha - TOL


def power(p1, p2, N1, N2, test, alt, alpha, ts=None):
    R = reject(N1, N2, test, alt, alpha, ts).astype(float)
    return binom.pmf(np.arange(N1 + 1), N1, p1) @ R @ binom.pmf(np.arange(N2 + 1), N2, p2)


# ---------------------------------------------------------------------------------------
# Sample size rules
def exact_n2(p1, p2, r, alpha, tp, test, alt, ts=None):
    pa = lambda n2: power(p1, p2, math.ceil(r * n2), n2, test, alt, alpha, ts)
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


# ---------------------------------------------------------------------------------------
# Non-inferiority tests of Blackwelder (1982) and Farrington and Manning (1990)
def restricted_mle(x1, n1, x2, n2, s0):
    """Maximum likelihood estimates under p1 - p2 = s0, by bisection on the score.

    The log likelihood is concave in p1 on the admissible interval, so its derivative is
    decreasing and its sign change is found by bisection. This shares nothing with the
    closed form of the package."""
    x1, x2 = np.broadcast_arrays(np.asarray(x1, float), np.asarray(x2, float))
    lo = np.full(x1.shape, max(0.0, s0)); hi = np.full(x1.shape, min(1.0, 1.0 + s0))
    def score(t):
        t2 = t - s0
        with np.errstate(divide="ignore", invalid="ignore"):
            a = np.where(x1 > 0, x1 / t, 0.0) - np.where(n1 - x1 > 0, (n1 - x1) / (1 - t), 0.0)
            b = np.where(x2 > 0, x2 / t2, 0.0) - np.where(n2 - x2 > 0, (n2 - x2) / (1 - t2), 0.0)
        return a + b
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        sc = score(mid)
        lo = np.where(sc > 0, mid, lo); hi = np.where(sc > 0, hi, mid)
    t = 0.5 * (lo + hi)
    return t, t - s0


def large_sample_restricted(p1, p2, theta, s0):
    """Large sample values: the restricted estimates with expected counts p1 and theta p2."""
    t1, t2 = restricted_mle(np.asarray(p1) * 1e6, 1e6, np.asarray(p2) * theta * 1e6, theta * 1e6, s0)
    return float(t1), float(t2)


@lru_cache(maxsize=None)
def ni_pvalues(N1, N2, test, margin):
    x1 = np.arange(N1 + 1)[:, None]; x2 = np.arange(N2 + 1)[None, :]
    h1, h2 = x1 / N1, x2 / N2
    num = h1 - h2 + margin + 0 * x2
    if test == "Blackwelder":
        var = h1 * (1 - h1) / N1 + h2 * (1 - h2) / N2
    else:
        t1, t2 = restricted_mle(x1 + 0 * x2, N1, x2 + 0 * x1, N2, -margin)
        var = t1 * (1 - t1) / N1 + t2 * (1 - t2) / N2
    var = var + 0 * num
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(var > 0, num / np.sqrt(np.where(var > 0, var, 1.0)),
                     np.where(num == 0, 0.0, np.sign(num) * np.inf))
    return norm.sf(z)


def ni_reject(N1, N2, test, margin, alpha):
    return ni_pvalues(N1, N2, test, margin) < alpha - TOL


def ni_power(p1, p2, N1, N2, test, margin, alpha):
    R = ni_reject(N1, N2, test, margin, alpha).astype(float)
    return binom.pmf(np.arange(N1 + 1), N1, p1) @ R @ binom.pmf(np.arange(N2 + 1), N2, p2)


def ni_raw_n2(p1, p2, r, alpha, tp, margin, method):
    """Unrounded size of group 2, alternative 'greater'."""
    za, zb = norm.ppf(1 - alpha), norm.ppf(tp)
    t1, t2 = large_sample_restricted(min(1, max(0, p1)), min(1, max(0, p2)), 1 / r, -margin)
    v0 = t1 * (1 - t1) / r + t2 * (1 - t2)
    v1 = max(p1 * (1 - p1), 0) / r + max(p2 * (1 - p2), 0)
    if method == "null.variance":
        v1 = v0
    if method == "alternative.variance":
        v0 = v1
    return (za * math.sqrt(v0) + zb * math.sqrt(v1)) ** 2 / (p1 - p2 + margin) ** 2


def ni_exact_n2(p1, p2, r, alpha, tp, test, margin):
    pa = lambda n2: ni_power(p1, p2, math.ceil(r * n2), n2, test, margin, alpha)
    ge = lambda a, b: a - b > -TOL
    lt = lambda a, b: b - a > TOL
    n2 = max(1, math.ceil(ni_raw_n2(p1, p2, r, alpha, tp, margin, "standard")))
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


def ni_bssr(n11, n12, r, DA, margin, alpha, tp, test, method, p1, p2):
    """Unrestricted design, nearest rounding of each group, power and expected size."""
    n = n11 + n12
    sizes = {}
    for s in range(n + 1):
        ph = s / n
        a, b = ph + DA / (1 + r), ph - r * DA / (1 + r)
        n2 = ni_raw_n2(a, b, r, alpha, tp, margin, method)
        tot = max(n, (1 + r) * n2)
        sizes[s] = (max(n11, math.floor(r * tot / (1 + r) + 0.5)), max(n12, math.floor(tot / (1 + r) + 0.5)))
    w1 = binom.pmf(np.arange(n11 + 1), n11, p1); w2 = binom.pmf(np.arange(n12 + 1), n12, p2)
    power = 0.0
    for x11 in range(n11 + 1):
        for x12 in range(n12 + 1):
            N1, N2 = sizes[x11 + x12]; m1, m2 = N1 - n11, N2 - n12
            R = ni_reject(N1, N2, test, margin, alpha)[x11:x11 + m1 + 1, x12:x12 + m2 + 1]
            power += w1[x11] * w2[x12] * (binom.pmf(np.arange(m1 + 1), m1, p1) @ R.astype(float)
                                           @ binom.pmf(np.arange(m2 + 1), m2, p2))
    ps = interim_total_pmf(n11, n12, p1, p2)
    EN = float(ps @ np.array([sum(sizes[s]) for s in range(n + 1)]))
    return power, EN, sizes


def block_ni():
    f = "test-fm_restricted.R"
    cases = [(3, 12, 7, 20, -0.15), (0, 18, 0, 18, -0.2), (12, 12, 20, 20, 0.1), (5, 10, 0, 15, 0.3),
             (8, 25, 9, 30, -0.05)]
    vals = []
    for x1, n1, x2, n2, s0 in cases:
        t1, t2 = restricted_mle(x1, n1, x2, n2, s0)
        vals += [float(t1), float(t2)]
    # Estimates on the boundary of the admissible range are 0 up to the bisection error
    vals = [0.0 if abs(v) < 1e-12 else v for v in vals]
    emit(f, "restricted estimates", vals, rtol=1e-9)
    emit(f, "large sample values", large_sample_restricted(0.7, 0.7, 1 / 3, -0.1), rtol=1e-6)
    f = "test-zstat_margin.R"
    z1 = (13 / 30 - 18 / 30 + 0) / math.sqrt((13 / 30) * (17 / 30) / 30 + 0.6 * 0.4 / 30)
    z2 = (18 / 30 - 21 / 30 + 0.2) / math.sqrt(0.6 * 0.4 / 30 + 0.7 * 0.3 / 30)
    emit(f, "Blackwelder 1982 examples", [-z1, z2])
    t1, t2 = restricted_mle(9, 25, 14, 30, -0.1)
    emit(f, "Farrington-Manning cell", (9 / 25 - 14 / 30 + 0.1) / math.sqrt(t1 * (1 - t1) / 25 + t2 * (1 - t2) / 30),
         rtol=1e-8)
    f = "test-binary-power.R"
    emit(f, "Farrington-Manning 1990 example", ni_power(0.4, 0.05, 80, 80, "FM", -0.2, 0.05), rtol=1e-8)
    emit(f, "Farrington-Manning 1990 Table I", [ni_power(0.2, 0.1, 57, 57, "FM", 0.1, 0.05),
                                                ni_power(0.5, 0.1, 67, 101, "FM", -0.2, 0.05),
                                                ni_power(0.1, 0.05, 168, 112, "FM", 0.05, 0.05)], rtol=1e-8)
    emit(f, "Blackwelder power", ni_power(0.65, 0.7, 120, 100, "Blackwelder", 0.15, 0.025), rtol=1e-8)
    f = "test-ss_raw_n2.R"
    emit(f, "Farrington-Manning raw n2", [ni_raw_n2(0.7, 0.7, 3, 0.025, 0.8, 0.1, "standard"),
                                          ni_raw_n2(0.3, 0.25, 0.5, 0.05, 0.9, 0.1, "null.variance")],
         rtol=1e-8)
    emit(f, "Blackwelder raw n2", ni_raw_n2(0.9, 0.9, 1, 0.05, 0.9, 0.1, "alternative.variance"), rtol=1e-10)
    f = "test-binary-sample-size.R"
    emit(f, "exact Farrington-Manning N2", [ni_exact_n2(0.8, 0.8, 1, 0.025, 0.8, "FM", 0.15),
                                            ni_exact_n2(0.6, 0.65, 2, 0.025, 0.9, "Blackwelder", 0.2)],
         rtol=0)
    f = "test-binary-power-bssr.R"
    pw, EN, _ = ni_bssr(30, 30, 1, 0, 0.2, 0.025, 0.8, "FM", "standard", 0.4, 0.4)
    pw0, EN0, _ = ni_bssr(30, 30, 1, 0, 0.2, 0.025, 0.8, "FM", "standard", 0.4, 0.6)
    emit(f, "non-inferiority power and E.N", [pw, EN])
    emit(f, "non-inferiority boundary", [pw0, EN0])
    f = "test-BinaryTypeIErrorBSSR.R"
    tie = []
    for th in (0.3, 0.5, 0.7):
        a, b = th - 0.2 / 2, th + 0.2 / 2
        tie.append(ni_bssr(30, 30, 1, 0, 0.2, 0.025, 0.8, "FM", "standard", a, b)[0])
        tie.append(ni_power(a, b, 54, 54, "FM", 0.2, 0.025))
    emit(f, "non-inferiority type I error", tie)


# ---------------------------------------------------------------------------------------
# Conditional rejection probabilities given the pooled responder counts (B3)
def hyper_exact(x, M, n, N):
    """Hypergeometric probabilities from exact integer binomial coefficients."""
    d = math.comb(M, N)
    return np.array([math.comb(n, int(k)) * math.comb(M - n, N - int(k)) / d
                     for k in np.atleast_1d(x)])


def crp_tables(sizes, n11, n12, test, alt, alpha):
    """Rows (s, s2, CRP, CRP.total) ordered by s and s2, by direct summation."""
    rows = []
    n1 = n11 + n12
    for s in range(n1 + 1):
        N1, N2 = sizes[s]
        m1, m2 = N1 - n11, N2 - n12
        R = reject(N1, N2, test, alt, alpha)
        x11 = np.arange(max(0, s - n12), min(n11, s) + 1)
        h1 = hyper_exact(x11, n1, n11, s)
        for s2 in range(m1 + m2 + 1):
            x21 = np.arange(max(0, s2 - m2), min(m1, s2) + 1)
            h2 = hyper_exact(x21, m1 + m2, m1, s2)
            X1 = x11[:, None] + x21[None, :]
            c = float(np.sum(h1[:, None] * h2[None, :] * R[X1, s + s2 - X1]))
            t = s + s2
            x1 = np.arange(max(0, t - N2), min(N1, t) + 1)
            ct = float(np.sum(hyper_exact(x1, N1 + N2, N1, t) * R[x1, t - x1]))
            rows.append((s, s2, c, ct))
    return rows


def block_crp():
    f = "test-BinaryCondRejectBSSR.R"
    # Fisher's exact test, Delta.A = 0.3, N1 = N2 = 12, interim 6 + 6, standard formula
    sizes = final_sizes_rd(0.3, 1, 6, 6, 0.025, 0.8, "Fisher", "greater", "standard")
    rows = crp_tables(sizes, 6, 6, "Fisher", "greater", 0.025)
    c = np.array([x[2] for x in rows]); ct = np.array([x[3] for x in rows])
    above = [x for x in rows if x[2] > 0.025 + TOL]
    emit(f, "Fisher 6 + 6 rows", len(rows), rtol=0)
    emit(f, "Fisher 6 + 6 sums", [c.sum(), ct.sum()])
    emit(f, "Fisher 6 + 6 cells above the level s", [x[0] for x in above], rtol=0)
    emit(f, "Fisher 6 + 6 cells above the level s2", [x[1] for x in above], rtol=0)
    emit(f, "Fisher 6 + 6 values above the level", [x[2] for x in above])
    emit(f, "Fisher 6 + 6 maxima", [c.max(), ct.max()])
    emit(f, "Fisher 6 + 6 counts above the level", [len(above), int((ct > 0.025 + TOL).sum())], rtol=0)

    def tie(t):
        tot = 0.0
        for (s, s2, cc, _) in rows:
            n2 = sum(sizes[s]) - 12
            tot += binom.pmf(s, 12, t) * binom.pmf(s2, n2, t) * cc
        return tot
    emit(f, "Fisher 6 + 6 type I error at 0.2 and 0.5", [tie(0.2), tie(0.5)])
    # The same design without re-estimation: CRP differs from CRP.total
    fixed = crp_tables({s: (12, 12) for s in range(13)}, 6, 6, "Fisher", "greater", 0.025)
    emit(f, "Fisher 12 + 12 fixed largest difference",
         max(abs(x[2] - x[3]) for x in fixed), rtol=1e-8)
    # Boschloo test, r = 2, Delta.A = 0.3, N2 = 10, interim 10 + 5, standard formula
    sizes = final_sizes_rd(0.3, 2, 10, 5, 0.025, 0.8, "Boschloo", "greater", "standard")
    rows = crp_tables(sizes, 10, 5, "Boschloo", "greater", 0.025)
    c = np.array([x[2] for x in rows]); ct = np.array([x[3] for x in rows])
    i, j = int(c.argmax()), int(ct.argmax())
    emit(f, "Boschloo r = 2 rows and locations of the maxima",
         [len(rows), rows[i][0], rows[i][1], rows[j][0], rows[j][1]], rtol=0)
    emit(f, "Boschloo r = 2 sums", [c.sum(), ct.sum()])
    emit(f, "Boschloo r = 2 maxima", [c.max(), ct.max()])
    emit(f, "Boschloo r = 2 counts above the level",
         [int((c > 0.025 + TOL).sum()), int((ct > 0.025 + TOL).sum())], rtol=0)


# ---------------------------------------------------------------------------------------
# Two-sided conventions of the conditional tests, including Blaker's
def block_blaker():
    f = "test-internal.R"
    vals = [fisher_two_sided(14, 7, ts, midp)[8, 1]
            for ts, midp in [("blaker", False), ("minlike", False), ("central", False),
                             ("blaker", True), ("minlike", True)]]
    emit(f, "Fisher two-sided 14 x 7 at (8, 1)", vals)
    f = "test-binary-rr.R"
    stat = fisher_two_sided(14, 7, "blaker")
    emit(f, "Boschloo blaker 14 x 7 grid sum", unconditional(stat, 14, 7, False).sum())
    cert = certified_pvalues(stat, 14, 7, False)
    emit(f, "Boschloo blaker 14 x 7 certified sum", cert.sum())
    emit(f, "blaker 14 x 7 rejected at 0.05",
         [(stat < 0.05 - TOL).sum(), (cert < 0.05 - TOL).sum()], rtol=0)
    # A configuration of Table 3 of Mehrotra, Chan and Berger (2003) in which the three
    # conventions give different powers
    f = "test-binary-power.R"
    emit(f, "Fisher two-sided 10 x 40 power at (0.5, 0.86)",
         [power(0.5, 0.86, 10, 40, "Fisher", "two.sided", 0.05, ts)
          for ts in ["blaker", "minlike", "central"]])
    f = "test-binary-sample-size.R"
    emit(f, "Fisher two-sided r = 4 N2 blaker minlike central",
         [exact_n2(0.1, 0.5, 4, 0.05, 0.8, "Fisher", "two.sided", ts)
          for ts in ["blaker", "minlike", "central"]], rtol=0)


# ---------------------------------------------------------------------------------------
# Grid of designs (B4)
def block_grid():
    f = "test-BinaryGridBSSR.R"
    # Two designs of test-binary-power-bssr.R evaluated through the grid: Delta.A = 0.3,
    # N1 = N2 = 10, omega = 0.5, exact re-estimation, p = 0.3 and 0.45, Delta.T = 0.3
    for test in ["Chisq", "Fisher"]:
        sizes = final_sizes_rd(0.3, 1, 5, 5, 0.025, 0.8, test, "greater", "exact")
        pw, en, sd = [], [], []
        for p in [0.3, 0.45]:
            p1, p2 = p + 0.15, p - 0.15
            pw.append(bssr_reject_prob(sizes, 5, 5, p1, p2, test, "greater", 0.025))
            ps = interim_total_pmf(5, 5, p1, p2)
            N = np.array([sum(sizes[s]) for s in range(11)])
            e = float(ps @ N)
            en.append(e)
            sd.append(math.sqrt(float(ps @ (N - e) ** 2)))
        emit(f, f"{test} grid power", pw)
        emit(f, f"{test} grid E.N", en)
        emit(f, f"{test} grid SD.N", sd)
    # Initial sample sizes from a planning proportion of 0.4 and Delta.A = 0.3
    n2 = exact_n2(0.55, 0.25, 1, 0.025, 0.8, "Chisq", "greater")
    emit(f, "p.plan exact Chisq r = 1 N1 N2", [math.ceil(n2), n2], rtol=0)
    emit(f, "p.plan exact Chisq r = 1 fixed power at p = 0.4",
         power(0.55, 0.25, math.ceil(n2), n2, "Chisq", "greater", 0.025))
    n2s = ceil_tol(normal_n2(0.5, 0.2, 2, 0.025, 0.8, "greater", "standard"))
    emit(f, "p.plan standard r = 2 N1 N2", [math.ceil(2 * n2s), n2s], rtol=0)


# ---------------------------------------------------------------------------------------
# Certified maxima of the type I error rate over an interval of theta (item C)
def interim_hyper(sizes, n11, n12):
    """For each pair of final sizes (N1, N2), the matrix H with H[X1, X2] the sum, over the
    interim cells (x11, x12) leading to (N1, N2), of P(x11 | X1) P(x12 | X2). Given the
    final responder count X1 of N1 patients, the count among the first n11 is
    hypergeometric, so the rejection probability of the design is
    sum over (N1, N2) of sum_{X1, X2} R[X1, X2] H[X1, X2] b(X1; N1, p1) b(X2; N2, p2)."""
    groups = {}
    for s, NN in sizes.items():
        groups.setdefault(NN, []).append(s)
    out = []
    for (N1, N2), ss in sorted(groups.items()):
        H = np.zeros((N1 + 1, N2 + 1))
        X1, X2 = np.arange(N1 + 1), np.arange(N2 + 1)
        for s in ss:
            for x11 in range(max(0, s - n12), min(n11, s) + 1):
                H += np.outer(hypergeom.pmf(x11, N1, n11, X1), hypergeom.pmf(s - x11, N2, n12, X2))
        out.append((N1, N2, H))
    return out


def tie_coefficients(hyper, rej, a1=0.0, c1=1.0, a2=0.0, c2=1.0):
    """Bernstein coefficients in t of the rejection probability when p1 = a1 (1 - t) + c1 t
    and p2 = a2 (1 - t) + c2 t. b(X; N, a (1 - t) + c t) = sum_j M[X, j] B_{j,N}(t), where
    column j of M is the distribution of Bin(N - j, a) + Bin(j, c); the product of the
    bases of degrees N1 and N2 is a basis of degree N1 + N2 with the factor
    hyper(j; N1 + N2, N1, j + l); every degree is raised to the largest by the matrix of
    hypergeometric probabilities."""
    def conv(N, a, c):
        M = np.zeros((N + 1, N + 1))
        for j in range(N + 1):
            M[:, j] = np.convolve(binom.pmf(np.arange(N - j + 1), N - j, a),
                                  binom.pmf(np.arange(j + 1), j, c))
        return M
    Nmax = max(N1 + N2 for N1, N2, _ in hyper)
    b = np.zeros(Nmax + 1)
    for N1, N2, H in hyper:
        D = conv(N1, a1, c1).T @ (rej(N1, N2) * H) @ conv(N2, a2, c2)
        N = N1 + N2
        h = np.zeros(N + 1)
        for j in range(N1 + 1):
            h[j:j + N2 + 1] += D[j] * hypergeom.pmf(j, N, N1, j + np.arange(N2 + 1))
        E = hypergeom.pmf(np.arange(N + 1)[None, :], Nmax, N, np.arange(Nmax + 1)[:, None])
        b += E @ h
    return b


def bernstein_value(b, t):
    return float(binom.pmf(np.arange(len(b)), len(b) - 1, t) @ b)


def checked_max(b, a=0.0, c=1.0):
    """Certified maximum over [a, c], compared with a grid of 2001 points refined by a
    bounded optimization around its three largest local maxima."""
    best, bound = certified_max(b, a, c)
    assert 0 <= bound - best < 1e-12, (best, bound)
    g = np.linspace(a, c, 2001)
    v = np.array([bernstein_value(b, t) for t in g])
    pk = [i for i in range(len(g)) if (i == 0 or v[i] >= v[i - 1]) and (i == len(g) - 1 or v[i] >= v[i + 1])]
    ref = v.max()
    for i in sorted(pk, key=lambda i: -v[i])[:3]:
        o = minimize_scalar(lambda t: -bernstein_value(b, t), method="bounded",
                            bounds=(g[max(0, i - 1)], g[min(len(g) - 1, i + 1)]),
                            options={"xatol": 1e-12})
        ref = max(ref, -o.fun)
    assert abs(best - ref) < 1e-12, (best, ref)
    return best


def block_certified():
    f = "test-BinaryTypeIErrorBSSR.R"
    # Design of block_type1: Delta.A = 0.3, N1 = N2 = 39, interim 20 + 20, chi-squared
    sizes = final_sizes_rd(0.3, 1, 20, 20, 0.025, 0.8, "Chisq", "greater", "standard")
    hyp = {"BSSR": interim_hyper(sizes, 20, 20), "TRAD": interim_hyper({0: (39, 39)}, 0, 0)}
    chisq = lambda a: (lambda N1, N2: reject(N1, N2, "Chisq", "greater", a))
    coef = {d: tie_coefficients(hyp[d], chisq(0.025)) for d in hyp}
    # The coefficients reproduce the direct sums at two values of theta
    for t in (0.2, 0.47):
        assert abs(bernstein_value(coef["BSSR"], t)
                   - bssr_reject_prob(sizes, 20, 20, t, t, "Chisq", "greater", 0.025)) < 1e-14
        assert abs(bernstein_value(coef["TRAD"], t) - power(t, t, 39, 39, "Chisq", "greater", 0.025)) < 1e-14
    emit(f, "certified max", [checked_max(coef[d]) for d in ("BSSR", "TRAD")])
    # On a sub-interval, by de Casteljau subdivision of the coefficients on [0, 1]
    emit(f, "certified max on [0.15, 0.85]", [checked_max(coef[d], 0.15, 0.85) for d in ("BSSR", "TRAD")])
    # Boundary p1 = theta - 0.1, p2 = theta + 0.1 of the non-inferiority design of block_ni,
    # over its whole range [0.1, 0.9]: t = (theta - 0.1) / 0.8, p1 from 0 to 0.8 and p2
    # from 0.2 to 1
    _, _, ni_sizes = ni_bssr(30, 30, 1, 0, 0.2, 0.025, 0.8, "FM", "standard", 0.4, 0.4)
    fm = lambda N1, N2: ni_reject(N1, N2, "FM", 0.2, 0.025)
    ni = {"BSSR": tie_coefficients(interim_hyper(ni_sizes, 30, 30), fm, 0.0, 0.8, 0.2, 1.0),
          "TRAD": tie_coefficients(interim_hyper({0: (54, 54)}, 0, 0), fm, 0.0, 0.8, 0.2, 1.0)}
    for th in (0.3, 0.5, 0.7):
        t = (th - 0.1) / 0.8
        a, b = th - 0.1, th + 0.1
        assert abs(bernstein_value(ni["BSSR"], t)
                   - ni_bssr(30, 30, 1, 0, 0.2, 0.025, 0.8, "FM", "standard", a, b)[0]) < 1e-14
        assert abs(bernstein_value(ni["TRAD"], t) - ni_power(a, b, 54, 54, "FM", 0.2, 0.025)) < 1e-14
    emit(f, "non-inferiority certified max", [checked_max(ni[d]) for d in ("BSSR", "TRAD")])

    # Adjusted levels with every decision certified over [0.1, 0.9]
    f = "test-BinaryAlphaAdjBSSR.R"
    lev, mx = [], []
    for d in ("BSSR", "TRAD"):
        cert = lambda a: certified_max(tie_coefficients(hyp[d], chisq(a)), 0.1, 0.9)
        lo, hi = 0.0, 0.025
        m = cert(0.025)
        if m[1] > 0.025:
            while hi - lo > 1e-10 * 0.025:
                mid = (lo + hi) / 2
                if cert(mid)[1] <= 0.025:
                    lo = mid
                else:
                    hi = mid
            m = cert(lo)
        else:
            lo = 0.025
        lev.append(lo)
        mx.append(m[0])
    emit(f, "adjusted level on [0.1, 0.9]", lev, rtol=1e-7)
    emit(f, "max at adjusted level on [0.1, 0.9]", mx, rtol=1e-7)


if __name__ == "__main__":
    print("# test file\tkey\trelative tolerance\tvalues (15 significant digits)")
    block_split_pooled()
    block_ss_raw_n2()
    block_type1()
    block_bssr_exact()
    block_refined()
    block_ni()
    block_crp()
    block_grid()
    block_blaker()
    block_certified()
