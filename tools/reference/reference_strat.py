"""Independent reference values for the stratified functions of bbssr.

Called from reference_values.py, which passes its emit() and reported_level(). The
computations follow the definitions directly and share no code with the package: the
conditional distribution of the stratified exact test is obtained by enumerating every
combination of the numbers of responders of group 1 in the strata, the power of a fixed
design by summing the binomial probabilities of every table that the p-value rejects, the
power of a re-estimation design by enumerating the interim and second-stage numbers of
responders of both groups in every stratum, and the largest type I error rate by a grid
with local optimization. The power of a design with random stratum and group sizes sums
the power of the fixed design of every configuration of sizes, without leaving any out.
"""
import itertools
import math
from fractions import Fraction

import numpy as np
from scipy.optimize import brentq, minimize
from scipy.special import gammaln
from scipy.stats import binom, hypergeom, norm

TOL = 1.4901161193847656e-08  # tolerance of fpCompare


def lchoose(n, k):
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


# ---------------------------------------------------------------------------------------
# p-values
_NULL = {}


def null_dist(N1, N2, z):
    """Distribution of S = sum_j x_j given the stratum totals z: the convolution of the
    hypergeometric distributions of the strata."""
    key = (tuple(N1), tuple(N2), tuple(z))
    if key in _NULL:
        return _NULL[key]
    lo = 0
    f = np.array([1.0])
    for j in range(len(N1)):
        a = max(0, z[j] - N2[j])
        b = min(z[j], N1[j])
        x = np.arange(a, b + 1)
        g = hypergeom.pmf(x, N1[j] + N2[j], N1[j], z[j])
        f = np.convolve(f, g)
        lo += a
    S = np.arange(lo, lo + len(f))
    _NULL[key] = (S, f)
    return S, f


def fisher_p(N1, N2, z, S, alt, ts):
    Ss, f = null_dist(N1, N2, z)
    up = f[Ss >= S].sum()
    lo = f[Ss <= S].sum()
    if alt == "greater":
        return min(1.0, up)
    if alt == "less":
        return min(1.0, lo)
    if ts == "central":
        return min(1.0, 2 * min(lo, up))
    if ts == "minlike":
        key = f
    else:
        cl = np.cumsum(f)
        cu = np.cumsum(f[::-1])[::-1]
        key = np.minimum(cl, cu)
    k0 = key[Ss == S][0]
    return min(1.0, f[key <= k0 * (1 + 1e-10)].sum())


def asym_p(N1, N2, z, S, test, alt):
    E = V = 0.0
    for j in range(len(N1)):
        n = N1[j] + N2[j]
        E += N1[j] * z[j] / n
        # A stratum of one patient has an empty group and adds nothing to V
        if n < 2:
            continue
        if test.startswith("Cochran"):
            V += N1[j] * N2[j] * z[j] * (n - z[j]) / n ** 3
        else:
            V += N1[j] * N2[j] * z[j] * (n - z[j]) / (n ** 2 * (n - 1))
    if V <= 0:
        return 1.0
    d = S - E
    y = 0.5 if (test.endswith("-cc") and abs(d) >= 0.5) else 0.0
    if alt == "two.sided":
        return min(1.0, 2 * norm.sf((abs(d) - y) / math.sqrt(V)))
    Z = math.copysign(1.0, d) * (abs(d) - y) / math.sqrt(V) if d != 0 else 0.0
    return norm.sf(Z) if alt == "greater" else norm.cdf(Z)


UNSTRAT = ("Unstrat-Chisq", "Unstrat-Chisq-cc", "Unstrat-Fisher")


def pearson_p(n1, n2, X1, Z, alt, cc):
    """Comparison of two proportions with the pooled variance on the table of totals,
    X1 of n1 and Z - X1 of n2 responders, written on the scale of the count X1: the
    deviation d = X1 - n1 Z / n (exact), reduced by 1/2 when |d| >= 1/2 with cc, over
    sqrt(pbar qbar n1 n2 / n)."""
    n = n1 + n2
    if Z == 0 or Z == n or n1 == 0 or n2 == 0:
        return 1.0
    d = Fraction(X1 * n - n1 * Z, n)
    y = Fraction(1, 2) if (cc and abs(d) >= Fraction(1, 2)) else Fraction(0)
    pbar = Z / n
    sd = math.sqrt(pbar * (1 - pbar) * n1 * n2 / n)
    zv = float(abs(d) - y) / sd
    if alt == "two.sided":
        return min(1.0, 2 * norm.sf(zv))
    zs = zv if d > 0 else (-zv if d < 0 else 0.0)
    return norm.sf(zs) if alt == "greater" else norm.cdf(zs)


def p_of(N1, N2, z, S, test, alt, ts):
    """p-value given the stratum totals z and S = sum of the group 1 responders. The
    unstratified tests use the table of totals only."""
    if test in UNSTRAT:
        n1, n2, Z = sum(N1), sum(N2), sum(z)
        if test == "Unstrat-Fisher":
            return fisher_p([n1], [n2], [Z], S, alt, ts)
        return pearson_p(n1, n2, S, Z, alt, test.endswith("-cc"))
    if test == "Strat-Fisher":
        return fisher_p(N1, N2, z, S, alt, ts)
    return asym_p(N1, N2, z, S, test, alt)


def pvalue(test, alt, ts, N1, N2, x1, x2):
    z = [x1[j] + x2[j] for j in range(len(N1))]
    return p_of(N1, N2, z, sum(x1), test, alt, ts)


# ---------------------------------------------------------------------------------------
# Rejection region and power of a fixed design
_REJ = {}


def reject_array(N1, N2, test, alt, ts, alpha):
    """Boolean array over the tables (x1_1, x2_1, x1_2, x2_2, ...) of the design, from the
    p-value of every pair of stratum totals z and statistic S."""
    key = (tuple(N1), tuple(N2), test, alt, ts, alpha)
    if key in _REJ:
        return _REJ[key]
    J = len(N1)
    Smax = sum(N1)
    rej = np.zeros([N1[j] + N2[j] + 1 for j in range(J)] + [Smax + 1], dtype=bool)
    for z in itertools.product(*[range(N1[j] + N2[j] + 1) for j in range(J)]):
        lo = sum(max(0, z[j] - N2[j]) for j in range(J))
        hi = sum(min(z[j], N1[j]) for j in range(J))
        for S in range(lo, hi + 1):
            p = p_of(N1, N2, list(z), S, test, alt, ts)
            rej[z + (S,)] = p < alpha - TOL
    grids = np.meshgrid(*[np.arange(n + 1) for j in range(J) for n in (N1[j], N2[j])],
                        indexing="ij")
    x1 = grids[0::2]
    x2 = grids[1::2]
    idx = tuple(x1[j] + x2[j] for j in range(J)) + (sum(x1),)
    R = rej[idx]
    _REJ[key] = R
    return R


def pvalue_set(N1, N2, test, alt, ts):
    """All p-values of the design, for reported_level()."""
    J = len(N1)
    out = []
    for z in itertools.product(*[range(N1[j] + N2[j] + 1) for j in range(J)]):
        lo = sum(max(0, z[j] - N2[j]) for j in range(J))
        hi = sum(min(z[j], N1[j]) for j in range(J))
        for S in range(lo, hi + 1):
            out.append(p_of(N1, N2, list(z), S, test, alt, ts))
    return out


def table_probs(N1, N2, p1, p2):
    P = np.array(1.0)
    for j in range(len(N1)):
        P = np.multiply.outer(P, binom.pmf(np.arange(N1[j] + 1), N1[j], p1[j]))
        P = np.multiply.outer(P, binom.pmf(np.arange(N2[j] + 1), N2[j], p2[j]))
    return P


def power(N1, N2, p1, p2, test, alt, ts, alpha):
    R = reject_array(N1, N2, test, alt, ts, alpha)
    return float(table_probs(N1, N2, p1, p2)[R].sum())


# ---------------------------------------------------------------------------------------
# Sample size
def raw_n_unstrat(p1, p2, frac, r, alpha, tp, alt, cc):
    """Formula (3.3) of Woolson et al. (1986) for the table of totals, written out as
    printed: group 1 holds the fraction a = sum_j t_j s_j, and the probabilities of the
    groups are the averages over the strata weighted by the group sizes. With cc, formula
    (2.3) of Nam (1992) with one stratum, D = a (1 - a) (pi1 - pi2)."""
    t = np.asarray(frac, float) / np.sum(frac)
    r = np.resize(np.asarray(r, float), len(t))
    s = r / (1 + r)
    a = float(np.sum(t * s))
    b = 1 - a
    pi1 = float(np.sum(t * s * np.asarray(p1, float))) / a
    pi2 = float(np.sum(t * (1 - s) * np.asarray(p2, float))) / b
    diff = pi1 - pi2
    if alt == "less":
        diff = -diff
    if alt == "two.sided":
        diff = abs(diff)
    if not diff > 1e-12:
        return float("nan")
    pipp = pi1 * a + pi2 * b
    za = norm.isf(alpha / 2 if alt == "two.sided" else alpha)
    zb = norm.ppf(tp)
    num = (za * math.sqrt(a * b * max(pipp * (1 - pipp), 0)) +
           zb * math.sqrt(a ** 2 * b ** 2 * (max(pi1 * (1 - pi1), 0) / a +
                                             max(pi2 * (1 - pi2), 0) / b))) ** 2
    n = num / (diff * a * b) ** 2
    if cc:
        D = a * b * diff
        n = (1 + math.sqrt(1 + 2 / (D * n))) ** 2 * n / 4
    return float(n)


def raw_n(p1, p2, frac, r, alpha, tp, alt, cc, pooled=False):
    """Formula (2.8) of Woolson et al. (1986), and (2.3) of Nam (1992) with cc; with
    pooled, the size of the table of totals."""
    if pooled:
        return raw_n_unstrat(p1, p2, frac, r, alpha, tp, alt, cc)
    frac = np.asarray(frac, float) / np.sum(frac)
    r = np.resize(np.asarray(r, float), len(frac))
    s = r / (1 + r)
    t = frac
    p1 = np.asarray(p1, float)
    p2 = np.asarray(p2, float)
    D = np.sum(t * s * (1 - s) * (p1 - p2))
    if alt == "less":
        D = -D
    if alt == "two.sided":
        D = abs(D)
    pbar = s * p1 + (1 - s) * p2
    v0 = np.sum(t * s * (1 - s) * np.maximum(pbar * (1 - pbar), 0))
    v1 = np.sum(t * s * (1 - s) * ((1 - s) * np.maximum(p1 * (1 - p1), 0) +
                                   s * np.maximum(p2 * (1 - p2), 0)))
    za = norm.isf(alpha / 2 if alt == "two.sided" else alpha)
    zb = norm.ppf(tp)
    if not D > 1e-12:
        return float("nan")
    n = (za * math.sqrt(v0) + zb * math.sqrt(v1)) ** 2 / D ** 2
    if cc:
        n = (1 + math.sqrt(1 + 2 / (D * n))) ** 2 * n / 4
    return float(n)


def ceil_tol(x):
    return math.ceil(x - 1e-9)


def sizes(n, frac, r):
    frac = np.asarray(frac, float) / np.sum(frac)
    r = np.resize(np.asarray(r, float), len(frac))
    N2 = [max(1, ceil_tol(f * n / (1 + rr))) for f, rr in zip(frac, r)]
    N1 = [max(1, ceil_tol(rr * m)) for rr, m in zip(r, N2)]
    return N1, N2


def exact_n(p1, p2, frac, r, alpha, tp, test, alt, ts="minlike", search="crossing"):
    cc = test.endswith("-cc")
    n0 = max(1, ceil_tol(raw_n(p1, p2, frac, r, alpha, tp, alt, cc, test in UNSTRAT)))
    pw = lambda n: power(*sizes(n, frac, r), p1, p2, test, alt, ts, alpha)
    if search == "smallest":
        n = 1
        while pw(n) < tp - TOL:
            n += 1
        return n
    n = n0
    P = pw(n)
    if P >= tp - TOL:
        while P >= tp - TOL and n > 1:
            n -= 1
            P = pw(n)
        if P < tp - TOL:
            n += 1
    else:
        while P < tp - TOL:
            n += 1
            P = pw(n)
    return n


# ---------------------------------------------------------------------------------------
# Re-estimation
def split(p, delta, r, effect):
    """Group probabilities with pooled probability p = (r p1 + p2) / (1 + r) and the effect
    delta on the given scale, by a root search for RR and OR."""
    if effect == "RD":
        return p + delta / (1 + r), p - r * delta / (1 + r)
    if p <= 0:
        return 0.0, 0.0
    if effect == "RR":
        p2 = (1 + r) * p / (1 + r * delta)
        return delta * p2, p2
    if p >= 1:
        return 1.0, 1.0
    g = lambda q2: (r * (delta * q2 / (1 - q2 + delta * q2)) + q2) / (1 + r) - p
    q2 = brentq(g, 0.0, 1.0, xtol=1e-15, rtol=8.9e-16)
    return delta * q2 / (1 - q2 + delta * q2), q2


def reestimate(s, n1, n2, DA, effect, frac, r, method, test, alpha, tp, alt, ts,
               N1p, N2p, restricted=False, Nmin=None, Nmax=None, ss_strata="stratified"):
    """Final group sizes after the interim pooled responders s (one per stratum). With
    ss_strata = 'pooled', every stratum uses the proportion of responders among all the
    interim patients."""
    J = len(s)
    r = list(np.resize(np.asarray(r, float), J))
    raw1, raw2 = [], []
    for j in range(J):
        hat = (sum(s) / sum(n1[k] + n2[k] for k in range(J)) if ss_strata == "pooled"
               else s[j] / (n1[j] + n2[j]))
        a, b = split(hat, DA, r[j], effect)
        raw1.append(a)
        raw2.append(b)
    q1 = [min(1, max(0, v)) for v in raw1]
    q2 = [min(1, max(0, v)) for v in raw2]
    cc = test.endswith("-cc")
    un = test in UNSTRAT
    if method == "standard":
        nr = raw_n(raw1, raw2, frac, r, alpha, tp, alt, cc, un)
        ok = math.isfinite(nr)
        n = max(1, ceil_tol(nr)) if ok else None
    else:
        ok = math.isfinite(raw_n(q1, q2, frac, r, alpha, tp, alt, cc, un))
        n = exact_n(q1, q2, frac, r, alpha, tp, test, alt, ts) if ok else None
    low1 = [max(n1[j], N1p[j]) if restricted else n1[j] for j in range(J)]
    low2 = [max(n2[j], N2p[j]) if restricted else n2[j] for j in range(J)]
    if not ok:
        return [max(N1p[j], n1[j]) for j in range(J)], [max(N2p[j], n2[j]) for j in range(J)]
    # smallest nominal size meeting the lower bounds, largest within N.max
    m = 1
    while True:
        a, b = sizes(m, frac, r)
        if all(a[j] >= low1[j] and b[j] >= low2[j] for j in range(J)) and \
                (Nmin is None or sum(a) + sum(b) >= Nmin):
            break
        m += 1
    n = max(n, m)
    if Nmax is not None:
        # the total of sizes(n) is at least n, so the size is at most Nmax
        n = min(n, int(Nmax))
        while sum(sizes(n, frac, r)[0]) + sum(sizes(n, frac, r)[1]) > Nmax:
            n -= 1
    return sizes(n, frac, r)


# ---------------------------------------------------------------------------------------
# Power of a re-estimation design, by enumeration
def bssr_power(n1, n2, fin, p1, p2, test, alt, ts, alpha):
    """fin(s) returns the final group sizes for the interim pooled responders s."""
    J = len(n1)
    tot = 0.0
    EN = 0.0
    ranges = []
    for j in range(J):
        ranges += [range(n1[j] + 1), range(n2[j] + 1)]
    finals = {}
    for idx in itertools.product(*ranges):
        a = idx[0::2]
        b = idx[1::2]
        w = 1.0
        for j in range(J):
            w *= binom.pmf(a[j], n1[j], p1[j]) * binom.pmf(b[j], n2[j], p2[j])
        if w == 0:
            continue
        s = tuple(a[j] + b[j] for j in range(J))
        if s not in finals:
            finals[s] = fin(s)
        N1, N2 = finals[s]
        EN += w * (sum(N1) + sum(N2))
        R = reject_array(N1, N2, test, alt, ts, alpha)
        m1 = [N1[j] - n1[j] for j in range(J)]
        m2 = [N2[j] - n2[j] for j in range(J)]
        P2 = table_probs(m1, m2, p1, p2)
        sl = []
        for j in range(J):
            sl += [slice(a[j], a[j] + m1[j] + 1), slice(b[j], b[j] + m2[j] + 1)]
        tot += w * float((P2 * R[tuple(sl)]).sum())
    return tot, EN


def sup_tie(f, lo=(0.0, 0.0), hi=(1.0, 1.0), G=21, starts=6):
    """Largest value of f over a rectangle: grid, then local optimization from the best
    grid points."""
    xs = np.linspace(lo[0], hi[0], G)
    ys = np.linspace(lo[1], hi[1], G)
    vals = {(x, y): f(x, y) for x in xs for y in ys}
    best = max(vals.items(), key=lambda kv: kv[1])
    for (x0, y0), _ in sorted(vals.items(), key=lambda kv: -kv[1])[:starts]:
        o = minimize(lambda v: -f(v[0], v[1]), [x0, y0], method="L-BFGS-B",
                     bounds=[(lo[0], hi[0]), (lo[1], hi[1])],
                     options={"ftol": 1e-15, "gtol": 1e-12, "maxiter": 500})
        if -o.fun > best[1]:
            best = ((o.x[0], o.x[1]), -o.fun)
    return best


# ---------------------------------------------------------------------------------------
def blocks(emit, reported_level):
    tests = ["CMH", "CMH-cc", "Cochran", "Cochran-cc", "Strat-Fisher"]
    # p-values of individual tables
    f = "test-strat_pvalue_cases.R"
    N1, N2 = [6, 9], [7, 5]
    tabs = [([4, 7], [1, 2]), ([2, 3], [5, 4]), ([6, 9], [0, 5]), ([3, 5], [3, 2]),
            ([2, 9], [1, 2])]
    for test in tests:
        for alt in ("greater", "less", "two.sided"):
            ts_list = ["minlike", "central", "blaker"] if (test == "Strat-Fisher" and
                                                           alt == "two.sided") else ["minlike"]
            for ts in ts_list:
                emit(f, f"{test} {alt} {ts}",
                     [pvalue(test, alt, ts, N1, N2, x1, x2) for x1, x2 in tabs])
    f = "test-StratRR.R"
    # Rejected outcomes over the grid of stratum totals: number of pairs (z, S)
    N1, N2 = [5, 7], [6, 4]
    for test in tests:
        for alt in ("greater", "two.sided"):
            ts = "blaker" if test == "Strat-Fisher" and alt == "two.sided" else "minlike"
            count = 0
            for z in itertools.product(range(N1[0] + N2[0] + 1), range(N1[1] + N2[1] + 1)):
                lo = sum(max(0, z[j] - N2[j]) for j in range(2))
                hi = sum(min(z[j], N1[j]) for j in range(2))
                for S in range(lo, hi + 1):
                    p = p_of(N1, N2, list(z), S, test, alt, ts)
                    count += p < 0.1 - TOL
            emit(f, f"rejected outcomes {test} {alt}", [count], rtol=0)
    # Power of fixed designs with two and three strata
    f = "test-StratPower.R"
    for test in tests:
        for alt in ("greater", "less", "two.sided"):
            ts = "central" if test == "Strat-Fisher" and alt == "two.sided" else "minlike"
            v2 = power([8, 6], [7, 9], [0.6, 0.45], [0.35, 0.3], test, alt, ts, 0.05)
            v3 = power([4, 5, 3], [5, 3, 4], [0.7, 0.3, 0.55], [0.4, 0.25, 0.5], test, alt, ts,
                       0.1)
            emit(f, f"power {test} {alt}", [v2, v3])
    # Normal approximation
    f = "test-strat_raw_n.R"
    emit(f, "Woolson and Nam", [
        raw_n([0.5, 0.7], [0.3, 0.5], [0.4, 0.6], 1, 0.025, 0.8, "greater", False),
        raw_n([0.5, 0.7], [0.3, 0.5], [0.4, 0.6], 1, 0.025, 0.8, "greater", True),
        raw_n([0.3, 0.5, 0.2], [0.45, 0.6, 0.35], [1, 2, 1], 2, 0.05, 0.9, "less", False),
        raw_n([0.6, 0.4], [0.4, 0.3], [0.5, 0.5], [0.5, 2], 0.05, 0.8, "two.sided", True)])
    # Sample size
    f = "test-StratSampleSize.R"
    emit(f, "standard", [ceil_tol(raw_n([0.5, 0.7], [0.3, 0.5], [0.4, 0.6], 1, 0.025, 0.8,
                                         "greater", False)),
                         ceil_tol(raw_n([0.5, 0.7], [0.3, 0.5], [0.4, 0.6], 1, 0.025, 0.8,
                                         "greater", True))], rtol=0)
    emit(f, "exact", [exact_n([0.65, 0.8], [0.3, 0.45], [0.5, 0.5], 1, 0.025, 0.8, t,
                              "greater") for t in ("CMH", "Strat-Fisher")], rtol=0)
    emit(f, "smallest", [exact_n([0.65, 0.8], [0.3, 0.45], [0.5, 0.5], 1, 0.025, 0.8,
                                 "CMH", "greater", search="smallest")], rtol=0)


# ---------------------------------------------------------------------------------------
# Type I error rate profiles for the adjusted level
def outcome_pvalues(n1, n2, fin, theta_pts, test, alt, ts):
    """For each point theta, the p-values of the final tables with their probabilities
    under the null hypothesis p1 = p2 = theta, aggregated by p-value."""
    J = len(n1)
    ranges = []
    for j in range(J):
        ranges += [range(n1[j] + 1), range(n2[j] + 1)]
    pv_cache = {}
    profiles = []
    finals = {}
    for th in theta_pts:
        acc = {}
        for idx in itertools.product(*ranges):
            a = idx[0::2]
            b = idx[1::2]
            w = 1.0
            for j in range(J):
                w *= binom.pmf(a[j], n1[j], th[j]) * binom.pmf(b[j], n2[j], th[j])
            if w == 0:
                continue
            s = tuple(a[j] + b[j] for j in range(J))
            if s not in finals:
                finals[s] = fin(s)
            N1, N2 = finals[s]
            m1 = [N1[j] - n1[j] for j in range(J)]
            m2 = [N2[j] - n2[j] for j in range(J)]
            P2 = table_probs(m1, m2, th, th)
            for sub in itertools.product(*[range(v + 1) for pair in zip(m1, m2) for v in pair]):
                pr = P2[sub]
                if pr == 0:
                    continue
                x1 = [a[j] + sub[2 * j] for j in range(J)]
                x2 = [b[j] + sub[2 * j + 1] for j in range(J)]
                key = (tuple(N1), tuple(N2), tuple(x1[j] + x2[j] for j in range(J)), sum(x1))
                if key not in pv_cache:
                    pv_cache[key] = pvalue(test, alt, ts, list(N1), list(N2), x1, x2)
                p = pv_cache[key]
                acc[p] = acc.get(p, 0.0) + w * pr
        ps = np.array(sorted(acc))
        profiles.append((ps, np.cumsum([acc[p] for p in ps])))
    return profiles, finals


def max_rate(profiles, level):
    best = 0.0
    for ps, cum in profiles:
        k = np.searchsorted(ps, level - TOL, side="left")  # p < level - TOL
        best = max(best, float(cum[k - 1]) if k > 0 else 0.0)
    return best


def bisect(profiles, alpha, tol=1e-8):
    if max_rate(profiles, alpha) <= alpha:
        return alpha
    lo, hi = 0.0, alpha
    while hi - lo > tol * alpha:
        mid = (lo + hi) / 2
        if max_rate(profiles, mid) <= alpha:
            lo = mid
        else:
            hi = mid
    return lo


def blocks_bssr(emit, reported_level):
    # Re-estimation designs with two strata: power and expected size
    f = "test-StratPowerBSSR.R"
    designs = {
        "CMH standard": dict(n1=[5, 5], n2=[5, 5], N1p=[10, 10], N2p=[10, 10],
                             frac=[0.5, 0.5], r=[1, 1], DA=0.3, effect="RD", test="CMH",
                             alt="greater", ts="minlike", method="standard", alpha=0.05),
        "Fisher exact": dict(n1=[4, 5], n2=[4, 5], N1p=[8, 10], N2p=[8, 10],
                             frac=[8, 10], r=[1, 1], DA=0.4, effect="RD",
                             test="Strat-Fisher", alt="greater", ts="minlike",
                             method="exact", alpha=0.1, Nmax=48),
    }
    fins = {}
    for key, d in designs.items():
        fin = (lambda s, d=d: reestimate(s, d["n1"], d["n2"], d["DA"], d["effect"],
                                         d["frac"], d["r"], d["method"], d["test"],
                                         d["alpha"], 0.8, d["alt"], d["ts"], d["N1p"],
                                         d["N2p"], Nmax=d.get("Nmax")))
        fins[key] = fin
        res = []
        for p1, p2 in (([0.55, 0.75], [0.25, 0.45]), ([0.4, 0.6], [0.4, 0.6])):
            pw, EN = bssr_power(d["n1"], d["n2"], fin, p1, p2, d["test"], d["alt"],
                                d["ts"], d["alpha"])
            res += [pw, EN]
        emit(f, "power and E.N " + key, res)
        # Final sizes of a few interim outcomes
        # Interim outcomes that are possible in both designs (stratum 1 of the second
        # design has 8 patients at the interim analysis)
        cells = [(0, 0), (3, 7), (6, 4), (8, 10)]
        emit("test-strat_reestimate.R", "final sizes " + key,
             [v for s in cells for part in fin(s) for v in part], rtol=0)
    # Final sizes with the assumed effect on the scale of the risk ratio and of the odds
    # ratio, where the pooled probability is divided by a root search
    for key, DA, effect in (("CMH standard RR", 1.5, "RR"), ("CMH standard OR", 3.0, "OR")):
        d = dict(designs["CMH standard"], DA=DA, effect=effect)
        emit("test-strat_reestimate.R", "final sizes " + key,
             [v for s in cells
              for part in reestimate(s, d["n1"], d["n2"], d["DA"], d["effect"], d["frac"],
                                     d["r"], d["method"], d["test"], d["alpha"], 0.8,
                                     d["alt"], d["ts"], d["N1p"], d["N2p"])
              for v in part], rtol=0)
    # Type I error rate on a grid and its largest value
    f = "test-StratTypeIErrorBSSR.R"
    d = designs["CMH standard"]
    fin = fins["CMH standard"]
    g = lambda x, y: bssr_power(d["n1"], d["n2"], fin, [x, y], [x, y], d["test"], d["alt"],
                                d["ts"], d["alpha"])[0]
    emit(f, "grid CMH standard", [g(x, y) for x, y in ((0.1, 0.3), (0.5, 0.5), (0.9, 0.2))])
    emit(f, "fixed CMH standard", [power([10, 10], [10, 10], [x, y], [x, y], "CMH",
                                         "greater", "minlike", 0.05)
                                   for x, y in ((0.1, 0.3), (0.5, 0.5), (0.9, 0.2))])
    (bx, by), bv = sup_tie(g, (0.0, 0.0), (1.0, 1.0), G=21, starts=8)
    t = lambda x, y: power([10, 10], [10, 10], [x, y], [x, y], "CMH", "greater", "minlike",
                           0.05)
    (tx, ty), tv = sup_tie(t, (0.0, 0.0), (1.0, 1.0), G=21, starts=8)
    emit(f, "certified max CMH standard", [bv, tv], rtol=1e-8)
    # Re-estimation of an ongoing trial
    f = "test-StratBSSR.R"
    N1, N2 = reestimate((9, 16), [10, 12], [10, 12], 0.2, "RD", [20, 24], [1, 1],
                        "standard", "CMH", 0.025, 0.8, "greater", "minlike", [20, 24],
                        [20, 24])
    emit(f, "final sizes CMH", list(N1) + list(N2), rtol=0)
    q1 = [9 / 20 + 0.1, 16 / 24 + 0.1]
    q2 = [9 / 20 - 0.1, 16 / 24 - 0.1]
    emit(f, "power CMH", [power(list(N1), list(N2), q1, q2, "CMH", "greater", "minlike",
                                0.025)])
    # Adjusted level, maximum over the grid
    f = "test-StratAlphaAdjBSSR.R"
    d = dict(n1=[3, 3], n2=[3, 3], N1p=[6, 6], N2p=[6, 6], frac=[0.5, 0.5], r=[1, 1],
             DA=0.35, effect="RD", test="CMH", alt="greater", ts="minlike",
             method="standard", alpha=0.1)
    fin = (lambda s: reestimate(s, d["n1"], d["n2"], d["DA"], d["effect"], d["frac"],
                                d["r"], d["method"], d["test"], d["alpha"], 0.8, d["alt"],
                                d["ts"], d["N1p"], d["N2p"]))
    th = np.round(np.arange(0.1, 0.95, 0.2), 10)
    pts = [(x, y) for y in th for x in th]
    prof, finals = outcome_pvalues(d["n1"], d["n2"], fin, pts, "CMH", "greater", "minlike")
    m0 = max_rate(prof, 0.1)
    lev = bisect(prof, 0.1)
    pv = []
    for N1, N2 in set((tuple(a), tuple(b)) for a, b in finals.values()):
        pv.append(pvalue_set(list(N1), list(N2), "CMH", "greater", "minlike"))
    rep = reported_level(pv, lev) if lev < 0.1 else lev
    fixed = [(np.array([0]), np.array([0.0]))]
    tprof = []
    for x, y in pts:
        acc = {}
        P = table_probs([6, 6], [6, 6], [x, y], [x, y])
        for idx in itertools.product(*[range(7)] * 4):
            p = pvalue("CMH", "greater", "minlike", [6, 6], [6, 6], [idx[0], idx[2]],
                       [idx[1], idx[3]])
            acc[p] = acc.get(p, 0.0) + P[idx]
        ps = np.array(sorted(acc))
        tprof.append((ps, np.cumsum([acc[p] for p in ps])))
    t0 = max_rate(tprof, 0.1)
    tlev = bisect(tprof, 0.1)
    trep = reported_level([pvalue_set([6, 6], [6, 6], "CMH", "greater", "minlike")], tlev) \
        if tlev < 0.1 else tlev
    del fixed
    emit(f, "grid max at nominal level", [m0, t0], rtol=1e-10)
    emit(f, "grid adjusted level", [rep, trep], rtol=1e-12)
    emit(f, "grid max at adjusted level", [max_rate(prof, rep), max_rate(tprof, trep)],
         rtol=1e-10)
    # Three strata: exact power of a small re-estimation design, for the simulation test
    f = "test-strat_simulate.R"
    d = dict(n1=[2, 2, 2], n2=[2, 2, 2], N1p=[4, 4, 4], N2p=[4, 4, 4], frac=[1, 1, 1],
             r=[1, 1, 1], DA=0.4, effect="RD", test="CMH", alt="greater", ts="minlike",
             method="standard", alpha=0.1)
    fin = (lambda s: reestimate(s, d["n1"], d["n2"], d["DA"], d["effect"], d["frac"],
                                d["r"], d["method"], d["test"], d["alpha"], 0.8, d["alt"],
                                d["ts"], d["N1p"], d["N2p"], Nmax=36))
    pw, EN = bssr_power(d["n1"], d["n2"], fin, [0.6, 0.7, 0.5], [0.3, 0.4, 0.2], "CMH",
                        "greater", "minlike", 0.1)
    emit(f, "exact power and E.N three strata", [pw, EN])


# ---------------------------------------------------------------------------------------
# Unstratified tests and the re-estimation from the pooled estimate (added for 3.0.0, A2)
def blocks_unstrat(emit, reported_level):
    alts = ("greater", "less", "two.sided")
    # p-values of individual tables
    f = "test-strat_pvalue_cases.R"
    N1, N2 = [6, 9], [7, 5]
    tabs = [([4, 7], [1, 2]), ([2, 3], [5, 4]), ([6, 9], [0, 5]), ([3, 5], [3, 2]),
            ([2, 9], [1, 2])]
    for test in UNSTRAT:
        for alt in alts:
            ts_list = ["minlike", "central", "blaker"] if (test == "Unstrat-Fisher" and
                                                           alt == "two.sided") else ["minlike"]
            for ts in ts_list:
                emit(f, f"{test} {alt} {ts}",
                     [pvalue(test, alt, ts, N1, N2, x1, x2) for x1, x2 in tabs])
    # Rejected outcomes over the grid of stratum totals
    f = "test-StratRR.R"
    N1, N2 = [5, 7], [6, 4]
    for test in UNSTRAT:
        for alt in ("greater", "two.sided"):
            ts = "blaker" if test == "Unstrat-Fisher" and alt == "two.sided" else "minlike"
            count = 0
            for z in itertools.product(range(N1[0] + N2[0] + 1), range(N1[1] + N2[1] + 1)):
                lo = sum(max(0, z[j] - N2[j]) for j in range(2))
                hi = sum(min(z[j], N1[j]) for j in range(2))
                for S in range(lo, hi + 1):
                    count += p_of(N1, N2, list(z), S, test, alt, ts) < 0.1 - TOL
            emit(f, f"rejected outcomes {test} {alt}", [count], rtol=0)
    # Power of fixed designs with two and three strata
    f = "test-StratPower.R"
    for test in UNSTRAT:
        for alt in alts:
            ts = "central" if test == "Unstrat-Fisher" and alt == "two.sided" else "minlike"
            v2 = power([8, 6], [7, 9], [0.6, 0.45], [0.35, 0.3], test, alt, ts, 0.05)
            v3 = power([4, 5, 3], [5, 3, 4], [0.7, 0.3, 0.55], [0.4, 0.25, 0.5], test, alt, ts,
                       0.1)
            emit(f, f"power {test} {alt}", [v2, v3])
    # Normal approximation of the table of totals
    f = "test-strat_raw_n.R"
    emit(f, "pooled", [
        raw_n([0.5, 0.7], [0.3, 0.5], [0.4, 0.6], 1, 0.025, 0.8, "greater", False, True),
        raw_n([0.5, 0.7], [0.3, 0.5], [0.4, 0.6], 1, 0.025, 0.8, "greater", True, True),
        raw_n([0.3, 0.5, 0.2], [0.45, 0.6, 0.35], [1, 2, 1], 2, 0.05, 0.9, "less", False,
              True),
        raw_n([0.6, 0.4], [0.4, 0.3], [0.5, 0.5], [0.5, 2], 0.05, 0.8, "two.sided", True,
              True)])
    # Sample size
    f = "test-StratSampleSize.R"
    emit(f, "standard unstratified",
         [ceil_tol(raw_n([0.5, 0.7], [0.3, 0.5], [0.4, 0.6], [1, 2], 0.025, 0.8, "greater",
                         cc, True)) for cc in (False, True)], rtol=0)
    emit(f, "exact unstratified", [exact_n([0.65, 0.8], [0.3, 0.45], [0.5, 0.5], 1, 0.025,
                                           0.8, t, "greater")
                                   for t in ("Unstrat-Chisq", "Unstrat-Fisher")], rtol=0)
    # Re-estimation from the pooled estimate, and with the unstratified test
    d = dict(n1=[5, 5], n2=[5, 5], N1p=[10, 10], N2p=[10, 10], frac=[0.5, 0.5], r=[1, 1],
             DA=0.3, effect="RD", alt="greater", ts="minlike", method="standard",
             alpha=0.05)
    cells = [(0, 0), (3, 7), (6, 4), (8, 10)]
    for key, test, strata in (("CMH standard pooled", "CMH", "pooled"),
                              ("Unstrat-Chisq standard", "Unstrat-Chisq", "stratified"),
                              ("Unstrat-Chisq standard pooled", "Unstrat-Chisq", "pooled")):
        emit("test-strat_reestimate.R", "final sizes " + key,
             [v for s in cells
              for part in reestimate(s, d["n1"], d["n2"], d["DA"], d["effect"], d["frac"],
                                     d["r"], d["method"], test, d["alpha"], 0.8, d["alt"],
                                     d["ts"], d["N1p"], d["N2p"], ss_strata=strata)
              for v in part], rtol=0)
    # Power of re-estimation designs
    f = "test-StratPowerBSSR.R"
    designs = {
        "CMH standard pooled": dict(d, test="CMH", strata="pooled"),
        "Unstrat-Fisher exact pooled": dict(n1=[4, 5], n2=[4, 5], N1p=[8, 10],
                                            N2p=[8, 10], frac=[8, 10], r=[1, 1], DA=0.4,
                                            effect="RD", test="Unstrat-Fisher",
                                            alt="greater", ts="minlike", method="exact",
                                            alpha=0.1, Nmax=48, strata="pooled"),
    }
    for key, e in designs.items():
        fin = (lambda s, e=e: reestimate(s, e["n1"], e["n2"], e["DA"], e["effect"],
                                         e["frac"], e["r"], e["method"], e["test"],
                                         e["alpha"], 0.8, e["alt"], e["ts"], e["N1p"],
                                         e["N2p"], Nmax=e.get("Nmax"),
                                         ss_strata=e["strata"]))
        res = []
        for p1, p2 in (([0.55, 0.75], [0.25, 0.45]), ([0.4, 0.6], [0.4, 0.6])):
            pw, EN = bssr_power(e["n1"], e["n2"], fin, p1, p2, e["test"], e["alt"],
                                e["ts"], e["alpha"])
            res += [pw, EN]
        emit(f, "power and E.N " + key, res)
    # Type I error rate of the unstratified test when the allocation ratio differs between
    # the strata (r = 1 and 3): the table of totals mixes strata with different response
    # probabilities in different proportions
    f = "test-StratTypeIErrorBSSR.R"
    e = dict(n1=[3, 6], n2=[3, 2], N1p=[6, 9], N2p=[6, 3], frac=[12, 12], r=[1, 3],
             DA=0.4, effect="RD", test="Unstrat-Chisq", alt="greater", ts="minlike",
             method="standard", alpha=0.05, Nmax=72)
    fin = (lambda s: reestimate(s, e["n1"], e["n2"], e["DA"], e["effect"], e["frac"],
                                e["r"], e["method"], e["test"], e["alpha"], 0.8, e["alt"],
                                e["ts"], e["N1p"], e["N2p"], Nmax=e["Nmax"]))
    pts = ((0.1, 0.3), (0.5, 0.5), (0.2, 0.8))
    emit(f, "grid Unstrat-Chisq unequal allocation",
         [bssr_power(e["n1"], e["n2"], fin, [x, y], [x, y], e["test"], e["alt"], e["ts"],
                     e["alpha"])[0] for x, y in pts])
    emit(f, "fixed Unstrat-Chisq unequal allocation",
         [power(e["N1p"], e["N2p"], [x, y], [x, y], e["test"], e["alt"], e["ts"],
                e["alpha"]) for x, y in pts])


# ---------------------------------------------------------------------------------------
# Conditional rejection probabilities of a two-stratum re-estimation design (added for
# 3.0.0, A1)
def hyp_pmf(m1, m2, z):
    """Distribution of the group 1 responders among z responders of m1 + m2 patients."""
    lo = max(0, z - m2)
    if m1 + m2 == 0:
        return 0, np.array([1.0])
    x = np.arange(lo, min(z, m1) + 1)
    return lo, hypergeom.pmf(x, m1 + m2, m1, z)


def cond_reject(n1, n2, fin, test, alt, ts, alpha, pts=()):
    """For every interim cell s = (s1, s2): the final sizes, and for every second-stage
    outcome w = (w1, w2) the conditional rejection probability CRP(s, w) and the
    conditional size CRP.total(s + w) of the final design, by direct summation of the
    hypergeometric distributions; with points (theta1, theta2), the type I error rate given
    s. Returns a list of cells in the order s1 fastest, each a dict."""
    out = []
    for s2 in range(n1[1] + n2[1] + 1):
        for s1 in range(n1[0] + n2[0] + 1):
            N1, N2 = fin((s1, s2))
            m1 = [N1[j] - n1[j] for j in range(2)]
            m2 = [N2[j] - n2[j] for j in range(2)]
            rows = []
            for w2 in range(m1[1] + m2[1] + 1):
                for w1 in range(m1[0] + m2[0] + 1):
                    z = [s1 + w1, s2 + w2]
                    # indicator of rejection for every S given z
                    Smax = sum(N1)
                    R = np.array([p_of(list(N1), list(N2), z, S, test, alt, ts) < alpha - TOL
                                  for S in range(Smax + 1)], dtype=float)
                    # interim and second-stage group 1 counts of each stratum
                    dists = []
                    for j, (s, w) in enumerate(((s1, w1), (s2, w2))):
                        la, fa = hyp_pmf(n1[j], n2[j], s)
                        lu, fu = hyp_pmf(m1[j], m2[j], w)
                        dists.append((la + lu, np.convolve(fa, fu)))
                    (l0, f0), (l1, f1) = dists
                    crp = 0.0
                    for i0, p0 in enumerate(f0):
                        for i1, p1 in enumerate(f1):
                            crp += p0 * p1 * R[l0 + i0 + l1 + i1]
                    (t0, g0), (t1, g1) = hyp_pmf(N1[0], N2[0], z[0]), hyp_pmf(N1[1], N2[1], z[1])
                    tot = 0.0
                    for i0, p0 in enumerate(g0):
                        for i1, p1 in enumerate(g1):
                            tot += p0 * p1 * R[t0 + i0 + t1 + i1]
                    rows.append((w1, w2, crp, tot))
            tie_s = []
            for th1, th2 in pts:
                acc = 0.0
                for w1, w2, crp, tot in rows:
                    acc += (binom.pmf(w1, m1[0] + m2[0], th1) *
                            binom.pmf(w2, m1[1] + m2[1], th2) * crp)
                tie_s.append(acc)
            out.append(dict(s=(s1, s2), N1=N1, N2=N2, rows=rows, tie_s=tie_s))
    return out


def blocks_crp(emit, reported_level):
    f = "test-StratCondRejectBSSR.R"
    d = dict(n1=[3, 3], n2=[3, 3], N1p=[6, 6], N2p=[6, 6], frac=[0.5, 0.5], r=[1, 1],
             DA=0.35, effect="RD", alt="greater", ts="minlike", method="standard",
             alpha=0.05)
    pts = [(0.3, 0.5), (0.5, 0.5)]
    for key, test in (("CMH", "CMH"), ("Strat-Fisher", "Strat-Fisher"),
                      ("Unstrat-Chisq", "Unstrat-Chisq")):
        fin = (lambda s, test=test: reestimate(s, d["n1"], d["n2"], d["DA"], d["effect"],
                                               d["frac"], d["r"], d["method"], test,
                                               d["alpha"], 0.8, d["alt"], d["ts"], d["N1p"],
                                               d["N2p"], Nmax=36))
        cells = cond_reject(d["n1"], d["n2"], fin, test, d["alt"], d["ts"], d["alpha"], pts)
        cmax = [max(r[2] for r in c["rows"]) for c in cells]
        cmin = [min(r[2] for r in c["rows"]) for c in cells]
        tmax = [max(r[3] for r in c["rows"]) for c in cells]
        emit(f, f"sums of CRP.max CRP.min CRP.total.max {key}",
             [sum(cmax), sum(cmin), sum(tmax)])
        # largest CRP over all outcomes, first in the order of the cells and of (w1, w2)
        best = max(cmax)
        loc = None
        for c in cells:
            for w1, w2, crp, tot in c["rows"]:
                if crp >= best * (1 - 1e-12):
                    loc = list(c["s"]) + [w1, w2]
                    break
            if loc is not None:
                break
        emit(f, f"largest CRP and CRP.total {key}", [best, max(tmax)])
        emit(f, f"location of the largest CRP {key}", loc, rtol=0)
        exceed = sum(r[2] - d["alpha"] > TOL for c in cells for r in c["rows"])
        exceed_t = sum(r[3] - d["alpha"] > TOL for c in cells for r in c["rows"])
        emit(f, f"exceed {key}", [exceed, exceed_t], rtol=0)
        n = [d["n1"][j] + d["n2"][j] for j in range(2)]
        tie = []
        for p, (th1, th2) in enumerate(pts):
            tie.append(sum(binom.pmf(c["s"][0], n[0], th1) * binom.pmf(c["s"][1], n[1], th2) *
                           c["tie_s"][p] for c in cells))
        emit(f, f"TIE at the points {key}", tie)
        cell = next(c for c in cells if c["s"] == (2, 3))
        emit(f, f"TIE.s of the cell (2, 3) {key}", cell["tie_s"])
    # Number of rows of the full listing of a design with large re-estimated totals,
    # which StratCondRejectBSSR refuses (more than 1e6)
    rows = 0
    for s2 in range(d["n1"][1] + d["n2"][1] + 1):
        for s1 in range(d["n1"][0] + d["n2"][0] + 1):
            a, b = reestimate((s1, s2), d["n1"], d["n2"], 0.1, d["effect"], d["frac"],
                              d["r"], d["method"], "CMH", d["alpha"], 0.8, d["alt"],
                              d["ts"], d["N1p"], d["N2p"], Nmax=1000)
            rows += ((a[0] + b[0] - d["n1"][0] - d["n2"][0] + 1) *
                     (a[1] + b[1] - d["n1"][1] - d["n2"][1] + 1))
    emit(f, "rows of the full listing with Delta.A = 0.1 and N.max = 1000", [rows], rtol=0)


# ---------------------------------------------------------------------------------------
# Random stratum and group sizes (two strata)
def alloc_group_sizes(model, n, b, B):
    """Probabilities of the size m of group 1 among the n patients of a stratum: binomial
    with the probability b, or from permuted blocks of size B holding round(B b) patients
    of group 1, of which the last block is incomplete when n is not a multiple of B."""
    if model == "stratum":
        g = round(B * b)
        k, e = divmod(n, B)
        return {k * g + h: float(hypergeom.pmf(h, B, g, e))
                for h in range(max(0, e - (B - g)), min(e, g) + 1)}
    return {m: float(binom.pmf(m, n, b)) for m in range(n + 1)}


def alloc_power(N, a1, b, model, p1, p2, test, alt, ts, alpha, B=4):
    """Exact power with random sizes: the sum over every stratum size and group size of
    the probability times the power of the fixed design. Model 'group' fixes
    n_1 = floor(N a1); 'stratum' and 'both' draw n_1 from Bin(N, a1). A stratum without
    patients is left out of the design."""
    if model == "group":
        strata = {math.floor(N * a1 + 1e-9): 1.0}
    else:
        strata = {n1: float(binom.pmf(n1, N, a1)) for n1 in range(N + 1)}
    total = 0.0
    for n1, wn in strata.items():
        n = (n1, N - n1)
        dists = [alloc_group_sizes(model, n[j], b[j], B) for j in range(2)]
        keep = [j for j in range(2) if n[j] > 0]
        for m0, w0 in dists[0].items():
            for m1, w1 in dists[1].items():
                m = (m0, m1)
                N1 = [m[j] for j in keep]
                N2 = [n[j] - m[j] for j in keep]
                total += wn * w0 * w1 * power(N1, N2, [p1[j] for j in keep],
                                              [p2[j] for j in keep], test, alt, ts, alpha)
    return total


def alloc_exact_n(p1, p2, frac, r, alpha, tp, test, alt, model, B=4):
    """Total sample size of the search 'crossing' with the power of alloc_power: from the
    rounded-up normal approximation, down while the power attains the target, or up
    until it does."""
    t = np.asarray(frac, float) / np.sum(frac)
    rr = np.resize(np.asarray(r, float), 2)
    b = rr / (1 + rr)
    cc = test.endswith("-cc")
    pw = lambda N: alloc_power(N, t[0], b, model, p1, p2, test, alt, "minlike", alpha, B)
    n = max(1, ceil_tol(raw_n(p1, p2, frac, r, alpha, tp, alt, cc, test in UNSTRAT)))
    P = pw(n)
    if P >= tp - TOL:
        while P >= tp - TOL and n > 1:
            n -= 1
            P = pw(n)
        if P < tp - TOL:
            n += 1
    else:
        while P < tp - TOL:
            n += 1
            P = pw(n)
    return n


def blocks_alloc(emit, reported_level):
    f = "test-StratPower.R"
    sc1 = ([0.5, 0.8], [0.2, 0.4])
    sc0 = ([0.3, 0.6], [0.3, 0.6])
    for model in ("stratum", "group", "both"):
        for test in ("CMH", "Strat-Fisher", "Unstrat-Chisq"):
            emit(f, f"random allocation {model} {test}",
                 [alloc_power(16, 0.5, (0.5, 0.5), model, p1, p2, test, "greater",
                              "minlike", 0.05) for p1, p2 in (sc1, sc0)])
    # Unequal fractions and allocation ratio 3: n_1 = floor(6.8) = 6 with 'group', and
    # blocks of 4 with 3 patients of group 1
    for model in ("stratum", "group"):
        emit(f, f"random allocation {model} Strat-Fisher, frac (0.4, 0.6), r 3",
             alloc_power(17, 0.4, (0.75, 0.75), model, [0.5, 0.8], [0.2, 0.4],
                         "Strat-Fisher", "greater", "minlike", 0.05))
    # One key for both sizes, so that a shifted size cannot match the other one
    emit("test-StratSampleSize.R", "random allocation stratum and group CMH, crossing",
         [alloc_exact_n([0.6, 0.8], [0.2, 0.4], [0.5, 0.5], 1, 0.05, 0.8, "CMH", "greater",
                        model) for model in ("stratum", "group")], rtol=0)


# ---------------------------------------------------------------------------------------
# Re-estimation designs given an allocation sequence (two strata)
# Two sequences of 26 patients, (stratum, group) in the order of enrolment; the first
# has no patient of group 2 in stratum 2 among the first 12
SEQ_A = [(1, 1), (2, 1), (1, 2), (2, 1), (1, 1), (1, 2), (2, 1), (1, 2), (2, 1), (1, 1),
         (1, 2), (2, 1), (2, 2), (1, 1), (2, 2), (1, 2), (2, 2), (1, 1), (2, 1), (2, 2),
         (1, 2), (1, 1), (2, 2), (2, 1), (1, 1), (2, 2)]
SEQ_B = [(2, 2), (1, 1), (1, 2), (2, 1), (1, 1), (2, 1), (1, 2), (1, 2), (2, 2), (1, 1),
         (2, 1), (1, 2), (2, 1), (1, 1), (2, 2), (2, 2), (1, 2), (2, 1), (1, 1), (2, 2),
         (1, 1), (2, 1), (1, 2), (2, 2), (2, 1), (1, 2)]


def seq_counts(seq):
    """Sizes of the groups (N1 and N2 of stratum 1, N1 and N2 of stratum 2) after every
    number of patients of a sequence of (stratum, group)."""
    counts = [[0, 0, 0, 0]]
    for st, g in seq:
        c = list(counts[-1])
        c[2 * (st - 1) + g - 1] += 1
        counts.append(c)
    return counts


def seq_nfin(s):
    """Final total of the interim cell s, the same rule for both sequences."""
    return min(26, 14 + 2 * (s[0] % 3) + s[1])


def seq_power(seq, n, p1, p2, test, alt, ts, alpha):
    """Power of the re-estimation design whose interim analysis has the first n patients of
    the sequence and whose final analysis has the first seq_nfin(s) patients."""
    counts = seq_counts(seq)
    c = counts[n]
    fin = lambda s: ([counts[seq_nfin(s)][0], counts[seq_nfin(s)][2]],
                     [counts[seq_nfin(s)][1], counts[seq_nfin(s)][3]])
    return bssr_power([c[0], c[2]], [c[1], c[3]], fin, p1, p2, test, alt, ts, alpha)[0]


def blocks_alloc_seq(emit, reported_level):
    f = "test-strat_seq_power.R"
    for test in ("CMH", "Strat-Fisher", "Unstrat-Chisq"):
        emit(f, f"power given the sequences {test}",
             [seq_power(sq, 12, p1, p2, test, "greater", "minlike", 0.05)
              for sq in (SEQ_A, SEQ_B)
              for p1, p2 in (([0.5, 0.8], [0.2, 0.4]), ([0.3, 0.6], [0.3, 0.6]))])
    f = "test-strat_seq_tie_coef.R"
    for test in ("CMH", "Strat-Fisher"):
        emit(f, f"type I error rate given the sequences {test}",
             [seq_power(sq, 12, th, th, test, "greater", "minlike", 0.05)
              for sq in (SEQ_A, SEQ_B) for th in ([0.3, 0.6], [0.7, 0.2])])
