"""Independent reference values for the paired (matched-pairs) functions of bbssr.

Called from reference_values.py, which passes its emit() and reported_level(). The
computations follow the definitions directly and share no code with the package: the
restricted maximum likelihood estimate of the score test is found by a root search of the
score equation (not by the closed form of Tango), the exact unconditional p-values by a
fine grid with local optimization, the power of a re-estimation design by enumerating the
interim and second-stage multinomial outcomes, and the largest type I error rate by a fine
grid with local optimization.
"""
import math

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import gammaln
from scipy.stats import beta as beta_dist
from scipy.stats import binom, norm

TOL = 1.4901161193847656e-08  # tolerance of fpCompare


# ---------------------------------------------------------------------------------------
# Outcomes, statistics and p-values
def outcome_arrays(N):
    a, b = np.meshgrid(np.arange(N + 1), np.arange(N + 1), indexing="ij")
    ok = a + b <= N
    return a, b, ok


def restricted_p01(a, b, N, d0):
    """Maximum likelihood estimate of p01 under p10 - p01 = d0, by a root search of the
    score equation over the feasible interval, or at an end point."""
    c = N - a - b
    lo, hi = max(0.0, -d0), (1.0 - d0) / 2.0

    def loglik(p):
        val = 0.0
        for k, q in ((a, p + d0), (b, p), (c, 1.0 - 2.0 * p - d0)):
            if k > 0:
                if q <= 0:
                    return -math.inf
                val += k * math.log(q)
        return val

    def score(p):
        s = 0.0
        if a > 0:
            s += a / (p + d0)
        if b > 0:
            s += b / p
        if c > 0:
            s -= 2.0 * c / (1.0 - 2.0 * p - d0)
        return s

    eps = 1e-15
    x0, x1 = lo + eps * (1 + abs(lo)), hi - eps * (1 + abs(hi))
    f0, f1 = score(x0), score(x1)
    if f0 <= 0:
        cand = lo
    elif f1 >= 0:
        cand = hi
    else:
        cand = brentq(score, x0, x1, xtol=1e-15, rtol=8.9e-16, maxiter=500)
    # Guard: compare with the end points
    best = max([lo, hi, cand], key=loglik)
    return best


def score_z(N, margin):
    """Score statistic of the null hypothesis p10 - p01 = -margin at every outcome."""
    Z = np.full((N + 1, N + 1), np.nan)
    d0 = -margin
    for a in range(N + 1):
        for b in range(N + 1 - a):
            if margin == 0:
                Z[a, b] = 0.0 if a + b == 0 else (a - b) / math.sqrt(a + b)
                continue
            p01 = restricted_p01(a, b, N, d0)
            p10 = p01 + d0
            Z[a, b] = (a - b - N * d0) / math.sqrt(N * (p10 + p01 - d0 ** 2))
    return Z


def cond_pvalues(N, alt, midp):
    P = np.full((N + 1, N + 1), np.nan)
    for a in range(N + 1):
        for b in range(N + 1 - a):
            D = a + b
            f = binom.pmf(a, D, 0.5)
            up = binom.sf(a - 1, D, 0.5)
            lo = binom.cdf(a, D, 0.5)
            if alt == "greater":
                P[a, b] = up - f / 2 if midp else up
            elif midp:
                P[a, b] = 1 - f / 2 if a == b else 2 * min(lo, up) - f
            else:
                P[a, b] = min(1.0, 2 * min(lo, up))
    return P


def tail_sup(h, N, lo=0.0, hi=1.0, G=4001):
    """Supremum over psi in [lo, hi] of sum_D h[D] b(D; N, psi), by a fine grid and a local
    optimization around the largest grid values."""
    grid = np.linspace(lo, hi, G)
    B = binom.pmf(np.arange(N + 1)[:, None], N, grid[None, :])
    vals = h @ B
    f = lambda t: float(h @ binom.pmf(np.arange(N + 1), N, t))
    best = vals.max()
    order = np.argsort(vals)[::-1][:5]
    for i in order:
        l, r = grid[max(0, i - 1)], grid[min(G - 1, i + 1)]
        if r > l:
            o = minimize_scalar(lambda t: -f(t), bounds=(l, r), method="bounded",
                                options={"xatol": 1e-13})
            best = max(best, -o.fun)
    return best


def uncond_pvalues(N, alt, gamma=0.0):
    """Exact unconditional p-values ordered by the McNemar statistic (its absolute value for
    the two-sided test), maximized over the probability of a discordant pair."""
    Z = score_z(N, 0)
    stat = Z if alt == "greater" else np.abs(Z)
    a, b, ok = outcome_arrays(N)
    cells = [(int(x), int(y)) for x, y in zip(a[ok], b[ok])]
    P = np.full((N + 1, N + 1), np.nan)
    for (x, y) in cells:
        s = stat[x, y]
        h = np.zeros(N + 1)
        for (u, v) in cells:
            t = stat[u, v]
            if t >= s - 1e-10 * max(abs(s), abs(t)):
                h[u + v] += binom.pmf(u, u + v, 0.5)
        if gamma > 0:
            D = x + y
            lo = 0.0 if D == 0 else beta_dist.ppf(gamma / 2, D, N - D + 1)
            hi = 1.0 if D == N else beta_dist.ppf(1 - gamma / 2, D + 1, N - D)
            P[x, y] = min(1.0, tail_sup(h, N, lo, hi, G=801) + gamma)
        else:
            P[x, y] = min(1.0, tail_sup(h, N))
    return P


def pvalues(N, test, alt, margin=0.0, gamma=0.0):
    if alt == "equivalence":
        return np.fmax(pvalues(N, test, "greater", margin[0]),
                       pvalues(N, test, "less", margin[1]))
    if alt == "less":
        return pvalues(N, test, "greater", margin).T
    if test in ("McNemar", "Score"):
        Z = score_z(N, 0 if test == "McNemar" else margin)
        return norm.sf(Z) if alt == "greater" else np.minimum(1.0, 2 * norm.sf(np.abs(Z)))
    if test == "exact":
        return cond_pvalues(N, alt, False)
    if test == "midP":
        return cond_pvalues(N, alt, True)
    return uncond_pvalues(N, alt, gamma)


_PV = {}


def pvalues_cached(N, test, alt, margin=0.0, gamma=0.0):
    key = (N, test, alt, tuple(np.atleast_1d(margin)), gamma)
    if key not in _PV:
        _PV[key] = pvalues(N, test, alt, margin, gamma)
    return _PV[key]


def reject(N, test, alt, alpha, margin=0.0, gamma=0.0):
    P = pvalues_cached(N, test, alt, margin, gamma)
    with np.errstate(invalid="ignore"):
        R = P < alpha - TOL
    return R & ~np.isnan(P)


def trinomial(N, p10, p01):
    a, b, ok = outcome_arrays(N)
    c = N - a - b
    pc = max(0.0, 1.0 - p10 - p01)
    with np.errstate(divide="ignore", invalid="ignore"):
        lp = (gammaln(N + 1) - gammaln(a + 1) - gammaln(b + 1) - gammaln(np.where(ok, c, 0) + 1)
              + np.where(a > 0, a * np.log(p10 if p10 > 0 else 1e-300), 0.0)
              + np.where(b > 0, b * np.log(p01 if p01 > 0 else 1e-300), 0.0)
              + np.where(c > 0, np.where(ok, c, 0) * np.log(pc if pc > 0 else 1e-300), 0.0))
    P = np.where(ok, np.exp(lp), 0.0)
    if p10 == 0:
        P[a > 0] = 0.0
    if p01 == 0:
        P[b > 0] = 0.0
    if pc == 0:
        P[ok & (c > 0)] = 0.0
    return P


def power(N, test, alt, alpha, p10, p01, margin=0.0, gamma=0.0):
    return float((trinomial(N, p10, p01) * reject(N, test, alt, alpha, margin, gamma)).sum())


# ---------------------------------------------------------------------------------------
# Sample size
def p01_bar(p10, p01, d0):
    """Large-sample limit of the restricted estimate of p01 under p10 - p01 = d0, the
    maximizer of the expected log likelihood, by a root search."""
    pc = 1.0 - p10 - p01
    lo, hi = max(0.0, -d0), (1.0 - d0) / 2.0

    def score(p):
        s = 0.0
        if p10 > 0:
            s += p10 / (p + d0)
        if p01 > 0:
            s += p01 / p
        if pc > 0:
            s -= 2.0 * pc / (1.0 - 2.0 * p - d0)
        return s

    eps = 1e-15
    x0, x1 = lo + eps * (1 + abs(lo)), hi - eps * (1 + abs(hi))
    if score(x0) <= 0:
        return lo
    if score(x1) >= 0:
        return hi
    return brentq(score, x0, x1, xtol=1e-16, rtol=8.9e-16, maxiter=500)


def raw_n(p10, p01, alpha, tp, alt, margin=0.0, za=None, zb=None):
    if alt == "less":
        p10, p01 = p01, p10
        alt = "greater"
    d0 = -margin
    d1 = p10 - p01
    if za is None:
        za = norm.isf(alpha / 2 if alt == "two.sided" else alpha)
    if zb is None:
        zb = norm.ppf(tp)
    pb = p01_bar(p10, p01, d0)
    v0 = (pb + d0) + pb - d0 ** 2
    v1 = p10 + p01 - d1 ** 2
    dist = abs(d1 - d0) if alt == "two.sided" else d1 - d0
    return (za * math.sqrt(v0) + zb * math.sqrt(v1)) ** 2 / dist ** 2


def exact_n(p10, p01, alpha, tp, test, alt, margin=0.0):
    N = max(1, math.ceil(raw_n(p10, p01, alpha, tp, alt, margin) - 1e-9))
    pw = lambda n: power(n, test, alt, alpha, p10, p01, margin)
    P = pw(N)
    if P >= tp - TOL:
        while P >= tp - TOL and N > 1:
            N -= 1
            P = pw(N)
        if P < tp - TOL:
            N += 1
    else:
        while P < tp - TOL:
            N += 1
            P = pw(N)
    return N


def smallest_n(p10, p01, alpha, tp, test, alt, margin=0.0):
    N = 1
    while power(N, test, alt, alpha, p10, p01, margin) < tp - TOL:
        N += 1
    return N


# ---------------------------------------------------------------------------------------
# Blinded sample size re-estimation
def plan(D1, n1, DA):
    psi = max(D1 / n1, abs(DA))
    return psi, (psi + DA) / 2, (psi - DA) / 2


def cond_power_bf(N, n1, D1, psi, pi_, test, alt, alpha, margin=0.0):
    """Conditional power by enumeration of the interim directions and the second-stage
    multinomial outcomes."""
    R = reject(N, test, alt, alpha, margin)
    T = trinomial(N - n1, psi * pi_, psi * (1 - pi_))
    s = 0.0
    for a1 in range(D1 + 1):
        w = binom.pmf(a1, D1, pi_)
        if w == 0:
            continue
        sub = R[a1:a1 + N - n1 + 1, D1 - a1:D1 - a1 + N - n1 + 1]
        s += w * float((T * sub).sum())
    return s


def final_sizes(n1, DA, N0, alpha, tp, test, alt, rule, margin=0.0, Nmax=None,
                restricted=False):
    lo = max(n1, N0 if restricted else 0)
    out = []
    for D1 in range(n1 + 1):
        psi, p10, p01 = plan(D1, n1, DA)
        if rule == "standard":
            N = max(1, math.ceil(raw_n(p10, p01, alpha, tp, alt, margin) - 1e-9))
        elif rule == "exact":
            N = exact_n(p10, p01, alpha, tp, test, alt, margin)
        else:
            pi_ = p10 / psi if psi > 0 else 0.5
            N = Nmax
            for n in range(lo, Nmax + 1):
                if cond_power_bf(n, n1, D1, psi, pi_, test, alt, alpha, margin) >= tp - TOL:
                    N = n
                    break
        N = max(lo, N)
        if Nmax is not None:
            N = min(Nmax, N)
        out.append(N)
    return out


def bssr_power_bf(sizes, n1, test, alt, alpha, p10, p01, margin=0.0):
    """Rejection probability by enumeration of the interim and second-stage outcomes."""
    T1 = trinomial(n1, p10, p01)
    tot = 0.0
    for a1 in range(n1 + 1):
        for b1 in range(n1 + 1 - a1):
            w = T1[a1, b1]
            if w == 0:
                continue
            N = sizes[a1 + b1]
            R = reject(N, test, alt, alpha, margin)
            T2 = trinomial(N - n1, p10, p01)
            sub = R[a1:a1 + N - n1 + 1, b1:b1 + N - n1 + 1]
            tot += w * float((T2 * sub).sum())
    return tot


def expected_n(sizes, n1, psi):
    return float(sum(binom.pmf(D1, n1, psi) * sizes[D1] for D1 in range(n1 + 1)))


def boundary(psi, d0):
    return (psi + d0) / 2, (psi - d0) / 2


def sup_over(f, lo, hi, G=2001):
    grid = np.linspace(lo, hi, G)
    vals = np.array([f(t) for t in grid])
    best_i = int(np.argmax(vals))
    best = (vals[best_i], grid[best_i])
    for i in np.argsort(vals)[::-1][:4]:
        l, r = grid[max(0, i - 1)], grid[min(G - 1, i + 1)]
        o = minimize_scalar(lambda t: -f(t), bounds=(l, r), method="bounded",
                            options={"xatol": 1e-12})
        if -o.fun > best[0]:
            best = (-o.fun, o.x)
    return best


# ---------------------------------------------------------------------------------------
def blocks(emit, reported_level):
    # p-values of PairedRR at N = 12
    f = "test-PairedRR.R"
    cells = [(5, 1), (6, 2), (3, 0), (0, 0), (2, 3), (4, 4)]
    for key, (test, alt, margin, gamma) in {
        "McNemar greater": ("McNemar", "greater", 0.0, 0.0),
        "McNemar two-sided": ("McNemar", "two.sided", 0.0, 0.0),
        "exact greater": ("exact", "greater", 0.0, 0.0),
        "exact two-sided": ("exact", "two.sided", 0.0, 0.0),
        "midP greater": ("midP", "greater", 0.0, 0.0),
        "midP two-sided": ("midP", "two.sided", 0.0, 0.0),
        "uncond greater": ("uncond", "greater", 0.0, 0.0),
        "uncond two-sided": ("uncond", "two.sided", 0.0, 0.0),
        "uncond Berger-Boos": ("uncond", "greater", 0.0, 0.001),
        "score margin 0.1": ("Score", "greater", 0.1, 0.0),
        "score margin -0.05": ("Score", "greater", -0.05, 0.0),
        "score less margin 0.1": ("Score", "less", 0.1, 0.0),
        "score equivalence": ("Score", "equivalence", (0.1, 0.15), 0.0),
    }.items():
        P = pvalues_cached(12, test, alt, margin, gamma)
        rtol = 1e-9 if test == "uncond" else 1e-10
        emit(f, "p-values " + key, [P[a, b] for a, b in cells], rtol=rtol)
        emit(f, "rejected " + key, [int(reject(12, test, alt, 0.05, margin, gamma).sum())],
             rtol=0)

    # Score statistic at outcomes with zero counts and a negative margin
    f = "test-paired_score_z.R"
    Z = score_z(20, 0.1)
    emit(f, "score margin 0.1", [Z[0, 0], Z[0, 3], Z[4, 0], Z[7, 2], Z[1, 19]])
    Z = score_z(20, -0.15)
    emit(f, "score margin -0.15", [Z[0, 0], Z[0, 3], Z[4, 0], Z[7, 2], Z[19, 1]])

    # Power at N = 25
    f = "test-PairedPower.R"
    pts = [(0.3, 0.1), (0.25, 0.25), (0.05, 0.0), (0.6, 0.4)]
    for key, (test, alt, margin) in {
        "McNemar": ("McNemar", "greater", 0.0),
        "exact two-sided": ("exact", "two.sided", 0.0),
        "midP": ("midP", "greater", 0.0),
        "uncond": ("uncond", "greater", 0.0),
        "score margin 0.1": ("Score", "greater", 0.1),
    }.items():
        emit(f, "power " + key, [power(25, test, alt, 0.05, p, q, margin) for p, q in pts],
             rtol=1e-9 if test == "uncond" else 1e-10)
    emit(f, "power equivalence",
         [power(40, "Score", "equivalence", 0.05, p, q, (0.2, 0.2))
          for p, q in [(0.15, 0.15), (0.2, 0.1), (0.05, 0.0)]])

    # Sample size
    f = "test-PairedSampleSize.R"
    emit(f, "standard",
         [math.ceil(raw_n(0.3, 0.1, 0.025, 0.8, "greater") - 1e-9),
          math.ceil(raw_n(0.2, 0.1, 0.05, 0.9, "two.sided") - 1e-9),
          math.ceil(raw_n(0.1, 0.1, 0.025, 0.8, "greater", 0.1) - 1e-9),
          math.ceil(raw_n(0.1, 0.25, 0.05, 0.8, "less") - 1e-9)], rtol=0)
    emit("test-paired_raw_n.R", "raw", [raw_n(0.3, 0.1, 0.025, 0.8, "greater"),
                    raw_n(0.1, 0.1, 0.025, 0.8, "greater", 0.1),
                    raw_n(0.12, 0.1, 0.05, 0.8, "greater", -0.01)])
    emit(f, "exact", [exact_n(0.3, 0.1, 0.025, 0.8, "exact", "greater"),
                      exact_n(0.3, 0.1, 0.025, 0.8, "McNemar", "greater"),
                      exact_n(0.1, 0.1, 0.025, 0.8, "Score", "greater", 0.1)], rtol=0)
    emit(f, "smallest", [smallest_n(0.3, 0.1, 0.025, 0.8, "midP", "greater")], rtol=0)

    # Equivalence: root of the approximate power
    f = "test-paired_raw_n.R"
    za = norm.isf(0.05)
    def eq_root(p10, p01, lo, up, tp):
        d1 = p10 - p01
        v1 = p10 + p01 - d1 ** 2
        def pw(n):
            out = -1.0
            for d0, dist in ((-lo, d1 + lo), (up, up - d1)):
                if d0 < 0:
                    pb = p01_bar(p10, p01, d0)
                else:
                    pb = p01_bar(p01, p10, -d0)
                v0 = 2 * pb - abs(d0) - d0 ** 2 if d0 > 0 else 2 * pb + d0 - d0 ** 2
                out += norm.cdf((dist * math.sqrt(n) - za * math.sqrt(v0)) / math.sqrt(v1))
            return out - tp
        return brentq(pw, 1e-6, 1e6, xtol=1e-12)
    emit(f, "equivalence", [eq_root(0.15, 0.15, 0.2, 0.2, 0.8),
                            eq_root(0.2, 0.1, 0.2, 0.25, 0.9)])

    # Re-estimation designs: 10 interim pairs, 30 planned
    f = "test-PairedPowerBSSR.R"
    psis = [0.3, 0.5]
    designs = {
        "McNemar standard": dict(test="McNemar", alt="greater", DA=0.25, rule="standard",
                                 margin=0.0, Nmax=None),
        "McNemar exact": dict(test="McNemar", alt="greater", DA=0.25, rule="exact",
                              margin=0.0, Nmax=None),
        "exact conditional power": dict(test="exact", alt="two.sided", DA=0.25, rule="cp",
                                        margin=0.0, Nmax=45),
        "score standard": dict(test="Score", alt="greater", DA=0.0, rule="standard",
                               margin=0.2, Nmax=None),
    }
    sizes_of = {}
    for key, d in designs.items():
        sizes = final_sizes(10, d["DA"], 30, 0.05, 0.8, d["test"], d["alt"], d["rule"],
                            d["margin"], d["Nmax"])
        sizes_of[key] = sizes
        emit("test-paired_reestimate.R", "sizes " + key, sizes, rtol=0)
        DT = d["DA"]
        res = []
        for psi in psis:
            p10, p01 = (psi + DT) / 2, (psi - DT) / 2
            res += [bssr_power_bf(sizes, 10, d["test"], d["alt"], 0.05, p10, p01, d["margin"]),
                    expected_n(sizes, 10, psi)]
        emit(f, "power and E.N " + key, res)

    # Type I error rate of the McNemar design and of the score design
    f = "test-PairedTypeIErrorBSSR.R"
    for key, d0 in (("McNemar standard", 0.0), ("score standard", -0.2)):
        d = designs[key]
        sizes = sizes_of[key]
        lo = abs(d0)
        tb = lambda s: bssr_power_bf(sizes, 10, d["test"], d["alt"], 0.05,
                                     *boundary(s, d0), d["margin"])
        tt = lambda s: power(30, d["test"], d["alt"], 0.05, *boundary(s, d0), d["margin"])
        grid = [lo + (1 - lo) * k / 10 for k in range(11)]
        emit(f, "grid " + key, [tb(s) for s in grid[1:10:3]] + [tt(s) for s in grid[1:10:3]])
        mb = sup_over(tb, lo, 1.0, G=401)
        mt = sup_over(tt, lo, 1.0, G=401)
        emit(f, "certified max " + key, [mb[0], mt[0]], rtol=1e-9)

    # Adjusted level of the McNemar design
    f = "test-PairedAlphaAdjBSSR.R"
    d = designs["McNemar standard"]
    sizes = sizes_of["McNemar standard"]
    def max_rate(level, which):
        if which == "BSSR":
            g = lambda s: bssr_power_bf(sizes, 10, "McNemar", "greater", level, s / 2, s / 2)
        else:
            g = lambda s: power(30, "McNemar", "greater", level, s / 2, s / 2)
        return sup_over(g, 0.0, 1.0, G=401)[0]
    res = []
    for which in ("BSSR", "TRAD"):
        m0 = max_rate(0.05, which)
        if m0 <= 0.05:
            res.append((0.05, m0, m0))
            continue
        lo_, hi_ = 0.0, 0.05
        while hi_ - lo_ > 1e-9 * 0.05:
            mid = (lo_ + hi_) / 2
            if max_rate(mid, which) <= 0.05:
                lo_ = mid
            else:
                hi_ = mid
        Ns = sorted(set(sizes)) if which == "BSSR" else [30]
        pv = [pvalues_cached(n, "McNemar", "greater") for n in Ns]
        res.append((reported_level(pv, lo_), max_rate(lo_, which), m0))
    emit(f, "McNemar max at nominal level", [r[2] for r in res], rtol=1e-9)
    emit(f, "McNemar adjusted level", [r[0] for r in res], rtol=1e-12)
    emit(f, "McNemar max at adjusted level", [r[1] for r in res], rtol=1e-8)
