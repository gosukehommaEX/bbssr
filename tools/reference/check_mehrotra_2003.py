# Independent check of the items of Mehrotra, Chan and Berger (2003), Berger and Boos
# (1994) and Fay and Hunsberger (2021) in inst/reproduce/reproduce-published.R.
# The numbers are computed without the code of bbssr. The conditional p-values come from
# the hypergeometric tails of scipy (fisher_two_sided of reference_values.py). The tail
# probability of an unconditional test is the polynomial sum_s h[s] b(s; N, theta), whose
# maximum is taken on a grid and, for the outcomes whose grid maximum lies within 0.01
# below the level, certified by branch and bound on its Bernstein coefficients
# (tail_coefficients and certified_max of reference_values.py). Sizes are certified in the
# same way. The predicted verdicts follow the rules of the reproduction script, and the
# recomputed values of bbssr in reproduce-output/published-comparison.csv are compared
# with the predictions side by side. Run from the package root where scipy is available:
#   python3 tools/reference/check_mehrotra_2003.py [path to published-comparison.csv]
# With the argument 'perturb' in place of the path, three published values are shifted by
# one unit in their last digit, and the predicted verdicts must then include FAIL.
import csv
import os
import sys

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import beta, binom, hypergeom, norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reference_values import (certified_max, certified_size, fisher_two_sided,  # noqa: E402
                              tail_coefficients)

TOL = 1.4901161193847656e-08  # tolerance of fpCompare, used by BinaryRR()
GAM = 0.001
COLS = ['F', 'D', 'D*', 'B', 'B*', 'ZP', 'ZP*', 'ZU', 'ZU*', 'ZP asymptotic',
        'ZU asymptotic']


def cp_bounds(N, gam):
    s = np.arange(N + 1)
    lo = beta.ppf(gam / 2, np.maximum(s, 1), N - s + 1)
    up = beta.ppf(1 - gam / 2, s + 1, np.maximum(N - s, 1))
    lo[s == 0], up[s == N] = 0.0, 1.0
    return lo, up


def statistics(N1, N2):
    """Absolute values of the difference, of the pooled and of the unpooled Z statistic."""
    h1 = np.arange(N1 + 1)[:, None] / N1
    h2 = np.arange(N2 + 1)[None, :] / N2
    d = h1 - h2
    t = (np.arange(N1 + 1)[:, None] + np.arange(N2 + 1)[None, :]) / (N1 + N2)
    vp = t * (1 - t) * (1 / N1 + 1 / N2)
    vu = h1 * (1 - h1) / N1 + h2 * (1 - h2) / N2
    with np.errstate(divide='ignore', invalid='ignore'):
        zp = np.where(vp > 0, d / np.sqrt(vp), 0.0)
        zu = np.where(vu > 0, d / np.sqrt(vu), np.where(d == 0, 0.0, np.inf))
    zu = np.abs(zu)
    zu_exact = zu.copy()
    zu[np.isinf(zu)] = zu[np.isfinite(zu)].max() + 1
    return np.abs(d), np.abs(zp), zu, zu_exact


def unconditional(stat, N1, N2, decreasing, gam=0.0):
    """Unconditional p-values of every outcome: grid maximum, certified near the level."""
    N = N1 + N2
    H, tot = tail_coefficients(stat, N1, N2, decreasing)
    lo, up = cp_bounds(N, gam) if gam > 0 else (np.zeros(N + 1), np.ones(N + 1))
    phi = np.linspace(0, np.pi / 2, int(np.ceil(4 * np.pi * np.sqrt(N))) + 1)
    theta = np.unique(np.r_[np.linspace(0, 1, 2001), np.sin(phi) ** 2, lo, up])
    best = np.full(len(H), -np.inf)
    s = np.arange(N + 1)
    for k in range(0, len(theta), 200):
        th = theta[k:k + 200]
        val = H @ binom.pmf(s[:, None], N, th[None, :])
        inside = ((th[None, :] >= lo[tot][:, None] - 1e-15)
                  & (th[None, :] <= up[tot][:, None] + 1e-15))
        best = np.maximum(best, np.where(inside, val, -np.inf).max(axis=1))
    near = np.where((best + gam > 0.04) & (best + gam < 0.05))[0]
    memo = {}
    for k in near:
        key = (H[k].tobytes(), lo[tot[k]], up[tot[k]])
        if key not in memo:
            b, bound = certified_max(H[k], lo[tot[k]], up[tot[k]])
            assert bound - b < 1e-12
            memo[key] = b
        best[k] = memo[key]
    return np.minimum(1, best + gam).reshape(N1 + 1, N2 + 1)


def single(stat, N1, N2, decreasing, cell, gam=0.0):
    """Certified unconditional p-value of one outcome."""
    N = N1 + N2
    H, tot = tail_coefficients(stat, N1, N2, decreasing)
    k = cell[0] * (N2 + 1) + cell[1]
    lo, up = cp_bounds(N, gam) if gam > 0 else (np.zeros(N + 1), np.ones(N + 1))
    b, bound = certified_max(H[k], lo[tot[k]], up[tot[k]])
    assert bound - b < 1e-12
    return min(1.0, b + gam), H[k]


def pvalues(N1, N2, ts):
    d, zp, zu, zu_exact = statistics(N1, N2)
    f = fisher_two_sided(N1, N2, ts)
    return {'F': f, 'D': unconditional(d, N1, N2, True),
            'D*': unconditional(d, N1, N2, True, GAM),
            'B': unconditional(f, N1, N2, False), 'B*': unconditional(f, N1, N2, False, GAM),
            'ZP': unconditional(zp, N1, N2, True), 'ZP*': unconditional(zp, N1, N2, True, GAM),
            'ZU': unconditional(zu, N1, N2, True), 'ZU*': unconditional(zu, N1, N2, True, GAM),
            'ZP asymptotic': np.minimum(1, 2 * norm.sf(zp)),
            'ZU asymptotic': np.minimum(1, 2 * norm.sf(zu_exact))}


def rate(rr, N1, N2, t1, t2):
    return float(binom.pmf(np.arange(N1 + 1), N1, t1) @ rr.astype(float)
                 @ binom.pmf(np.arange(N2 + 1), N2, t2))


def round2(x):
    return (np.floor(x * 100 + 0.5 + 1e-9) + 5) // 10 / 10


def round1(x):
    return np.floor(x * 10 + 0.5 + 1e-9) / 10


# Published values. Table 1: for each (N1, N2), the rows theta = 0.02, 0.10, 0.25, 0.50
# and the size, each with the eleven columns of COLS. Table 3: theta1, theta2 and the
# powers of the nine exact tests
T1 = {
    (10, 10): [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1,
               0.1, 0.1, 0.1, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 5.0,
               1.0, 1.8, 1.8, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 9.4,
               1.3, 4.1, 4.1, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 8.8,
               1.3, 4.1, 4.1, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 9.5],
    (25, 25): [0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
               0.1, 0.1, 0.1, 1.9, 1.9, 3.9, 3.9, 3.9, 3.9, 4.9, 4.9,
               2.2, 1.4, 1.4, 4.1, 4.1, 4.2, 4.2, 4.2, 4.2, 5.4, 6.4,
               3.3, 3.3, 3.3, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 6.5, 6.5,
               3.3, 3.3, 3.3, 4.6, 4.6, 4.6, 4.6, 4.6, 4.6, 6.5, 6.5],
    (50, 50): [0.0, 0.0, 0.0, 0.2, 0.2, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3,
               1.8, 0.1, 0.3, 3.2, 3.6, 3.8, 4.8, 3.8, 4.8, 5.1, 5.9,
               3.0, 1.5, 1.6, 4.1, 4.1, 4.5, 4.7, 4.5, 4.7, 5.1, 5.7,
               3.5, 3.5, 3.5, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 5.7, 5.7,
               3.5, 3.5, 3.5, 4.9, 4.9, 4.9, 4.9, 4.9, 4.9, 5.7, 6.1],
    (150, 150): [1.2, 0.0, 0.0, 2.3, 2.9, 4.6, 4.6, 4.6, 4.6, 4.6, 4.6,
                 3.3, 0.1, 1.0, 3.9, 4.5, 4.8, 4.8, 4.8, 4.8, 5.0, 5.2,
                 3.7, 2.0, 2.7, 4.6, 4.8, 4.8, 4.8, 4.8, 4.8, 5.0, 5.2,
                 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 5.7, 5.7,
                 4.3, 4.3, 4.3, 5.0, 4.9, 5.0, 4.9, 5.0, 4.9, 5.7, 5.7],
    (16, 4): [0.2, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 5.7, 0.2,
              1.2, 0.3, 0.3, 2.9, 2.9, 2.9, 2.9, 0.0, 0.0, 8.3, 5.8,
              1.5, 1.3, 1.3, 3.7, 3.7, 3.7, 3.7, 0.4, 0.4, 4.3, 22.4,
              1.4, 3.0, 3.0, 3.4, 3.4, 3.4, 3.4, 2.8, 2.8, 5.6, 14.3,
              1.5, 3.0, 3.0, 3.9, 3.9, 3.9, 3.9, 2.8, 2.8, 9.0, 22.8],
    (40, 10): [0.8, 0.0, 0.0, 0.8, 0.8, 1.3, 1.3, 0.0, 0.0, 8.8, 0.7,
               2.4, 0.2, 0.4, 2.6, 2.6, 3.9, 3.9, 0.6, 0.6, 4.5, 20.8,
               3.4, 1.6, 1.6, 4.3, 4.3, 4.1, 4.1, 4.1, 4.1, 4.7, 9.5,
               2.9, 3.9, 3.9, 3.5, 3.5, 4.1, 4.6, 1.5, 1.5, 5.5, 8.9,
               3.4, 3.9, 3.9, 4.8, 4.8, 4.3, 4.7, 4.5, 4.5, 9.1, 21.2],
    (80, 20): [1.5, 0.0, 0.0, 1.6, 1.6, 3.3, 3.3, 0.0, 0.0, 8.7, 5.2,
               2.7, 3.4, 0.2, 3.7, 3.7, 3.2, 3.2, 3.4, 3.4, 3.7, 13.0,
               3.3, 1.7, 2.1, 4.6, 4.3, 4.6, 4.8, 1.2, 2.0, 5.1, 6.9,
               4.3, 4.0, 4.0, 4.3, 4.3, 4.1, 4.6, 0.6, 4.1, 5.1, 6.6,
               4.3, 4.0, 4.0, 5.0, 4.7, 4.6, 4.9, 3.9, 4.1, 9.2, 21.3],
    (240, 60): [1.9, 0.0, 0.2, 2.4, 3.2, 4.0, 4.0, 0.7, 0.7, 4.3, 21.3,
                3.5, 0.1, 1.2, 3.9, 4.9, 4.4, 4.7, 1.2, 2.5, 4.7, 7.0,
                4.1, 2.1, 2.9, 4.6, 4.7, 4.4, 4.8, 0.5, 4.6, 4.9, 5.7,
                4.3, 4.6, 4.6, 4.3, 4.3, 4.6, 4.6, 0.4, 4.9, 5.3, 5.5,
                4.3, 4.6, 4.6, 4.9, 4.9, 4.6, 4.9, 4.3, 4.9, 9.2, 21.4],
}
T3 = {
    (10, 10): [(0.02, 0.54, [62.8, 66.9, 66.9, 80.8, 80.8, 80.8, 80.8, 80.8, 80.8]),
               (0.10, 0.68, [62.9, 77.7, 77.7, 79.4, 79.4, 79.4, 79.4, 79.4, 79.4]),
               (0.25, 0.84, [61.9, 78.9, 78.9, 79.2, 79.2, 79.2, 79.2, 79.2, 79.2]),
               (0.50, 0.99, [57.9, 59.9, 59.9, 78.4, 78.4, 78.4, 78.4, 78.4, 78.4])],
    (25, 25): [(0.02, 0.29, [68.1, 36.8, 36.8, 76.4, 76.4, 80.4, 80.4, 80.4, 80.4]),
               (0.10, 0.45, [72.3, 66.9, 66.9, 80.9, 80.9, 80.9, 80.9, 80.9, 80.9]),
               (0.25, 0.65, [78.4, 78.4, 78.4, 80.4, 80.4, 80.4, 80.4, 80.4, 80.4]),
               (0.50, 0.86, [70.8, 69.3, 69.3, 79.7, 79.7, 79.7, 79.7, 79.7, 79.7])],
    (50, 50): [(0.02, 0.18, [68.6, 19.1, 39.6, 77.8, 78.6, 79.5, 82.2, 79.5, 82.2]),
               (0.10, 0.33, [75.7, 60.0, 65.3, 80.1, 80.3, 80.7, 82.1, 80.7, 82.1]),
               (0.25, 0.52, [74.5, 74.1, 74.1, 80.1, 80.1, 80.1, 80.1, 80.1, 80.1]),
               (0.50, 0.77, [75.3, 74.3, 74.3, 80.7, 80.7, 80.8, 80.8, 80.8, 80.8])],
    (150, 150): [(0.02, 0.09, [70.4, 4.0, 43.1, 75.0, 77.7, 79.4, 79.4, 79.4, 79.4]),
                 (0.10, 0.22, [77.6, 53.1, 69.0, 80.6, 81.5, 81.5, 81.5, 81.5, 81.5]),
                 (0.25, 0.40, [76.1, 73.4, 73.8, 78.9, 78.5, 78.9, 78.9, 78.9, 78.9]),
                 (0.50, 0.66, [78.0, 78.0, 78.0, 80.0, 79.7, 80.0, 79.7, 80.0, 79.7])],
    (16, 4): [(0.02, 0.54, [64.2, 37.4, 37.4, 73.0, 73.0, 73.0, 73.0, 8.5, 8.5]),
              (0.10, 0.68, [58.3, 53.1, 53.1, 73.5, 73.5, 73.5, 73.5, 21.4, 21.4]),
              (0.25, 0.84, [47.9, 53.3, 53.3, 61.9, 61.9, 61.9, 61.9, 45.8, 45.8]),
              (0.50, 0.99, [10.1, 21.8, 21.8, 21.9, 21.9, 21.9, 21.9, 21.8, 21.8])],
    (4, 16): [(0.02, 0.54, [16.3, 31.0, 31.0, 31.2, 31.2, 31.2, 31.2, 31.0, 31.0]),
              (0.10, 0.68, [41.0, 52.9, 52.9, 56.6, 56.6, 56.6, 56.6, 50.8, 50.8]),
              (0.25, 0.84, [53.7, 53.2, 53.2, 68.4, 68.4, 68.4, 68.4, 31.4, 31.4]),
              (0.50, 0.99, [63.2, 31.2, 31.2, 68.3, 68.3, 68.3, 68.3, 6.3, 6.3])],
    (40, 10): [(0.02, 0.29, [68.7, 28.8, 31.5, 68.9, 68.9, 77.6, 77.6, 4.0, 4.0]),
               (0.10, 0.45, [67.3, 46.6, 50.0, 71.3, 71.3, 71.8, 71.8, 16.6, 16.6]),
               (0.25, 0.65, [60.0, 59.8, 59.8, 68.7, 68.7, 65.7, 66.9, 29.1, 29.1]),
               (0.50, 0.86, [50.2, 51.7, 51.7, 54.0, 54.0, 56.1, 56.4, 42.9, 42.9])],
    (10, 40): [(0.02, 0.29, [30.3, 12.9, 12.9, 30.5, 30.5, 30.5, 30.5, 70.4, 70.4]),
               (0.10, 0.45, [51.2, 48.6, 48.6, 56.1, 56.1, 56.8, 56.8, 47.1, 47.1]),
               (0.25, 0.65, [55.6, 60.9, 60.9, 60.9, 60.9, 61.9, 65.2, 36.8, 36.8]),
               (0.50, 0.86, [64.1, 49.5, 50.1, 68.7, 68.7, 68.5, 68.5, 18.5, 18.5])],
    (80, 20): [(0.02, 0.18, [64.6, 12.9, 38.5, 70.6, 70.6, 76.2, 76.2, 2.4, 2.4]),
               (0.10, 0.33, [65.0, 39.9, 51.1, 68.1, 68.1, 68.4, 68.4, 10.6, 10.6]),
               (0.25, 0.52, [57.4, 54.6, 55.7, 64.8, 63.3, 63.9, 63.9, 18.4, 41.8]),
               (0.50, 0.77, [59.5, 56.2, 56.2, 59.7, 59.7, 57.7, 59.9, 32.1, 60.7])],
    (20, 80): [(0.02, 0.18, [32.4, 2.9, 5.0, 40.6, 39.8, 40.7, 40.7, 62.7, 62.7]),
               (0.10, 0.33, [49.4, 39.5, 42.1, 56.5, 55.7, 53.2, 56.9, 42.6, 57.1]),
               (0.25, 0.52, [58.6, 56.0, 56.0, 58.6, 58.6, 56.9, 59.0, 30.1, 59.0]),
               (0.50, 0.77, [59.3, 54.5, 56.3, 65.9, 64.8, 65.7, 65.7, 18.2, 37.8])],
    (240, 60): [(0.02, 0.09, [64.1, 3.4, 38.2, 67.3, 68.5, 68.9, 69.2, 2.8, 2.8]),
                (0.10, 0.22, [62.7, 33.2, 53.4, 65.1, 67.4, 67.8, 67.8, 13.4, 41.4]),
                (0.25, 0.40, [59.8, 53.4, 56.2, 61.1, 61.2, 61.1, 62.3, 18.1, 54.5]),
                (0.50, 0.66, [59.2, 59.6, 59.6, 59.2, 59.2, 59.7, 60.0, 25.3, 61.5])],
    (60, 240): [(0.02, 0.09, [41.1, 0.1, 8.5, 41.1, 48.5, 43.2, 44.6, 45.0, 46.8]),
                (0.10, 0.22, [54.1, 31.5, 43.5, 57.5, 57.8, 55.5, 58.3, 37.0, 66.2]),
                (0.25, 0.40, [56.0, 54.5, 54.6, 58.6, 58.6, 57.5, 58.5, 27.1, 61.8]),
                (0.50, 0.66, [58.7, 59.0, 59.0, 62.7, 62.0, 61.3, 62.2, 21.0, 59.0])],
}


items = []  # (item, published, recomputed, digits, rule, tolerance of the comparison)
near_count = 0
# Section 3.1: 8 of 148 against 1 of 132, with the blaker convention for F, B and B*
lab = 'MCB2003 Section 3.1: 8 / 148 against 1 / 132'
lo, up = cp_bounds(280, GAM)
items += [(lab + ': 99.9% confidence interval, lower', 0.0080, lo[9], 4, 'round', 1e-10),
          (lab + ': 99.9% confidence interval, upper', 0.0826, up[9], 4, 'round', 1e-10)]
d, zp, zu, _ = statistics(148, 132)
fb = fisher_two_sided(148, 132, 'blaker')
ex = {'F': fb[8, 1]}
for name, stat, dec in [('D', d, True), ('B', fb, False), ('ZP', zp, True), ('ZU', zu, True)]:
    ex[name] = single(stat, 148, 132, dec, (8, 1))[0]
    ex[name + '*'] = single(stat, 148, 132, dec, (8, 1), GAM)[0]
pub = [0.0388, 0.4386, 0.1603, 0.0347, 0.0325, 0.0291, 0.0282, 0.0229, 0.0215]
for k, c in enumerate(COLS[:9]):
    items.append((f'{lab}: two-sided p-value, {c}', pub[k], ex[c], 4, 'round', 1e-8))
# Tables 1 and 3
tab = []
for key in list(T1) + [k for k in T3 if k not in T1]:
    N1, N2 = key
    P = pvalues(N1, N2, 'minlike')
    near_count += sum(int((np.abs(p - 0.05) < 1e-6).sum()) for p in P.values())
    R = {c: p < 0.05 - TOL for c, p in P.items()}
    if N1 != N2:
        f = fisher_two_sided(N1, N2, 'blaker')
        Rb = {'F': f < 0.05 - TOL, 'B': unconditional(f, N1, N2, False) < 0.05 - TOL,
              'B*': unconditional(f, N1, N2, False, GAM) < 0.05 - TOL}
    else:
        Rb = {c: R[c] for c in ['F', 'B', 'B*']}
    if key in T1:
        rows = [f'theta = {t:.2f}' for t in (0.02, 0.10, 0.25, 0.50)] + ['size']
        for j, c in enumerate(COLS):
            for i, r in enumerate(rows):
                v = (rate(R[c], N1, N2, [0.02, 0.10, 0.25, 0.50][i], [0.02, 0.10, 0.25, 0.50][i])
                     if i < 4 else certified_size(R[c], N1, N2)) * 100
                tab.append((f'MCB2003 Table 1: ({N1}, {N2}), {r}: {c}', T1[key][i * 11 + j], v))
    if key in T3:
        for (t1, t2, vals) in T3[key]:
            for j, c in enumerate(COLS[:9]):
                v = rate(R[c], N1, N2, t1, t2) * 100
                name = f'MCB2003 Table 3: ({N1}, {N2}), ({t1:.2f}, {t2:.2f}): {c}'
                tab.append((name, vals[j], v))
    print(f'({N1}, {N2}) done', flush=True)
lost = 0
for name, p, v in tab:
    items.append((name, p, v, 1, 'round2', 1e-6))
    lost += int(abs(round1(v) - p) < 1e-9 and abs(round2(v) - p) > 1e-9)
items.append(('MCB2003 Tables 1 and 3: p-values within 1e-6 of the level 0.05', None,
              near_count, None, 'info', 0))
items.append(('MCB2003 Tables 1 and 3: items that ordinary rounding reproduces and the rule '
              'round2 does not', None, lost, None, 'info', 0))
# Berger and Boos (1994), Example 2
lab = 'BB1994 Example 2: 14 / 47 against 48 / 283'
d, zp, zu, _ = statistics(47, 283)
lo, up = cp_bounds(330, 0.001)
psup, h = single(zp, 47, 283, True, (14, 48))
pbb = single(zp, 47, 283, True, (14, 48), 0.001)[0]
# Location of the maximum in [0, 1/2]: grid of step 1e-4 refined by a bounded search
th = np.round(np.arange(0, 5001) * 1e-4, 10)
tail = binom.pmf(np.arange(331)[:, None], 330, th[None, :]).T @ h
kmax = int(np.argmax(tail))
loc = minimize_scalar(lambda t: -float(binom.pmf(np.arange(331), 330, t) @ h),
                      bounds=(th[max(0, kmax - 1)], th[min(5000, kmax + 1)]),
                      method='bounded', options={'xatol': 1e-12}).x
items += [(lab + ': chi-squared statistic', 4.346, zp[14, 48] ** 2, 3, 'round', 1e-10),
          (lab + ': .999 confidence interval, lower', 0.123, lo[62], 3, 'round', 1e-10),
          (lab + ': .999 confidence interval, upper', 0.267, up[62], 3, 'round', 1e-10),
          (lab + ': p-value maximized over [0, 1]', 0.061, psup, 3, 'truncate', 1e-8),
          (lab + ': maximum over the confidence interval', 0.036, pbb - 0.001, 3, 'truncate',
           1e-8),
          (lab + ': p-value p_.001', 0.037, pbb, 3, 'truncate', 1e-8),
          (lab + ': location of the maximum in [0, 1/2]', 0.003, loc, None, 'info', 1e-6)]
# Fay and Hunsberger (2021), Section 8 and Table 1
fh = {ts: fisher_two_sided(14, 7, ts) for ts in ['blaker', 'minlike', 'central']}
for ts, p in zip(fh, [0.087, 0.159, 0.157]):
    items.append((f'FH2021 Section 8: 8 / 14 against 1 / 7: two-sided p-value, {ts}', p,
                  fh[ts][8, 1], 3, 'round', 1e-10))
for x2, p in enumerate([0.007, 0.087, 0.642, 1.000, 0.397, 0.159, 0.016, 0.000]):
    items.append((f'FH2021 Table 1: T_B(x, 1) at x2 = {x2}', p, fh['blaker'][9 - x2, x2], 3,
                  'round', 1e-10))

if len(sys.argv) > 1 and sys.argv[1] == 'perturb':  # self-test
    for k, shift in [(2, 1e-4), (20, 0.1), (len(items) - 5, 1e-3)]:
        n, p, v, dg, rule, tol = items[k]
        items[k] = (n, p + shift, v, dg, rule, tol)


def verdict(p, v, digits, rule):
    if rule == 'info':
        return 'INFO'
    if rule == 'truncate':
        shown = np.floor(v * 10 ** digits + 1e-9) / 10 ** digits
    elif rule == 'round2':
        shown = (np.floor(v * 10 ** (digits + 1) + 0.5 + 1e-9) + 5) // 10 / 10 ** digits
    else:
        shown = round(v, digits)
    return 'PASS' if abs(shown - p) < 1e-9 else 'FAIL'


results = {}
counts = {}
for n, p, v, dg, rule, tol in items:
    vd = verdict(p, v, dg, rule)
    results[n] = (vd, p, v, tol)
    counts[vd] = counts.get(vd, 0) + 1
    if vd != 'PASS':
        print(f'{vd:5s} | {n} | published {p} | recomputed {v:.10g}')
print('predicted items:', len(items), counts)

SOURCES = ['Mehrotra, Chan and Berger (2003)', 'Berger and Boos (1994)',
           'Fay and Hunsberger (2021)']


def compare(path):
    rows = [r for r in csv.DictReader(open(path, encoding='utf-8')) if r['source'] in SOURCES]
    bad = 0
    seen = set()
    for r in rows:
        name = r['item']
        if name not in results:
            print('NOT PREDICTED', name)
            bad += 1
            continue
        seen.add(name)
        vd, p, v, tol = results[name]
        rr = float(r['recomputed'])
        pub_ok = (p is None and r['published'] in ('', 'NA')) or \
            (p is not None and abs(float(r['published']) - p) < 1e-12)
        ok = abs(rr - v) <= tol and r['verdict'] == vd and pub_ok
        if not ok:
            bad += 1
            print(f'MISMATCH | {name} | bbssr {rr:.12g} {r["verdict"]} | python {v:.12g} {vd}')
    missing = set(results) - seen
    for name in sorted(missing):
        print('NOT IN CSV', name)
    print(f'items in the CSV: {len(rows)}, predictions: {len(results)}, '
          f'mismatches: {bad + len(missing)}')
    return bad + len(missing)


if len(sys.argv) > 1 and sys.argv[1] != 'perturb':
    sys.exit(1 if compare(sys.argv[1]) else 0)
