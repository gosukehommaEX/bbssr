# Independent check of the values of Figure 21.1 of Kieser (2020) that
# inst/reproduce/reproduce-figures.R stores in
# inst/extdata/published-figures/kieser-2020-figure21-1.csv. The computation uses numpy and
# the standard library only, without the code of bbssr:
#   1. For the designs with an initial total of at most 100 patients, the type I error rate
#      of the recalculation design under the rule of bbssr (formula (21.5) of the book with
#      the Bernoulli variances truncated at zero, groups rounded up) is computed from
#      scratch and compared with the column ips.
#   2. For every design, the change of the type I error rate when the trial stops after
#      the pilot for the interim outcomes with a recovered rate outside the unit interval
#      is computed and compared with ips.stop - ips.
#   3. For the design of Example 21.1 (pilot of 79 + 79, |Delta| = 0.15), the same change
#      must vanish over the overall rates of Figure 21.2 a, [0.075, 0.925]; with r = 1,
#      exchanging the labels of the groups turns pE < pC with Delta = -0.15 into pE > pC
#      with Delta = 0.15.
# Run from the package root:
#   python3 tools/reference/check_kieser_2020_pilot_stop.py
# With the argument 'perturb', one stored value of each column is shifted by 1e-6, and the
# check must then fail.
import csv
import math
import sys
from statistics import NormalDist

import numpy as np

SRC = 'inst/extdata/published-figures/kieser-2020-figure21-1.csv'
TOL = 1e-10
ZA = NormalDist().inv_cdf(0.975)
ZB = NormalDist().inv_cdf(0.8)
# bbssr rejects when the one-sided p-value is below 0.025 - 1.49e-8
ZC = NormalDist().inv_cdf(1 - 0.025 + 1.49e-8)


def pmf(n, p):
    k = np.arange(n + 1)
    lg = np.array([math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
                   for i in k])
    return np.exp(lg + k * math.log(p) + (n - k) * math.log1p(-p))


_rr = {}


def reject(N1, N2):
    """Rejection region of the one-sided normal approximation test, E (group 1) > C."""
    if (N1, N2) not in _rr:
        x1 = np.arange(N1 + 1)[:, None]
        x2 = np.arange(N2 + 1)[None, :]
        p = (x1 + x2) / (N1 + N2)
        with np.errstate(divide='ignore', invalid='ignore'):
            u = math.sqrt(N1 * N2 / (N1 + N2)) * (x1 / N1 - x2 / N2) / np.sqrt(p * (1 - p))
        _rr[(N1, N2)] = np.where(np.isfinite(u), u > ZC, False).astype(float)
    return _rr[(N1, N2)]


def final_sizes(s, n1E, n1C, r, D):
    """Final sizes under the rule of bbssr, and whether a recovered rate is outside [0, 1]."""
    ph = s / (n1E + n1C)
    pE = ph + D / (1 + r)
    pC = ph - r * D / (1 + r)
    v0 = max(ph * (1 - ph), 0) * (1 + 1 / r)
    v1 = max(pE * (1 - pE), 0) / r + max(pC * (1 - pC), 0)
    n2 = (ZA * math.sqrt(v0) + ZB * math.sqrt(v1)) ** 2 / (pE - pC) ** 2
    N2 = max(n1C, math.ceil(n2 - 1e-9))
    N1 = max(n1E, math.ceil(r * N2))
    return N1, N2, pE > 1 + 1e-12 or pC < -1e-12


def cond_reject(a, b, N1, N2, n1E, n1C, p):
    m1 = N1 - n1E
    m2 = N2 - n1C
    return pmf(m1, p) @ reject(N1, N2)[a:a + m1 + 1, b:b + m2 + 1] @ pmf(m2, p)


def level_bbssr(p, D, r, n1E, n1C):
    pa = pmf(n1E, p)
    pb = pmf(n1C, p)
    tot = 0.0
    for s in range(n1E + n1C + 1):
        N1, N2, _ = final_sizes(s, n1E, n1C, r, D)
        for a in range(max(0, s - n1C), min(n1E, s) + 1):
            tot += pa[a] * pb[s - a] * cond_reject(a, s - a, N1, N2, n1E, n1C, p)
    return tot


def change_pilot_stop(p, D, r, n1E, n1C):
    pa = pmf(n1E, p)
    pb = pmf(n1C, p)
    R0 = reject(n1E, n1C)
    d = 0.0
    for s in range(n1E + n1C + 1):
        N1, N2, out = final_sizes(s, n1E, n1C, r, D)
        if not out:
            continue
        for a in range(max(0, s - n1C), min(n1E, s) + 1):
            b = s - a
            d += pa[a] * pb[b] * (R0[a, b] - cond_reject(a, b, N1, N2, n1E, n1C, p))
    return d


with open(SRC, newline='') as f:
    rows = list(csv.DictReader(f))
if len(sys.argv) > 1 and sys.argv[1] == 'perturb':
    small = [i for i, row in enumerate(rows) if float(row['nE']) + float(row['nC']) <= 100]
    rows[small[0]]['ips'] = str(float(rows[small[0]]['ips']) + 1e-6)
    rows[0]['ips.stop'] = str(float(rows[0]['ips.stop']) + 1e-6)
d1 = []
d2 = []
for row in rows:
    r = float(row['r'])
    D = float(row['Delta'])
    p = float(row['pA'])
    nE, nC, n1E, n1C = (int(float(row[k])) for k in ('nE', 'nC', 'n1E', 'n1C'))
    ips = float(row['ips'])
    if nE + nC <= 100:
        d1.append(abs(level_bbssr(p, D, r, n1E, n1C) - ips))
    d2.append(abs(change_pilot_stop(p, D, r, n1E, n1C) - (float(row['ips.stop']) - ips)))
print('rule of bbssr: {} designs, largest difference {:.2e}'.format(len(d1), max(d1)))
print('stop after the pilot: {} designs, largest difference {:.2e}'.format(len(d2), max(d2)))
d3 = [abs(change_pilot_stop(round(0.075 + 0.005 * k, 3), 0.15, 1, 79, 79))
      for k in range(171)]
print('Example 21.1: {} overall rates, largest change {:.2e}'.format(len(d3), max(d3)))
ok = max(d1) < TOL and max(d2) < TOL and max(d3) < TOL
print('PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
