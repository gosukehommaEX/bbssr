# Independent check of the levels of the fixed designs of Figure 1 of Friede, Mitchell and
# Mueller-Velten (2007) that inst/reproduce/reproduce-figures.R stores in
# inst/extdata/published-figures/friede-2007-figure1.csv. The computation uses numpy and
# the standard library only, without the code of bbssr. For each design, the level of the
# test of Blackwelder or of Farrington and Manning at the one-sided level 0.025 with the
# margin 0.1 is the probability of the rejection region on the boundary of the null
# hypothesis at the assumed overall rate pa, with the response rates pa - 0.1 q / (1 + q)
# in group 1 (experimental) and pa + 0.1 / (1 + q) in group 2. The restricted maximum
# likelihood estimates of the test of Farrington and Manning are found by bisection on the
# score equation instead of the closed form of the article of Farrington and Manning
# (1990), and an outcome is rejected when its p-value is below 0.025 by more than 1.49e-8,
# as in bbssr.
# The script also prints the level of the design q = 1/3, pa = 0.7 for the splits of its
# total of 815 patients that the vignette 'Validation by Reproducing Published Figures'
# quotes.
# Run from the package root:
#   python3 tools/reference/check_friede_2007_figure1.py
# With the argument 'perturb', one stored level of each test is shifted by 1e-6, and the
# check must then fail.
import csv
import math
import sys
from statistics import NormalDist

import numpy as np

SRC = 'inst/extdata/published-figures/friede-2007-figure1.csv'
TOL = 1e-10
DELTA = 0.1
ALPHA = 0.025
# An outcome is rejected when its p-value is below ALPHA - 1.49e-8, that is, when z
# exceeds ZC
ZC = NormalDist().inv_cdf(1 - ALPHA + 1.49e-8)


def pmf(n, p):
    k = np.arange(n + 1)
    lg = np.array([math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
                   for i in k])
    return np.exp(lg + k * math.log(p) + (n - k) * math.log1p(-p))


def reject(N1, N2, test):
    h1 = (np.arange(N1 + 1) / N1)[:, None] * np.ones((1, N2 + 1))
    h2 = (np.arange(N2 + 1) / N2)[None, :] * np.ones((N1 + 1, 1))
    num = h1 - h2 + DELTA
    if test == 'Blackwelder':
        v = h1 * (1 - h1) / N1 + h2 * (1 - h2) / N2
    else:
        # Score equation of the likelihood with p1 = p2 - DELTA in p2, decreasing on
        # (DELTA, 1); bisection to machine precision
        lo = np.full(h1.shape, DELTA)
        hi = np.ones(h1.shape)
        for _ in range(80):
            mid = (lo + hi) / 2
            q1 = mid - DELTA
            with np.errstate(divide='ignore', invalid='ignore'):
                s = (N1 * (np.where(h1 > 0, h1 / q1, 0) -
                           np.where(h1 < 1, (1 - h1) / (1 - q1), 0)) +
                     N2 * (np.where(h2 > 0, h2 / mid, 0) -
                           np.where(h2 < 1, (1 - h2) / (1 - mid), 0)))
            up = s > 0
            lo = np.where(up, mid, lo)
            hi = np.where(up, hi, mid)
        q2 = (lo + hi) / 2
        q1 = q2 - DELTA
        v = q1 * (1 - q1) / N1 + q2 * (1 - q2) / N2
    with np.errstate(divide='ignore', invalid='ignore'):
        z = num / np.sqrt(v)
    zero = v <= 0
    z[zero] = np.where(num[zero] == 0, 0, np.sign(num[zero]) * np.inf)
    return (z > ZC).astype(float)


def level(N1, N2, pa, q, test):
    p1 = pa - DELTA * q / (1 + q)
    p2 = pa + DELTA / (1 + q)
    return pmf(N1, p1) @ reject(N1, N2, test) @ pmf(N2, p2)


with open(SRC, newline='') as f:
    rows = list(csv.DictReader(f))
if len(sys.argv) > 1 and sys.argv[1] == 'perturb':
    for test in ('Blackwelder', 'Farrington-Manning'):
        k = [i for i, r in enumerate(rows) if r['Test'] == test][0]
        rows[k]['fixed'] = str(float(rows[k]['fixed']) + 1e-6)
ok = True
for test in ('Blackwelder', 'Farrington-Manning'):
    d = [abs(level(int(float(r['N1'])), int(float(r['N2'])), float(r['pa']),
                   float(r['q']), test) - float(r['fixed']))
         for r in rows if r['Test'] == test]
    print('{}: {} fixed designs, largest difference {:.2e}'.format(test, len(d), max(d)))
    ok = ok and max(d) < TOL
for n1, n2 in ((611, 204), (610, 204), (610, 205), (612, 203)):
    print('q = 1/3, pa = 0.7, {} + {}: level {:.6f}'.format(
        n1, n2, level(n1, n2, 0.7, 1 / 3, 'Farrington-Manning')))
print('PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
