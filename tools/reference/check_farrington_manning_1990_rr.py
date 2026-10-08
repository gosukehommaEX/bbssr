# Independent check of the values for the relative risk of Farrington and Manning (1990)
# that inst/reproduce/reproduce-published.R stores in reproduce-output/:
# farrington-manning-1990-table1-rr.csv (Table I), farrington-manning-1990-table2-rr.csv
# (Table II, Method 1) and farrington-manning-1990-example2-rr.csv (the second example).
# The computation uses numpy and the standard library only, without the code of bbssr.
# The restricted maximum likelihood estimates under p1 = R0 p2 are found by bisection on
# the score equation instead of the closed form of formula (13) of the article, the
# binomial probabilities come from math.lgamma, and an outcome is rejected when its
# p-value is below the level by more than 1.49e-8, as in bbssr.
# For Table I the script recomputes the sample sizes of formula (8) with the standard
# normal quantiles and with the quantiles 1.645 and 1.282, and the true powers at the
# published sample sizes, which bbssr gives in the columns N1, N2, N1.tab, N2.tab and
# power. For Table II it recomputes the sample sizes of Method 1, and for the second
# example the sample size of formula (8) rounded up and the true power, summed over the
# responder counts 0 to 450 of each group. It also prints how many published values each
# computation gives.
# Run from the package root after inst/reproduce/reproduce-published.R:
#   python3 tools/reference/check_farrington_manning_1990_rr.py
# With the argument 'perturb', one stored power is shifted by 1e-6 and one stored sample
# size by one patient, and the check must then fail.
import csv
import math
import sys
from statistics import NormalDist

import numpy as np

DIR = 'reproduce-output'
TOL = 1e-10
EPS = 1.49e-8
Z = NormalDist()


def read(name):
    with open(DIR + '/' + name, newline='') as f:
        return list(csv.DictReader(f))


def restricted(x1, n1, x2, n2, R0):
    """Estimates under p1 = R0 p2 by bisection on the score of the likelihood in p2,
    which is decreasing on (0, min(1, 1 / R0))."""
    x1, x2 = np.broadcast_arrays(np.asarray(x1, float), np.asarray(x2, float))
    lo = np.zeros(x1.shape)
    hi = np.full(x1.shape, min(1.0, 1.0 / R0))
    for _ in range(64):
        t = (lo + hi) / 2
        with np.errstate(divide='ignore', invalid='ignore'):
            s = (np.where(x1 > 0, x1 / t, 0) -
                 np.where(n1 - x1 > 0, (n1 - x1) * R0 / (1 - R0 * t), 0) +
                 np.where(x2 > 0, x2 / t, 0) -
                 np.where(n2 - x2 > 0, (n2 - x2) / (1 - t), 0))
        up = s > 0
        lo = np.where(up, t, lo)
        hi = np.where(up, hi, t)
    t = (lo + hi) / 2
    return R0 * t, t


def pmf(N, p, K=None):
    K = N if K is None else K
    lg = math.lgamma(N + 1)
    return np.array([math.exp(lg - math.lgamma(k + 1) - math.lgamma(N - k + 1) +
                              k * math.log(p) + (N - k) * math.log1p(-p))
                     for k in range(K + 1)])


def formula8(p1, p2, theta, R0, za, zb):
    """Size of group 1 from formula (8), with theta = N2 / N1. The large sample values
    are the estimates with the expected counts, scaled by 1e6, which is exact because
    the score is linear in the counts."""
    t1, t2 = restricted(p1 * 1e6, 1e6, p2 * theta * 1e6, theta * 1e6, R0)
    t1, t2 = float(t1), float(t2)
    v0 = t1 * (1 - t1) + R0 ** 2 / theta * t2 * (1 - t2)
    v1 = p1 * (1 - p1) + R0 ** 2 / theta * p2 * (1 - p2)
    return (za * math.sqrt(v0) + zb * math.sqrt(v1)) ** 2 / (p1 - R0 * p2) ** 2


def power(N1, N2, p1, p2, R0, alpha, lower, K1=None, K2=None):
    """Exact power over the responder counts 0 to K1 and 0 to K2. For lower = True the
    null hypothesis p1 / p2 >= R0 is rejected for a small statistic."""
    K1 = N1 if K1 is None else K1
    K2 = N2 if K2 is None else K2
    x1 = np.arange(K1 + 1)[:, None] + 0 * np.arange(K2 + 1)[None, :]
    x2 = np.arange(K2 + 1)[None, :] + 0 * np.arange(K1 + 1)[:, None]
    h1, h2 = x1 / N1, x2 / N2
    t1, t2 = restricted(x1, N1, x2, N2, R0)
    v = t1 * (1 - t1) / N1 + R0 ** 2 * t2 * (1 - t2) / N2
    num = h1 - R0 * h2
    with np.errstate(divide='ignore', invalid='ignore'):
        z = np.where(v > 0, num / np.sqrt(np.where(v > 0, v, 1)),
                     np.where(num == 0, 0, np.sign(num) * np.inf))
    rej = z < Z.inv_cdf(alpha - EPS) if lower else z > Z.inv_cdf(1 - alpha + EPS)
    return float(pmf(N1, p1, K1) @ rej.astype(float) @ pmf(N2, p2, K2))


t1 = read('farrington-manning-1990-table1-rr.csv')
t2 = read('farrington-manning-1990-table2-rr.csv')
ex = read('farrington-manning-1990-example2-rr.csv')[0]
if len(sys.argv) > 1 and sys.argv[1] == 'perturb':
    t1[0]['power'] = str(float(t1[0]['power']) + 1e-6)
    t2[0]['N1'] = str(float(t2[0]['N1']) + 1)
ok = True
za, zb = Z.inv_cdf(0.95), Z.inv_cdf(0.9)
d_n, d_p = 0, 0.0
pub_n = pub_tab = pub_pw = 0
for r in t1:
    p1, p2, R0, th = (float(r[k]) for k in ('p1', 'p2', 'R0', 'theta'))
    n1 = formula8(p1, p2, th, R0, za, zb)
    n1t = formula8(p1, p2, th, R0, 1.645, 1.282)
    N = (math.floor(n1 + 0.5), math.floor(th * n1 + 0.5))
    Nt = (math.floor(n1t + 0.5), math.floor(th * n1t + 0.5))
    got = (float(r['N1']), float(r['N2']), float(r['N1.tab']), float(r['N2.tab']))
    d_n = max(d_n, max(abs(a - b) for a, b in zip(got, N + Nt)))
    Np = (int(float(r['N1.pub'])), int(float(r['N2.pub'])))
    pw = power(Np[0], Np[1], p1, p2, R0, 0.05, False)
    d_p = max(d_p, abs(pw - float(r['power'])))
    pub_n += N == Np
    pub_tab += Nt == Np
    pub_pw += abs(math.floor(pw * 1e4 + 1e-9) / 1e4 - float(r['pw.pub'])) < 1e-9
print('Table I: largest difference in the sample sizes {:g}, in the powers {:.2e}'.format(
    d_n, d_p))
print('Table I: published sample sizes given by the standard normal quantiles {}, by '
      '1.645 and 1.282 {}, published powers after truncation {}, of {}'.format(
          pub_n, pub_tab, pub_pw, len(t1)))
ok = ok and d_n == 0 and d_p < TOL
d_n2 = 0
pub_n2 = 0
for r in t2:
    p1, p2, R0 = (float(r[k]) for k in ('p1', 'p2', 'R0'))
    n = (za + zb) ** 2 * (p1 * (1 - p1) + R0 ** 2 * p2 * (1 - p2)) / (p1 - R0 * p2) ** 2
    N = math.floor(n + 0.5)
    d_n2 = max(d_n2, abs(float(r['N1']) - N), abs(float(r['N2']) - N))
    pub_n2 += N == int(float(r['N.pub']))
print('Table II, Method 1: largest difference in the sample sizes {:g}, published sample '
      'sizes given {} of {}'.format(d_n2, pub_n2, len(t2)))
ok = ok and d_n2 == 0
n = formula8(0.01, 0.01, 1, 1.5, Z.inv_cdf(0.975), zb)
N = math.ceil(n - 1e-9)
pw = power(N, N, 0.01, 0.01, 1.5, 0.025, True, 450, 450)
print('Second example: sample size {} per group (stored {:g} and {:g}), power {:.6f} '
      '(stored {:.6f})'.format(N, float(ex['N1']), float(ex['N2']), pw, float(ex['power'])))
ok = (ok and float(ex['N1']) == N and float(ex['N2']) == N and
      abs(pw - float(ex['power'])) < TOL)
print('PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
