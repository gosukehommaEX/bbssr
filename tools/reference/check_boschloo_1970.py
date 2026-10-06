# Independent check of the Boschloo (1970) items of inst/reproduce/reproduce-published.R.
# The numbers of Sections 2, 3, 4 and 6 of Boschloo (1970), Statistica Neerlandica 24, 1-9,
# are computed without the code of bbssr: exact rational conditional p-values (fractions),
# binomial sums with numpy and the maximum over the common probability with scipy. The
# predicted verdicts follow the rules of the reproduction script, and the recomputed
# values of bbssr in reproduce-output/published-comparison.csv are compared with the
# predictions side by side. Run from the package root where scipy is available:
#   python3 tools/reference/check_boschloo_1970.py [path to published-comparison.csv]
# With the argument 'perturb' in place of the path, two published values and one count
# are shifted, and the predicted verdicts must then include FAIL.
import csv
import sys
from fractions import Fraction
from math import comb
import numpy as np
from scipy.optimize import minimize_scalar


def fisher_upper(m, n):
    """P(A >= a | R = a + b) under H0, as Fractions, indexed [a][b]."""
    p = [[None] * (n + 1) for _ in range(m + 1)]
    for a in range(m + 1):
        for b in range(n + 1):
            r = a + b
            den = comb(m + n, r)
            num = sum(comb(m, k) * comb(n, r - k) for k in range(a, min(m, r) + 1))
            p[a][b] = Fraction(num, den)
    return p


def fisher_lower(m, n):
    p = [[None] * (n + 1) for _ in range(m + 1)]
    for a in range(m + 1):
        for b in range(n + 1):
            r = a + b
            den = comb(m + n, r)
            num = sum(comb(m, k) * comb(n, r - k) for k in range(max(0, r - n), a + 1))
            p[a][b] = Fraction(num, den)
    return p


def central(m, n):
    up, lo = fisher_upper(m, n), fisher_lower(m, n)
    return [[min(Fraction(1), 2 * min(up[a][b], lo[a][b])) for b in range(n + 1)]
            for a in range(m + 1)]


def binom_pmf(N, t):
    k = np.arange(N + 1)
    c = np.array([comb(N, int(x)) for x in k], dtype=float)
    return c * t ** k * (1 - t) ** (N - k)


def prob_region(region, m, n, p1, p2):
    return float(binom_pmf(m, p1) @ region.astype(float) @ binom_pmf(n, p2))


def size(region, m, n, n_grid=4001):
    f = lambda t: prob_region(region, m, n, t, t)
    th = np.linspace(0, 1, n_grid)
    v = np.array([f(t) for t in th])
    best = v.max()
    for k in np.argsort(v)[-6:]:
        lo, hi = th[max(k - 1, 0)], th[min(k + 1, n_grid - 1)]
        if hi > lo:
            res = minimize_scalar(lambda t: -f(t), bounds=(lo, hi), method='bounded',
                                  options={'xatol': 1e-12})
            best = max(best, -res.fun)
    return best


def region_of(P, c):
    return np.array([[P[a][b] <= c for b in range(len(P[0]))] for a in range(len(P))])


def raised_level(P, m, n, alpha):
    """Largest attainable conditional level gamma with unconditional size <= alpha,
    and the next attainable value above it."""
    vals = sorted({P[a][b] for a in range(m + 1) for b in range(n + 1)})
    best = None
    for i, c in enumerate(vals):
        if size(region_of(P, c), m, n) <= alpha:
            best = i
        else:
            break
    nxt = vals[best + 1] if best + 1 < len(vals) else None
    return vals[best], nxt


def grid_region(P, m, n, alpha, n_grid=100):
    """Boschloo region with the maximum over seq(0, 1, length.out = 100), as in bbssr."""
    th = np.linspace(0, 1, n_grid)
    f1 = np.array([binom_pmf(m, t) for t in th])
    f2 = np.array([binom_pmf(n, t) for t in th])
    reg = np.zeros((m + 1, n + 1), bool)
    for a in range(m + 1):
        for b in range(n + 1):
            R = region_of(P, P[a][b]).astype(float)
            reg[a, b] = np.einsum('ga,ab,gb->g', f1, R, f2).max() <= alpha
    return reg


items = []  # (item, published, recomputed, digits)
rows = [(0.01, .3, .1), (0.01, .6, .1), (0.01, .7, .2), (0.01, .8, .2),
        (0.05, .3, .1), (0.05, .6, .1), (0.05, .7, .2), (0.05, .8, .2)]
fpub = [0.0755, 0.6087, 0.5268, 0.7647, 0.2558, 0.8451, 0.8066, 0.9391]
bpub = [0.1493, 0.7394, 0.6796, 0.8723, 0.3531, 0.9189, 0.8962, 0.9744]
up15 = fisher_upper(15, 15)
reg = {}
for alpha in (0.01, 0.05):
    reg[alpha] = (region_of(up15, Fraction(alpha).limit_denominator(1000)),
                  grid_region(up15, 15, 15, alpha))
for k, (alpha, p1, p2) in enumerate(rows):
    lab = f'Section 4: alpha = {alpha:g}, p1 = {p1:g}, p2 = {p2:g}'
    KF, KB = reg[alpha]
    items.append((lab + ': power, Fisher test (column I)', fpub[k],
                  prob_region(KF, 15, 15, p1, p2), 4))
    items.append((lab + ': power, raised level (column II)', bpub[k],
                  prob_region(KB, 15, 15, p1, p2), 4))
up = fisher_upper(15, 10)
cen = central(15, 10)
KF = region_of(up, Fraction(5, 100))
KB = grid_region(up, 15, 10, 0.05)
KB2 = grid_region(cen, 15, 10, 0.05)
th = np.arange(0, 1 + 1e-12, 1e-4)
sz = max(prob_region(KF, 15, 10, t, t) for t in th)
pf = np.array([[float(up[a][b]) for b in range(11)] for a in range(16)])
pf2 = np.array([[float(cen[a][b]) for b in range(11)] for a in range(16)])
items += [
    ('Section 6: Fisher p-value', 0.0565, pf[5, 0], 4),
    ('Section 6: rejected by the Fisher test', 0, float(KF[5, 0]), 0),
    ('Section 6: rejected at the raised level, one-sided', 1, float(KB[5, 0]), 0),
    ('Section 6: rejected at the raised level, two-sided', 1, float(KB2[5, 0]), 0),
    ('Section 2: size of the Fisher test', 0.02, sz, 2),
    ('Section 6: outcomes added', 8, float((KB & ~KF).sum()), 0),
    ('Section 6: outcomes removed', 0, float((KF & ~KB).sum()), 0),
    ('Sections 3 and 6: differ at 0.09', 0, float(((pf <= 0.09) != KB).sum()), 0),
    ('Section 6: differ at 0.114', 0, float(((pf2 <= 0.114) != KB2).sum()), 0),
]
explained = {'Section 4: alpha = 0.05, p1 = 0.6, p2 = 0.1: power, Fisher test (column I)':
             6e-5}
if len(sys.argv) > 1 and sys.argv[1] == 'perturb':  # self-test
    # Shift two published values by one unit in the last digit and one count by one
    items[0] = (items[0][0], items[0][1] + 1e-4, items[0][2], items[0][3])
    items[10] = (items[10][0], items[10][1] - 1e-4, items[10][2], items[10][3])
    items[21] = (items[21][0], 9, items[21][2], items[21][3])
counts = {}
results = []
for name, pub, rec, d in items:
    shown = round(rec, d)
    v = 'PASS' if abs(shown - pub) < 1e-9 else 'FAIL'
    if v == 'FAIL' and name in explained and abs(rec - pub) <= explained[name] + 1e-12:
        v = 'EXPLAINED'
    counts[v] = counts.get(v, 0) + 1
    results.append((v, name, pub, rec))
    print(f'{v:9s} | {name} | published {pub:g} | recomputed {rec:.17g}')
print('predicted items:', len(items), counts)


def compare(path):
    """Side-by-side comparison of bbssr (CSV) and the predictions; returns mismatches."""
    rows = [r for r in csv.DictReader(open(path, encoding='utf-8'))
            if r['source'] == 'Boschloo (1970)']
    pred = {name: (v, pub, rec) for v, name, pub, rec in results}
    def key(item):
        k = item.replace('Boschloo1970 ', '').replace(' (1 = yes)', '')
        k = k.replace('5 / 15 against 0 / 10: ', '')
        k = k.replace('rejected by the Fisher test at 0.05', 'rejected by the Fisher test')
        short_names = [
            ('Section 2: size of the Fisher test', 'Section 2: size of the Fisher test'),
            ('Section 6: outcomes added', 'Section 6: outcomes added'),
            ('Section 6: outcomes removed', 'Section 6: outcomes removed'),
            ('Sections 3 and 6: outcomes on which', 'Sections 3 and 6: differ at 0.09'),
            ('Section 6: outcomes on which', 'Section 6: differ at 0.114')]
        for start, short in short_names:
            if k.startswith(start):
                k = short
        return k
    bad = 0
    print(f'{"item":70s} {"published":>9s} {"bbssr":>12s} {"python":>12s} '
          f'{"verdict":>9s} {"predicted":>9s}')
    for r in rows:
        v, pub, rec = pred[key(r['item'])]
        rr = float(r['recomputed'])
        ok = (abs(rr - rec) <= 1e-9 * max(1.0, abs(rec)) and r['verdict'] == v
              and abs(float(r['published']) - pub) < 1e-12)
        bad += not ok
        print(f'{r["item"][13:83]:70s} {float(r["published"]):9.4g} {rr:12.8g} {rec:12.8g} '
              f'{r["verdict"]:>9s} {v:>9s} {"" if ok else "MISMATCH"}')
    print(f'Boschloo items in the CSV: {len(rows)}, predictions: {len(pred)}, mismatches: {bad}')
    return bad

if len(sys.argv) > 1 and sys.argv[1] != 'perturb':
    sys.exit(1 if compare(sys.argv[1]) else 0)
