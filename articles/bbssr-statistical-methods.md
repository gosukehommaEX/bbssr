# Statistical Methods in bbssr

``` r

library(bbssr)
```

## Notation

Let $`X_{j}`$ denote the number of responders in group $`j`$, so that
$`X_{1} \sim \mathrm{Bin}(N_{1}, p_{1})`$ and
$`X_{2} \sim \mathrm{Bin}(N_{2}, p_{2})`$ independently. Realized counts
are written $`x_{1}`$ and $`x_{2}`$, and $`s = x_{1} + x_{2}`$ is the
total number of responders. The null hypothesis is
$`H_{0}: p_{1} = p_{2}`$, and the common value under the null is denoted
$`\theta`$ and treated as a nuisance parameter.

Every test in the package is defined through its rejection region, a
subset of the $`(N_{1} + 1) \times (N_{2} + 1)`$ grid of possible
outcomes. Once that region $`\mathcal{R}`$ is available, the power at
any pair of response probabilities follows from a single sum,

``` math
1 - \beta = \sum_{(x_{1}, x_{2}) \in \mathcal{R}}
\binom{N_{1}}{x_{1}} p_{1}^{x_{1}} (1 - p_{1})^{N_{1} - x_{1}}
\binom{N_{2}}{x_{2}} p_{2}^{x_{2}} (1 - p_{2})^{N_{2} - x_{2}} .
```

The type I error rate is the same sum evaluated at
$`p_{1} = p_{2} = \theta`$, maximized over $`\theta`$. The tests of
non-inferiority with a margin, `Blackwelder` and `Farrington-Manning`,
have the null hypothesis $`p_{1} - p_{2} \le -\delta`$ instead, and the
`bbssr-non-inferiority` vignette describes them. The rest of this
vignette concerns the hypothesis of equal response probabilities.

## Conditional tests

Conditioning on $`s`$ removes the nuisance parameter. Under the null the
count $`X_{1}`$ then follows a hypergeometric distribution, and the
one-sided Fisher p-value is

``` math
p_{F}(x_{1}, x_{2}) = \Pr(X_{1} \ge x_{1} \mid s)
= \sum_{k \ge x_{1}} \frac{\binom{N_{1}}{k} \binom{N_{2}}{s - k}}{\binom{N_{1} + N_{2}}{s}} .
```

The mid-p variant replaces the contribution of the observed table by
half of it, giving
$`\Pr(X_{1} > x_{1} \mid s) + \tfrac{1}{2} \Pr(X_{1} = x_{1} \mid s)`$.
This is no longer a valid p-value in the strict sense, so the mid-p test
can exceed the nominal level, but it removes much of the conservatism
that conditioning introduces.

### Two-sided conventions

A two-sided version of a discrete conditional test is not unique. Three
conventions are available through the `tsmethod` argument.

The `minlike` convention sums the null probabilities of all tables that
are no more likely than the observed one,
``` math
p(x_{1}, x_{2}) = \sum_{k \,:\, f(k) \le f(x_{1})} f(k), \qquad f(k) = \Pr(X_{1} = k \mid s) .
```
This is the convention of
[`stats::fisher.test`](https://rdrr.io/r/stats/fisher.test.html).

The `central` convention doubles the smaller of the two one-sided tail
probabilities and truncates at one,
``` math
p(x_{1}, x_{2}) = \min\bigl\{1, \; 2 \min(\Pr(X_{1} \le x_{1} \mid s), \Pr(X_{1} \ge x_{1} \mid s))\bigr\} .
```

The `blaker` convention orders the tables by the smaller of their two
one-sided tail probabilities and sums the null probabilities of all
tables that are at least as extreme as the observed one in this
ordering,
``` math
p(x_{1}, x_{2}) = \sum_{k \,:\, g(k) \le g(x_{1})} f(k), \qquad
g(k) = \min\bigl\{\Pr(X_{1} \le k \mid s), \Pr(X_{1} \ge k \mid s)\bigr\} .
```
This is formula (2) of Mehrotra, Chan and Berger (2003), and the
convention `'blaker'` of
[`exact2x2::exact2x2()`](https://rdrr.io/pkg/exact2x2/man/exact2x2.html).
The tables of the other tail that it adds form a tail of their own,
whose probability is no larger than the tail probability of the observed
table. The p-value is therefore at most twice that tail probability and
never exceeds the central one. When the two groups are of equal size the
conditional distribution is symmetric, and the three conventions give
the same p-value.

Fay and Hunsberger (2021, Section 8) illustrate the differences with 8
responders out of 14 against 1 out of 7.

``` r

ex <- vapply(c('minlike', 'central', 'blaker'), function(ts) {
  attr(BinaryRR(14, 7, 0.05, 'Fisher', alternative = 'two.sided', tsmethod = ts),
       'p.value')[9, 2]
}, numeric(1))
round(ex, 3)
#> minlike central  blaker 
#>   0.159   0.157   0.087
```

The table with 4 responders out of 14 against 5 out of 7 is exactly as
likely as the observed one, so the `minlike` convention counts it. Its
lower tail probability, 0.0805, exceeds the upper tail probability of
the observed table, 0.0783, so the `blaker` convention leaves it out.

Under the mid-p correction, the `minlike` and `blaker` conventions count
the tables tied with the observed one in their ordering, the observed
one included, with half of their probability, which is the definition of
the mid-p value in Fay and Hunsberger (2021, Section 9). The `central`
convention doubles the smaller of the two one-sided mid-p values.

The bound by the central p-value holds for the exact Fisher p-value
only. It fails for the mid-p value, and for the Boschloo test, whose
ordering of the outcomes changes with the convention, as the two
outcomes below show: 5 responders out of 8 against none out of 5 for the
mid-p test, and 2 out of 4 against none out of 3 for the Boschloo test.

``` r

two_sided_p <- function(N1, N2, Test, x1, x2, ts) {
  attr(BinaryRR(N1, N2, 0.05, Test, alternative = 'two.sided', tsmethod = ts),
       'p.value')[x1 + 1, x2 + 1]
}
data.frame(
  Test = c('Fisher-midP', 'Boschloo'),
  blaker = c(two_sided_p(8, 5, 'Fisher-midP', 5, 0, 'blaker'),
             two_sided_p(4, 3, 'Boschloo', 2, 0, 'blaker')),
  central = c(two_sided_p(8, 5, 'Fisher-midP', 5, 0, 'central'),
              two_sided_p(4, 3, 'Boschloo', 2, 0, 'central'))
)
#>          Test     blaker    central
#> 1 Fisher-midP 0.05361305 0.04351204
#> 2    Boschloo 0.30541895 0.21874043
```

The `central` convention has a property that the `minlike` and `blaker`
conventions lack. Its two-sided rejection region at level $`2\alpha`$ is
exactly the union of the two one-sided rejection regions at level
$`\alpha`$.

``` r

N1 <- 9
N2 <- 7
alpha <- 0.02
two <- BinaryRR(N1, N2, 2 * alpha, 'Fisher',
                alternative = 'two.sided', tsmethod = 'central')
upper <- BinaryRR(N1, N2, alpha, 'Fisher')
lower <- t(BinaryRR(N2, N1, alpha, 'Fisher'))
identical(as.vector(two), as.vector(upper | lower))
#> [1] TRUE
```

The `minlike` and `central` conventions give different regions of the
same nominal size.

``` r

minlike <- BinaryRR(N1, N2, 0.05, 'Fisher', alternative = 'two.sided',
                    tsmethod = 'minlike')
central <- BinaryRR(N1, N2, 0.05, 'Fisher', alternative = 'two.sided',
                    tsmethod = 'central')
data.frame(
  tsmethod = c('minlike', 'central'),
  rejected = c(sum(minlike), sum(central)),
  rejected.by.this.only = c(sum(minlike & !central), sum(central & !minlike))
)
#>   tsmethod rejected rejected.by.this.only
#> 1  minlike       22                     2
#> 2  central       20                     0
```

The chi-squared and Z-pooled tests order outcomes by $`|Z|`$ when the
alternative is two-sided, so `tsmethod` does not apply to them.

## Unconditional tests

Conditioning is not the only way to eliminate $`\theta`$. An exact
unconditional test keeps the full binomial model and maximizes the null
tail probability over the nuisance parameter,

``` math
p(x_{1}, x_{2}) = \sup_{0 \le \theta \le 1}
\Pr_{\theta}\bigl(T(X_{1}, X_{2}) \text{ at least as extreme as } T(x_{1}, x_{2})\bigr) ,
```

where $`T`$ is an ordering statistic. The Z-pooled test uses the
two-sample Z statistic with a pooled variance estimator,

``` math
Z(x_{1}, x_{2}) = \frac{x_{1} / N_{1} - x_{2} / N_{2}}
{\sqrt{\hat{p}(1 - \hat{p})(1 / N_{1} + 1 / N_{2})}}, \qquad
\hat{p} = \frac{x_{1} + x_{2}}{N_{1} + N_{2}} ,
```

with larger values more extreme. The Boschloo test uses the Fisher
p-value itself as the ordering statistic, with smaller values more
extreme. Boschloo (1970) described this test as Fisher’s test at a
raised conditional level, the largest level at which the probability of
rejection under the null hypothesis does not exceed $`\alpha`$ at any
$`\theta`$. The probability $`\Pr_{\theta}(p_{F} \le c)`$ does not
decrease as $`c`$ increases, so rejecting when $`p_{F}`$ is at most the
raised level is the same as rejecting when the unconditional p-value
above is at most $`\alpha`$. For a two-sided alternative the article
gives each of the two parts of the rejection region half of the raised
level, which is the `central` convention. The default
`tsmethod = 'minlike'` orders the outcomes by the two-sided p-value of
[`stats::fisher.test`](https://rdrr.io/r/stats/fisher.test.html)
instead, and `tsmethod = 'blaker'` by the p-value of formula (2) of
Mehrotra, Chan and Berger (2003), which is how they define the two-sided
Boschloo test.

The supremum is approximated by a search over the grid of `n.grid`
equally spaced values of $`\theta`$ from 0 to 1, where `n.grid` defaults
to 100. A grid that contains a coarser grid can only find a maximum at
least as large, so refining the grid in this way can only raise the
p-values. The grids of 21 and 2001 points below are nested, since the
spacing of the first is a multiple of that of the second.

``` r

coarse <- attr(BinaryRR(12, 12, 0.025, 'Boschloo', n.grid = 21), 'p.value')
fine <- attr(BinaryRR(12, 12, 0.025, 'Boschloo', n.grid = 2001), 'p.value')
data.frame(all.p.values.increased = all(fine >= coarse - 1e-12),
           largest.increase = max(fine - coarse))
#>   all.p.values.increased largest.increase
#> 1                   TRUE       0.00222195
```

### Refining the maximum between the grid points

The tail probability of each outcome is a polynomial in $`\theta`$, and
its maximum over the grid is a lower bound of the maximum over the unit
interval. A p-value computed on the grid can therefore fall below the
exact p-value, and the test can then exceed its nominal level. With
`ref.pvalue = TRUE` the grid is extended by points equally spaced on the
arcsine square-root scale, which resolve the narrow local maxima near 0
and 1, and every local maximum on the extended grid is refined by a
safeguarded Newton iteration between its two neighbouring grid points. A
refined p-value is never smaller than the p-value on the grid.

The size of a fixed-sample design is the largest rejection probability
over $`\theta`$, which the function below evaluates on a fine grid and
refines with [`optimize()`](https://rdrr.io/r/stats/optimize.html)
around the largest grid value.

``` r

size_fixed <- function(RR) {
  N1 <- attr(RR, 'N1')
  N2 <- attr(RR, 'N2')
  m <- matrix(as.vector(RR), N1 + 1L, N2 + 1L)
  f <- function(t) sum(outer(dbinom(0:N1, N1, t), dbinom(0:N2, N2, t)) * m)
  theta <- seq(0.001, 0.999, by = 0.001)
  v <- vapply(theta, f, numeric(1))
  k <- which.max(v)
  lower <- theta[max(k - 1L, 1L)]
  upper <- theta[min(k + 1L, length(theta))]
  max(v[k], optimize(f, c(lower, upper), maximum = TRUE)$objective)
}
grid.rr <- BinaryRR(32, 32, 0.025, 'Z-pool')
refined.rr <- BinaryRR(32, 32, 0.025, 'Z-pool', ref.pvalue = TRUE)
data.frame(
  p.values = c('grid', 'refined'),
  rejected = c(sum(grid.rr), sum(refined.rr)),
  size = signif(c(size_fixed(grid.rr), size_fixed(refined.rr)), 4)
)
#>   p.values rejected    size
#> 1     grid      342 0.02501
#> 2  refined      340 0.02334
```

With 32 patients per group the refinement removes 2 outcomes from the
rejection region of the Z-pooled test at the one-sided level 0.025, and
the size of the test changes from 0.02501 to 0.02334. Every function
that computes a rejection region accepts `ref.pvalue`. The default is
`FALSE`, which keeps the results of earlier versions, and the refinement
takes longer than the grid alone.

### Ties in the ordering statistic

The ordering statistic takes the same value at several outcomes far more
often than one might expect. With $`N_{1} = N_{2} = 7`$ the outcomes
$`(x_{1}, x_{2}) = (5, 1)`$ and $`(6, 2)`$ both have a Fisher p-value of
$`2/39`$.

``` r

stat <- attr(BinaryRR(7, 7, 0.025, 'Fisher'), 'p.value')
c(cell_5_1 = stat[6, 2], cell_6_2 = stat[7, 3])
#>   cell_5_1   cell_6_2 
#> 0.05128205 0.05128205
```

The tail event is defined by “at least as extreme as”, so both outcomes
belong to each other’s tail set and must receive the same p-value.
Accumulating the null probabilities in an arbitrary order within a tie
group would give them different values, and the decision at those
outcomes would depend on how the sorting routine happens to break the
tie. The package groups tied values explicitly and assigns each group
the tail probability accumulated up to its last member.

``` r

p <- attr(BinaryRR(7, 7, 0.025, 'Boschloo'), 'p.value')
c(cell_5_1 = p[6, 2], cell_6_2 = p[7, 3])
#>   cell_5_1   cell_6_2 
#> 0.02867981 0.02867981
```

### The Berger-Boos procedure

Maximizing over the whole unit interval is wasteful, because values of
$`\theta`$ far from the observed pooled proportion are implausible.
Berger and Boos (1994) proposed maximizing over a $`100(1 - \gamma)`$
percent confidence interval $`C_{\gamma}`$ for $`\theta`$ instead, and
adding $`\gamma`$ to compensate,

``` math
p_{\gamma}(x_{1}, x_{2}) = \gamma + \sup_{\theta \in C_{\gamma}}
\Pr_{\theta}\bigl(T(X_{1}, X_{2}) \text{ at least as extreme as } T(x_{1}, x_{2})\bigr) .
```

The package uses an exact Clopper-Pearson interval based on $`s`$
responders among $`N_{1} + N_{2}`$ patients, so the interval differs
from outcome to outcome. Setting `bb.gamma` to a positive value
activates the procedure. Berger and Boos (1994) write $`\beta`$ for
$`\gamma`$ and suggest a small value such as 0.001 or 0.0001. The
resulting test still controls the type I error rate at the nominal
level.

``` r

plain <- BinaryRR(15, 15, 0.025, 'Boschloo', n.grid = 200)
bb <- BinaryRR(15, 15, 0.025, 'Boschloo', n.grid = 200, bb.gamma = 1e-4)
p.plain <- attr(plain, 'p.value')
p.bb <- attr(bb, 'p.value')
data.frame(rejected.plain = sum(plain), rejected.berger.boos = sum(bb),
           largest.decrease = max(p.plain - p.bb),
           largest.increase = max(p.bb - p.plain))
#>   rejected.plain rejected.berger.boos largest.decrease largest.increase
#> 1             61                   61      0.005134331     0.0001087237
```

Two forces act in opposite directions. Restricting the search lowers the
maximum, and the additive $`\gamma`$ raises the result. Which one wins
varies from outcome to outcome, so the rejection region can grow, shrink
or stay as it is. The gain is largest when the outcome is extreme,
because the confidence interval for $`\theta`$ then excludes the value
at which the unrestricted maximum is attained.

## Relationships between the tests

The conditional distribution of a p-value satisfies
$`\Pr(p_{F} \le c \mid s) \le c`$ for any fixed $`c`$. Averaging over
$`s`$ gives $`\Pr_{\theta}(p_{F} \le c) \le c`$ for every $`\theta`$, so
the Boschloo p-value never exceeds the Fisher p-value at the same
outcome. The Boschloo rejection region therefore contains the Fisher
rejection region, which is the sense in which Boschloo’s test is
uniformly more powerful.

``` r

fisher <- BinaryRR(15, 15, 0.025, 'Fisher')
boschloo <- BinaryRR(15, 15, 0.025, 'Boschloo', n.grid = 200)
data.frame(rejected.fisher = sum(fisher), rejected.boschloo = sum(boschloo),
           fisher.region.contained = all(as.vector(boschloo)[as.vector(fisher)]))
#>   rejected.fisher rejected.boschloo fisher.region.contained
#> 1              51                61                    TRUE
```

The type I error rates show how much of the nominal level each test
actually spends.

``` r

max_type1 <- function(RR, n.grid = 401) {
  N1 <- attr(RR, 'N1')
  N2 <- attr(RR, 'N2')
  m <- matrix(as.vector(RR), N1 + 1L, N2 + 1L)
  theta <- seq(0, 1, length.out = n.grid)
  max(vapply(theta, function(t) {
    sum(outer(dbinom(0:N1, N1, t), dbinom(0:N2, N2, t)) * m)
  }, numeric(1)))
}
tests <- c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')
t1 <- vapply(tests, function(tst) {
  max_type1(BinaryRR(30, 30, 0.025, tst, n.grid = 200), n.grid = 801)
}, numeric(1))
data.frame(
  Test = tests, max.type1 = round(t1, 5), exceeds.alpha = t1 > 0.025,
  row.names = NULL
)
#>          Test max.type1 exceeds.alpha
#> 1       Chisq   0.02770          TRUE
#> 2      Fisher   0.01370         FALSE
#> 3 Fisher-midP   0.02595          TRUE
#> 4      Z-pool   0.02346         FALSE
#> 5    Boschloo   0.02344         FALSE
```

The Fisher test spends the least, which is the cost of conditioning. The
two unconditional tests spend much more while staying below the level,
which is where their extra power comes from. The chi-squared and mid-p
tests carry no such guarantee, and the `exceeds.alpha` column shows what
that means at this configuration.

## Blinded sample size re-estimation

At the interim analysis, $`n_{1}`$ and $`n_{2}`$ patients have been
observed and the total number of responders $`S`$ is known. The blinded
estimate of the pooled response probability is
$`\hat{p} = S / (n_{1} + n_{2})`$. With an allocation ratio of $`r`$ to
1 and an assumed treatment effect $`\Delta_{A}`$, group-specific
probabilities are recovered as

``` math
\hat{p}_{1} = \hat{p} + \frac{\Delta_{A}}{1 + r}, \qquad
\hat{p}_{2} = \hat{p} - \frac{r \Delta_{A}}{1 + r} ,
```

and these enter the sample size calculation in place of the original
assumptions. Nothing in this chain requires knowledge of which patient
received which treatment. With `effect = 'RR'` or `effect = 'OR'` the
assumed effect is a risk ratio or an odds ratio, and the recovered
probabilities are those with the pooled value $`\hat{p}`$ and that
ratio.

The sample size is re-estimated in one of two ways, chosen by
`ss.method`. The default `'exact'` steps from the normal approximation,
one patient at a time, to the sample size at which the exact power of
the test first attains the target, as
[`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)
does, and this search needs probabilities inside the unit interval, so
the recovered probabilities are truncated to it. The other three methods
use the normal approximation, which gives the size of group 2 as
``` math
n_{2} = \frac{\bigl(z_{1 - \alpha} \sqrt{v_{0}} + z_{1 - \beta} \sqrt{v_{1}}\bigr)^{2}}{(\hat{p}_{1} - \hat{p}_{2})^{2}}, \qquad
v_{0} = \hat{p}(1 - \hat{p}) \Bigl(1 + \frac{1}{r}\Bigr) ,
```
where
$`v_{1} = \hat{p}_{1}(1 - \hat{p}_{1}) / r + \hat{p}_{2}(1 - \hat{p}_{2})`$
under `'standard'`, formula (21.3) of Kieser (2020), and
$`v_{1} = v_{0}`$ under `'null.variance'`, formula (1) of Friede and
Kieser (2004). Under `'alternative.variance'` both terms use $`v_{1}`$.
A two-sided alternative uses $`\alpha / 2`$ in place of $`\alpha`$, and
for the risk difference the denominator is $`\Delta_{A}^{2}`$. Group 1
receives $`r n_{2}`$ patients. With a non-inferiority margin the
difference and the variance $`v_{0}`$ change, as the
`bbssr-non-inferiority` vignette describes. The approximation uses the
recovered probabilities before truncation, with each Bernoulli variance
truncated at zero, so that the assumed effect $`\Delta_{A}`$ is kept as
in formula (2) of Friede and Kieser (2004). The
`bbssr-reestimation-rules` vignette compares the methods and the rules
that turn $`n_{2}`$ into whole numbers.

Two rules govern what happens next. The unrestricted rule takes the
re-estimated sample size as it stands, allowing the trial to end up
smaller than planned. The restricted rule raises it to the planned
sample size first, so the trial can only grow. The table below applies
both to a trial planned at 24 patients per group with an interim
analysis at 12, with the default exact search.

``` r

interim <- data.frame(S = c(2, 4, 6, 8, 10))
interim$pooled <- round(interim$S / 24, 3)
interim$unrestricted <- vapply(interim$S, function(s) {
  BinaryBSSR(n1 = 12, n2 = 12, S = s, Delta.A = 0.36, r = 1,
             alpha = 0.025, tar.power = 0.8, Test = 'Chisq')$N.final
}, numeric(1))
interim$restricted <- vapply(interim$S, function(s) {
  BinaryBSSR(n1 = 12, n2 = 12, S = s, Delta.A = 0.36, r = 1,
             alpha = 0.025, tar.power = 0.8, Test = 'Chisq',
             restricted = TRUE, N1 = 24, N2 = 24)$N.final
}, numeric(1))
interim
#>    S pooled unrestricted restricted
#> 1  2  0.083           40         48
#> 2  4  0.167           30         48
#> 3  6  0.250           40         48
#> 4  8  0.333           50         50
#> 5 10  0.417           58         58
```

The unrestricted column is not monotone in $`S`$. At a pooled rate of
0.083 the recovered control probability $`\hat{p} - \Delta_{A} / 2`$ is
negative and is truncated at zero, which shrinks the recovered risk
difference below $`\Delta_{A}`$ and pushes the sample size back up.

[`BinaryPowerBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryPowerBSSR.md)
evaluates a design by averaging the conditional power over the
distribution of the interim outcome,

``` math
1 - \beta_{\mathrm{BSSR}} = \sum_{x_{1}, x_{2}}
\Pr(x_{1} \mid n_{1}, p_{1}) \Pr(x_{2} \mid n_{2}, p_{2}) \,
\mathrm{CP}(x_{1}, x_{2}) ,
```

where the conditional power $`\mathrm{CP}`$ is computed from the
rejection region of the final sample size that the interim outcome leads
to. The sum runs over every possible interim outcome, so the same
rejection region is required many times and is cached.

Setting $`p_{1} = p_{2} = \theta`$ in this sum gives the type I error
rate of the design at $`\theta`$. The final sample size depends on the
interim data, so the rate is not bounded by the size of the test at any
one sample size, and
[`BinaryTypeIErrorBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryTypeIErrorBSSR.md)
evaluates it over $`\theta`$. The `bbssr-type1-error` vignette covers
this calculation, the adjusted nominal level and the decomposition of
the rate by the interim outcome.

## Choosing a test

`Boschloo` is the default recommendation when the computation is
affordable, since it is exact and uniformly more powerful than `Fisher`.
`Z-pool` is close behind and slightly cheaper. Mehrotra, Chan and Berger
(2003) recommend these two tests over Fisher’s test, and caution that an
exact unconditional test based on another statistic, such as the
difference in proportions, can be less powerful than Fisher’s test.
`Fisher` is exact but conservative, and is the conventional choice when
a regulator expects the classical procedure. `Chisq` is useful for
exploration and for the starting value of a sample size search, but does
not control the type I error rate exactly at small sample sizes.
`Fisher-midP` sits between `Fisher` and the unconditional tests and is
worth considering when exact control is not a formal requirement.
`Blackwelder` and `Farrington-Manning` are the tests to use with a
non-inferiority margin. Their rejection regions and power are computed
exactly, but the tests themselves are asymptotic, so their size at a
given sample size is worth checking.

## References

Berger, R. L. and Boos, D. D. (1994). P values maximized over a
confidence set for the nuisance parameter. *Journal of the American
Statistical Association*, 89, 1012-1016.

Boschloo, R. D. (1970). Raised conditional level of significance for the
2 × 2-table when testing the equality of two probabilities. *Statistica
Neerlandica*, 24, 1-9.

Fay, M. P. and Hunsberger, S. A. (2021). Practical valid inferences for
the two-sample binomial problem. *Statistics Surveys*, 15, 72-110.

Friede, T. and Kieser, M. (2004). Sample size recalculation for binary
data in internal pilot study designs. *Pharmaceutical Statistics*, 3,
269-279.

Kieser, M. (2020). *Methods and Applications of Sample Size Calculation
and Recalculation in Clinical Trials*. Springer.

Mehrotra, D. V., Chan, I. S. F. and Berger, R. L. (2003). A cautionary
note on exact unconditional inference for a difference between two
independent binomial proportions. *Biometrics*, 59, 441-450.
