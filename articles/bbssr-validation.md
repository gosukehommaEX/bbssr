# Validation of bbssr

``` r

library(bbssr)
has_exact <- requireNamespace('Exact', quietly = TRUE)
has_exact2x2 <- requireNamespace('exact2x2', quietly = TRUE)
has_bench <- requireNamespace('microbenchmark', quietly = TRUE)
c(Exact = has_exact, exact2x2 = has_exact2x2, microbenchmark = has_bench)
#>          Exact       exact2x2 microbenchmark 
#>           TRUE           TRUE           TRUE
```

Chunks that need one of the optional packages are skipped when it is
unavailable.

## Conditional tests against stats

The Fisher p-values are compared with
[`stats::fisher.test`](https://rdrr.io/r/stats/fisher.test.html) at
every outcome of a small grid, for the one-sided alternative and for the
two-sided `minlike` convention.

``` r

N1 <- 7
N2 <- 6
compare_fisher <- function(alternative, tsmethod) {
  got <- attr(BinaryRR(N1, N2, 0.05, 'Fisher', alternative = alternative,
                       tsmethod = tsmethod), 'p.value')
  want <- outer(0:N1, 0:N2, Vectorize(function(i, j) {
    tab <- matrix(c(i, j, N1 - i, N2 - j), nrow = 2)
    if (alternative == 'greater') {
      stats::fisher.test(tab, alternative = 'greater')$p.value
    } else {
      stats::fisher.test(tab)$p.value
    }
  }))
  max(abs(got - want))
}
data.frame(
  comparison = c('one-sided', 'two-sided, minlike'),
  max.absolute.difference = c(
    compare_fisher('greater', 'minlike'),
    compare_fisher('two.sided', 'minlike')
  )
)
#>           comparison max.absolute.difference
#> 1          one-sided            4.440892e-16
#> 2 two-sided, minlike            5.551115e-16
```

The `central` convention has no counterpart in `stats`, so it is checked
against its definition.

``` r

got <- attr(BinaryRR(N1, N2, 0.05, 'Fisher', alternative = 'two.sided',
                     tsmethod = 'central'), 'p.value')
want <- outer(0:N1, 0:N2, function(i, j) {
  s <- i + j
  pmin(1, 2 * pmin(stats::phyper(i, N1, N2, s),
                   stats::phyper(i - 1, N1, N2, s, lower.tail = FALSE)))
})
max(abs(got - want))
#> [1] 6.661338e-16
```

The `blaker` convention is checked against formula (2) of Mehrotra, Chan
and Berger (2003) in the same way, and against `exact2x2` further below.

``` r

got <- attr(BinaryRR(N1, N2, 0.05, 'Fisher', alternative = 'two.sided',
                     tsmethod = 'blaker'), 'p.value')
want <- outer(0:N1, 0:N2, Vectorize(function(i, j) {
  s <- i + j
  k <- max(0, s - N2):min(N1, s)
  g <- pmin(stats::phyper(k, N1, N2, s),
            stats::phyper(k - 1, N1, N2, s, lower.tail = FALSE))
  min(1, sum(stats::dhyper(k, N1, N2, s)[g <= g[k == i] * (1 + 1e-10)]))
}))
max(abs(got - want))
#> [1] 0
```

## Unconditional tests against a direct evaluation

The unconditional p-value is the null tail probability maximized over
the nuisance parameter. The reference below forms the tail set of each
outcome explicitly and searches the same grid, which is slow but follows
the definition. It shares with the package only the Fisher p-values that
order the outcomes for the Boschloo test, and these are checked against
[`stats::fisher.test`](https://rdrr.io/r/stats/fisher.test.html) above.

``` r

unconditional_ref <- function(stat, N1, N2, n.grid, decreasing) {
  theta <- seq(0, 1, length.out = n.grid)
  joint <- lapply(theta, function(t) outer(dbinom(0:N1, N1, t), dbinom(0:N2, N2, t)))
  p <- matrix(0, nrow = N1 + 1, ncol = N2 + 1)
  for (i in 0:N1) {
    for (j in 0:N2) {
      s0 <- stat[i + 1, j + 1]
      tol <- 1e-10 * pmax(abs(stat), abs(s0))
      mask <- if (decreasing) stat >= s0 - tol else stat <= s0 + tol
      p[i + 1, j + 1] <- max(vapply(joint, function(m) sum(m * mask), numeric(1)))
    }
  }
  pmin(p, 1)
}

N1 <- 8
N2 <- 7
n.grid <- 40
boschloo <- attr(BinaryRR(N1, N2, 0.025, 'Boschloo', n.grid = n.grid), 'p.value')
fisher <- attr(BinaryRR(N1, N2, 0.025, 'Fisher'), 'p.value')
zpool <- attr(BinaryRR(N1, N2, 0.025, 'Z-pool', n.grid = n.grid), 'p.value')
zstats <- outer(0:N1, 0:N2, function(i, j) {
  hat.p <- (i + j) / (N1 + N2)
  z <- (i / N1 - j / N2) / sqrt(hat.p * (1 - hat.p) * (1 / N1 + 1 / N2))
  ifelse(is.finite(z), z, 0)
})
data.frame(
  Test = c('Boschloo', 'Z-pool'),
  max.absolute.difference = c(
    max(abs(boschloo - unconditional_ref(fisher, N1, N2, n.grid, FALSE))),
    max(abs(zpool - unconditional_ref(zstats, N1, N2, n.grid, TRUE)))
  )
)
#>       Test max.absolute.difference
#> 1 Boschloo            2.220446e-16
#> 2   Z-pool            1.110223e-16
```

## Why the tie handling matters

Version 1 of the package accumulated the null probabilities along an
arbitrary ordering of the outcomes, which splits groups of tied values.
The function below reproduces that behaviour so the two can be compared.

``` r

legacy_boschloo <- function(N1, N2, n.grid = 100) {
  stat <- attr(BinaryRR(N1, N2, 0.5, 'Fisher'), 'p.value')
  ord <- order(c(stat))
  x1 <- c(row(stat))[ord] - 1L
  x2 <- c(col(stat))[ord] - 1L
  theta <- seq(0, 1, length.out = n.grid)
  out <- numeric(length(ord))
  for (t in theta) {
    out <- pmax(out, cumsum(dbinom(x1, N1, t) * dbinom(x2, N2, t)))
  }
  p <- stat
  p[ord] <- pmin(1, out)
  p
}

N1 <- 7
N2 <- 7
current <- attr(BinaryRR(N1, N2, 0.025, 'Boschloo'), 'p.value')
legacy <- legacy_boschloo(N1, N2)
c(cells.with.different.p = sum(abs(current - legacy) > 1e-12),
  cells.with.different.decision = sum((current < 0.025) != (legacy < 0.025)))
#>        cells.with.different.p cells.with.different.decision 
#>                            29                             1
```

The outcomes $`(x_{1}, x_{2}) = (5, 1)`$ and $`(6, 2)`$ share a Fisher
p-value of $`2/39`$, so no test whose ordering statistic is the Fisher
p-value can distinguish them.

``` r

fisher77 <- attr(BinaryRR(7, 7, 0.5, 'Fisher'), 'p.value')
data.frame(
  outcome = c('x1 = 5, x2 = 1', 'x1 = 6, x2 = 2'),
  fisher.p = c(fisher77[6, 2], fisher77[7, 3]),
  legacy.p = c(legacy[6, 2], legacy[7, 3]),
  current.p = c(current[6, 2], current[7, 3])
)
#>          outcome   fisher.p   legacy.p  current.p
#> 1 x1 = 5, x2 = 1 0.05128205 0.02867981 0.02867981
#> 2 x1 = 6, x2 = 2 0.05128205 0.02161924 0.02867981
```

The legacy calculation gives the two outcomes different p-values and
rejects at one of them but not the other. The rejection region it
produces still holds the type I error rate below the nominal level. The
legacy p-values do not decrease along the ordering, so the region is an
initial segment of it, and at every value of the nuisance parameter on
the grid its rejection probability is the running total at its last
member, which is below alpha. The decision nevertheless depends on how
the sorting routine breaks the tie rather than on the data.

``` r

max_type1 <- function(reject, N1, N2, n.grid = 401) {
  theta <- seq(0, 1, length.out = n.grid)
  max(vapply(theta, function(t) {
    sum(outer(dbinom(0:N1, N1, t), dbinom(0:N2, N2, t)) * reject)
  }, numeric(1)))
}
c(legacy = max_type1(legacy < 0.025, 7, 7),
  current = max_type1(current < 0.025, 7, 7))
#>     legacy    current 
#> 0.02162352 0.01182923
```

## Comparison with the Exact package

``` r

N1 <- 10
N2 <- 10
cells <- expand.grid(x1 = c(7, 8, 9), x2 = c(1, 2))
compare_exact <- function(method, Test) {
  ours <- attr(BinaryRR(N1, N2, 0.025, Test, n.grid = 500), 'p.value')
  theirs <- vapply(seq_len(nrow(cells)), function(k) {
    tab <- matrix(c(cells$x1[k], cells$x2[k],
                    N1 - cells$x1[k], N2 - cells$x2[k]), nrow = 2)
    Exact::exact.test(tab, alternative = 'greater', method = method,
                      npNumbers = 500, to.plot = FALSE)$p.value
  }, numeric(1))
  max(abs(ours[cbind(cells$x1 + 1, cells$x2 + 1)] - theirs))
}
data.frame(
  Test = c('Boschloo', 'Z-pool'),
  max.absolute.difference = c(
    compare_exact('boschloo', 'Boschloo'),
    compare_exact('z-pooled', 'Z-pool')
  )
)
#>       Test max.absolute.difference
#> 1 Boschloo            2.109562e-07
#> 2   Z-pool            2.109562e-07
```

Residual differences come from the search over the nuisance parameter,
which the two packages carry out differently.

## Comparison with exact2x2

[`exact2x2::exact2x2()`](https://rdrr.io/pkg/exact2x2/man/exact2x2.html)
computes the two-sided Fisher p-value of the `blaker` convention. It is
compared with `tsmethod = 'blaker'` at every outcome with at least one
responder and at least one non-responder in total; the two remaining
outcomes have the p-value 1.

``` r

N1 <- 9
N2 <- 6
ours <- attr(BinaryRR(N1, N2, 0.05, 'Fisher', alternative = 'two.sided',
                      tsmethod = 'blaker'), 'p.value')
theirs <- outer(0:N1, 0:N2, Vectorize(function(i, j) {
  if (i + j == 0 || i + j == N1 + N2) return(1)
  tab <- matrix(c(i, j, N1 - i, N2 - j), nrow = 2)
  exact2x2::exact2x2(tab, tsmethod = 'blaker', conf.int = FALSE)$p.value
}))
max(abs(ours - theirs))
#> [1] 4.440892e-16
```

[`exact2x2::boschloo()`](https://rdrr.io/pkg/exact2x2/man/boschloo.html)
offers the `central` and `minlike` conventions. The two-sided Boschloo
p-values of the `central` convention are compared below.

``` r

N1 <- 10
N2 <- 10
ours <- attr(BinaryRR(N1, N2, 0.05, 'Boschloo', alternative = 'two.sided',
                      tsmethod = 'central', n.grid = 1000), 'p.value')
cells <- list(c(8, 2), c(7, 3), c(9, 1))
data.frame(
  outcome = vapply(cells, function(c) sprintf('x1 = %d, x2 = %d', c[1], c[2]), ''),
  bbssr = vapply(cells, function(c) ours[c[1] + 1, c[2] + 1], numeric(1)),
  exact2x2 = vapply(cells, function(c) {
    exact2x2::boschloo(c[1], N1, c[2], N2, alternative = 'two.sided',
                       tsmethod = 'central')$p.value
  }, numeric(1))
)
#>          outcome       bbssr     exact2x2
#> 1 x1 = 8, x2 = 2 0.011817889 0.0118135656
#> 2 x1 = 7, x2 = 3 0.115318107 0.1152988452
#> 3 x1 = 9, x2 = 1 0.000402448 0.0004021879
```

## Type I error rate of every test

``` r

tests <- c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo', 'Blackwelder',
           'Farrington-Manning')
type1 <- function(Test, alternative, N = 30, alpha = 0.025) {
  RR <- BinaryRR(N, N, alpha, Test, alternative = alternative, n.grid = 200)
  max_type1(matrix(as.vector(RR), N + 1, N + 1), N, N, n.grid = 801)
}
one <- vapply(tests, type1, numeric(1), 'greater')
two <- vapply(tests, type1, numeric(1), 'two.sided')
data.frame(
  Test = tests,
  one.sided = round(one, 5),
  two.sided = round(two, 5),
  exceeds.alpha = one > 0.025 | two > 0.025,
  row.names = NULL
)
#>                 Test one.sided two.sided exceeds.alpha
#> 1              Chisq   0.02770   0.02741          TRUE
#> 2             Fisher   0.01370   0.01350         FALSE
#> 3        Fisher-midP   0.02595   0.02735          TRUE
#> 4             Z-pool   0.02346   0.02293         FALSE
#> 5           Boschloo   0.02344   0.02267         FALSE
#> 6        Blackwelder   0.03020   0.03230          TRUE
#> 7 Farrington-Manning   0.02770   0.02741          TRUE
```

The nominal level is 0.025. In this configuration the three exact tests
stay below it for both alternatives. Their p-values are maximized over a
grid of the nuisance parameter, which can understate the exact p-value,
and `ref.pvalue = TRUE` refines the maximum when the size has to be
guaranteed, as the `bbssr-statistical-methods` vignette shows. The
chi-squared, mid-p, Blackwelder and Farrington-Manning tests carry no
such guarantee, and the `exceeds.alpha` column records where they
overshoot at this configuration. Without a margin the Farrington-Manning
test coincides with the chi-squared test, and the Blackwelder test is
the Wald test with the unpooled standard error.

## Speed

The whole rejection region is computed in one pass, so obtaining a power
curve costs little more than a single p-value. The inner loop over the
nuisance parameter runs in compiled code. The p-values of each test are
normally kept for the rest of the session, so the timings below turn
this off with `options(bbssr.cache = FALSE)` and measure the computation
itself.

``` r

old <- options(bbssr.cache = FALSE)
microbenchmark::microbenchmark(
  Chisq = BinaryRR(50, 50, 0.025, 'Chisq'),
  Fisher = BinaryRR(50, 50, 0.025, 'Fisher'),
  `Z-pool` = BinaryRR(50, 50, 0.025, 'Z-pool'),
  Boschloo = BinaryRR(50, 50, 0.025, 'Boschloo'),
  `Boschloo, Berger-Boos` = BinaryRR(50, 50, 0.025, 'Boschloo', bb.gamma = 1e-4),
  times = 10L,
  unit = 'ms'
)
#> Unit: milliseconds
#>                   expr      min       lq       mean    median        uq
#>                  Chisq 0.218223 0.223931  0.2400154  0.239084  0.255138
#>                 Fisher 1.797993 1.837241  1.8634245  1.862513  1.887435
#>                 Z-pool 2.968446 3.008305  3.0212621  3.017979  3.033202
#>               Boschloo 4.723965 4.750605  4.7784060  4.778391  4.805947
#>  Boschloo, Berger-Boos 9.943775 9.992627 10.0779836 10.034916 10.131923
#>        max neval
#>   0.266224    10
#>   1.959912    10
#>   3.080241    10
#>   4.836011    10
#>  10.409735    10
options(old)
```

The Berger-Boos variant costs more than the plain Boschloo test, because
the confidence bounds of every possible responder total are added to the
grid over the nuisance parameter.

A cell-by-cell comparison against a package that computes one p-value at
a time is not like for like, since
[`BinaryRR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryRR.md)
returns the p-values of every outcome of its grid. The comparison below
therefore charges `Exact` only for the p-values of a single row.

``` r

N1 <- 20
N2 <- 20
old <- options(bbssr.cache = FALSE)
microbenchmark::microbenchmark(
  bbssr.whole.grid = BinaryRR(N1, N2, 0.025, 'Boschloo'),
  Exact.one.row = for (x1 in 0:N1) {
    Exact::exact.test(matrix(c(x1, 5, N1 - x1, N2 - 5), nrow = 2),
                      alternative = 'greater', method = 'boschloo',
                      npNumbers = 100, to.plot = FALSE)
  },
  times = 5L,
  unit = 'ms'
)
#> Unit: milliseconds
#>              expr        min        lq       mean     median         uq
#>  bbssr.whole.grid   1.486792   1.62106   1.676154   1.630004   1.646988
#>     Exact.one.row 199.249984 199.56055 201.125234 200.649077 202.927380
#>         max neval
#>    1.995926     5
#>  203.239172     5
options(old)
```

## Sample size search

The sample size returned by
[`BinarySampleSize()`](https://gosukehommaex.github.io/bbssr/reference/BinarySampleSize.md)
is checked against a direct scan of the power function.

``` r

check <- function(Test) {
  ss <- BinarySampleSize(0.6, 0.2, 1, 0.025, 0.8, Test)
  data.frame(
    Test = Test,
    N2 = ss$N2,
    power.at.N2 = round(ss$Power, 4),
    power.at.N2.minus.1 = round(
      BinaryPower(0.6, 0.2, ss$N2 - 1, ss$N2 - 1, 0.025, Test)$Power, 4)
  )
}
do.call(rbind, lapply(c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo'), check))
#>          Test N2 power.at.N2 power.at.N2.minus.1
#> 1       Chisq 23      0.8191              0.7911
#> 2      Fisher 27      0.8024              0.7993
#> 3 Fisher-midP 23      0.8155              0.7851
#> 4      Z-pool 23      0.8088              0.7851
#> 5    Boschloo 23      0.8088              0.7911
```

The power reaches the target at the returned sample size and falls short
one patient per group below it.

The two other searches are checked against a scan of the power function
at every size of group 2 up to the limit of the stable search, for
Fisher’s exact test with a target power of 0.85, where they return
different sizes.

``` r

found <- lapply(c('smallest', 'stable'), function(s) {
  BinarySampleSize(0.6, 0.3, 1, 0.025, 0.85, 'Fisher', search = s)
})
limit <- attr(found[[2]], 'search.limit')
scan <- vapply(seq_len(limit), function(n) {
  BinaryPower(0.6, 0.3, n, n, 0.025, 'Fisher')$Power
}, numeric(1))
data.frame(
  search = c('smallest', 'stable'),
  N2 = c(found[[1]]$N2, found[[2]]$N2),
  N2.from.scan = c(min(which(scan >= 0.85)), max(which(scan < 0.85)) + 1)
)
#>     search N2 N2.from.scan
#> 1 smallest 52           52
#> 2   stable 56           56
```

`N2.from.scan` is the first size whose power attains the target, and one
more than the last size below the limit whose power falls short of it.
Both agree with the searches.

## Re-estimation designs

A design whose final sample size does not depend on the interim data is
a fixed-sample design. Bounding the final total sample size from below
and from above by the initial total turns a re-estimation design into
one, and its power must then equal that of the fixed-sample design.

``` r

fixed <- BinaryPowerBSSR(
  p = seq(0.2, 0.6, by = 0.1), Delta.A = 0.3, Delta.T = 0.3,
  N1 = 30, N2 = 30, n.interim = c(15, 15), r = 1, alpha = 0.025, tar.power = 0.8,
  Test = 'Z-pool', ss.method = 'standard', N.min = 60, N.max = 60
)
max(abs(fixed$power.BSSR - fixed$power.TRAD))
#> [1] 3.330669e-16
```

The type I error rate of a re-estimation design is computed in two
separate ways.
[`BinaryTypeIErrorBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryTypeIErrorBSSR.md)
sums the rejection probability over the interim outcomes and the
outcomes of the second stage, and
[`BinaryCondRejectBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryCondRejectBSSR.md)
averages the conditional rejection probabilities given the pooled
numbers of responders of the two stages over their binomial
distributions. The two agree for every test and alternative below.

``` r

theta <- c(0.1, 0.3, 0.5, 0.7, 0.9)
routes <- function(Test, alternative, Delta.A) {
  args <- list(Delta.A = Delta.A, N1 = 16, N2 = 16, n.interim = c(8, 8), r = 1,
               alpha = 0.025, tar.power = 0.8, Test = Test, alternative = alternative,
               ss.method = 'standard', theta = theta)
  crp <- do.call(BinaryCondRejectBSSR, args)
  tie <- do.call(BinaryTypeIErrorBSSR, c(args, list(maximize = 'grid')))
  max(abs(attr(crp, 'TIE')$TIE - tie$TIE.BSSR))
}
data.frame(
  Test = c('Chisq', 'Fisher', 'Z-pool', 'Boschloo'),
  alternative = c('greater', 'less', 'two.sided', 'greater'),
  max.absolute.difference = c(routes('Chisq', 'greater', 0.3),
                              routes('Fisher', 'less', -0.3),
                              routes('Z-pool', 'two.sided', 0.3),
                              routes('Boschloo', 'greater', 0.3))
)
#>       Test alternative max.absolute.difference
#> 1    Chisq     greater            4.510281e-17
#> 2   Fisher        less            3.122502e-17
#> 3   Z-pool   two.sided            1.734723e-17
#> 4 Boschloo     greater            3.122502e-17
```

The largest type I error rate is also obtained in two separate ways.
With the default `maximize = 'certified'`,
[`BinaryTypeIErrorBSSR()`](https://gosukehommaex.github.io/bbssr/reference/BinaryTypeIErrorBSSR.md)
computes the coefficients of the rate in the Bernstein basis and halves
the interval until an upper bound is within $`10^{-12}`$ of the largest
value found, as the `bbssr-statistical-methods` vignette describes. With
`theta = c(0, 1)` the grid holds only the two ends of the interval, so
the summation above is used only there and the largest value in between
comes from the coefficients. With `maximize = 'refined'` the function
evaluates the rate by that summation on a grid, here of step 0.001, and
refines the largest local maxima by a one-dimensional optimization.

``` r

cert.args <- list(Delta.A = 0.3, N1 = 39, N2 = 39, n.interim = c(20, 20), r = 1,
                  alpha = 0.025, tar.power = 0.8, Test = 'Chisq', ss.method = 'standard')
certified <- attr(do.call(BinaryTypeIErrorBSSR, c(cert.args, list(theta = c(0, 1)))),
                  'max')
dense <- attr(do.call(BinaryTypeIErrorBSSR,
                      c(cert.args, list(theta = seq(0, 1, by = 0.001),
                                        maximize = 'refined'))), 'max')
data.frame(Design = certified$Design, certified = certified$TIE, bound = certified$bound,
           dense.refined = dense$TIE, difference = certified$TIE - dense$TIE)
#>   Design  certified      bound dense.refined    difference
#> 1   BSSR 0.02728643 0.02728643    0.02728643 -3.105051e-13
#> 2   TRAD 0.02937611 0.02937611    0.02937611 -1.157720e-13
```

For both designs the two maxima differ by at most 3.1e-13, and the
refined maximum does not exceed the bound by more than rounding error.

## Reproduction of published results

The script `reproduce-published.R`, installed with the package in the
folder given by `system.file('reproduce', package = 'bbssr')`,
recomputes the published numerical results of Berger and Boos (1994),
Blackwelder (1982), Boschloo (1970), Farrington and Manning (1990), Fay
and Hunsberger (2021), Friede and Kieser (2004), Friede, Mitchell and
Mueller-Velten (2007), Kieser (2020) and Mehrotra, Chan and Berger
(2003), and lists every published value next to the recomputed one. It
is run from the root of the package source after `devtools::load_all()`,
writes its results to the folder `reproduce-output`, and runs for
several minutes. The values below are a selection that can be recomputed
quickly, each compared after rounding to the digits of the publication.

``` r

# Kieser (2020), Example 21.1: initial sample size by formula (21.3), and the largest
# type I error rate over the overall rate in [0.15, 0.85] for Delta = 0.30
k.plan <- BinarySampleSize(0.27, 0.42, 1, 0.025, 0.8, 'Chisq', alternative = 'less',
                           method = 'standard')
k.grid <- seq(0.01, 0.99, by = 0.005)
k.tie <- BinaryTypeIErrorBSSR(
  Delta.A = -0.30, N1 = 39, N2 = 39, n.interim = c(20, 20), r = 1, alpha = 0.025,
  tar.power = 0.8, Test = 'Chisq', alternative = 'less', ss.method = 'standard',
  theta = k.grid[k.grid >= 0.15 - 1e-9 & k.grid <= 0.85 + 1e-9], maximize = 'grid'
)
# Friede and Kieser (2004), Table I: control rate 0.2, alternative 0.2, equal groups and
# an internal pilot study of 40 patients
fk.plan <- BinarySampleSize(0.4, 0.2, 1, 0.05, 0.8, 'Chisq', alternative = 'two.sided',
                            method = 'null.variance', rounding = 'friede-kieser')
fk <- BinaryPowerBSSR(
  p = 0.3, Delta.A = 0.2, Delta.T = 0.2, N1 = fk.plan$N1, N2 = fk.plan$N2,
  n.interim = c(20, 20), r = 1, alpha = 0.05, tar.power = 0.8, Test = 'Chisq',
  alternative = 'two.sided', ss.method = 'null.variance', rounding = 'friede-kieser',
  N.max = 2 * fk.plan$N
)
# Farrington and Manning (1990), first example: p1 = 0.4, p2 = 0.05 and s0 = 0.2, so the
# margin is -0.2 in the parametrization of bbssr
fm.plan <- BinarySampleSize(0.4, 0.05, 1, 0.05, 0.8, 'Farrington-Manning',
                            method = 'standard', rounding = 'nearest', margin = -0.2)
fm.power <- BinaryPower(0.4, 0.05, 80, 80, 0.05, 'Farrington-Manning', margin = -0.2)
# Friede, Mitchell and Mueller-Velten (2007), Section 5: margin 0.1 and equal groups
fmm.bw <- BinarySampleSize(0.7, 0.7, 1, 0.025, 0.8, 'Blackwelder',
                           method = 'alternative.variance', rounding = 'friede-kieser',
                           margin = 0.1)
fmm.fm <- BinarySampleSize(0.9, 0.9, 1, 0.025, 0.8, 'Farrington-Manning',
                           method = 'standard', rounding = 'friede-kieser', margin = 0.1)
# Boschloo (1970), Section 4: 15 patients per group, one-sided level 0.05, p1 = 0.3 and
# p2 = 0.1, for Fisher's test and for Fisher's test at the raised level (Boschloo)
bo.power <- c(BinaryPower(0.3, 0.1, 15, 15, 0.05, 'Fisher')$Power,
              BinaryPower(0.3, 0.1, 15, 15, 0.05, 'Boschloo')$Power)
published <- data.frame(
  source = c('Kieser (2020)', 'Kieser (2020)', 'Kieser (2020)',
             'Friede and Kieser (2004)', 'Friede and Kieser (2004)',
             'Friede and Kieser (2004)', 'Farrington and Manning (1990)',
             'Farrington and Manning (1990)', 'Friede et al. (2007)',
             'Friede et al. (2007)', 'Boschloo (1970)', 'Boschloo (1970)'),
  quantity = c('initial total sample size, Delta = 0.15',
               'largest level, fixed design, Delta = 0.30',
               'largest level, re-estimation, Delta = 0.30',
               'initial total sample size', 'expected total sample size', 'power',
               'sample size per group', 'power at 80 patients per group',
               'total sample size, Blackwelder, p = 0.7',
               'total sample size, Farrington-Manning, p = 0.9',
               'power of the Fisher test, 15 per group, p1 = 0.3',
               'power of the Boschloo test, 15 per group, p1 = 0.3'),
  published = c(314, 0.0268, 0.0273, 166, 162.0, 0.794, 80, 0.813, 660, 310, 0.2558,
                0.3531),
  digits = c(0, 4, 4, 0, 1, 3, 0, 3, 0, 0, 4, 4),
  recomputed = c(k.plan$N, max(k.tie$TIE.TRAD), max(k.tie$TIE.BSSR), fk.plan$N, fk$E.N,
                 fk$power.BSSR, fm.plan$N1, fm.power$Power, fmm.bw$N, fmm.fm$N,
                 bo.power)
)
published$agrees <- abs(round(published$recomputed, published$digits) -
                          published$published) < 1e-9
published
#>                           source
#> 1                  Kieser (2020)
#> 2                  Kieser (2020)
#> 3                  Kieser (2020)
#> 4       Friede and Kieser (2004)
#> 5       Friede and Kieser (2004)
#> 6       Friede and Kieser (2004)
#> 7  Farrington and Manning (1990)
#> 8  Farrington and Manning (1990)
#> 9           Friede et al. (2007)
#> 10          Friede et al. (2007)
#> 11               Boschloo (1970)
#> 12               Boschloo (1970)
#>                                              quantity published digits
#> 1             initial total sample size, Delta = 0.15  314.0000      0
#> 2           largest level, fixed design, Delta = 0.30    0.0268      4
#> 3          largest level, re-estimation, Delta = 0.30    0.0273      4
#> 4                           initial total sample size  166.0000      0
#> 5                          expected total sample size  162.0000      1
#> 6                                               power    0.7940      3
#> 7                               sample size per group   80.0000      0
#> 8                      power at 80 patients per group    0.8130      3
#> 9             total sample size, Blackwelder, p = 0.7  660.0000      0
#> 10     total sample size, Farrington-Manning, p = 0.9  310.0000      0
#> 11   power of the Fisher test, 15 per group, p1 = 0.3    0.2558      4
#> 12 power of the Boschloo test, 15 per group, p1 = 0.3    0.3531      4
#>      recomputed agrees
#> 1  314.00000000   TRUE
#> 2    0.02679047   TRUE
#> 3    0.02728569   TRUE
#> 4  166.00000000   TRUE
#> 5  162.04925552   TRUE
#> 6    0.79416905   TRUE
#> 7   80.00000000   TRUE
#> 8    0.81320091   TRUE
#> 9  660.00000000   TRUE
#> 10 310.00000000   TRUE
#> 11   0.25576535   TRUE
#> 12   0.35314961   TRUE
```

Of the 12 values above, 12 agree with the publication after rounding.
The script reports every value that it does not reproduce, together with
the reason when one is known. The vignette ‘Validation by Reproducing
Published Figures’ redraws twelve figures of the same publications from
values computed with the package.

The two-sided Fisher p-values of the conventions are compared with two
published examples: 8 responders out of 14 against 1 out of 7 in Fay and
Hunsberger (2021, Section 8), and 8 out of 148 against 1 out of 132 in
Mehrotra, Chan and Berger (2003, Section 3.1).

``` r

fisher_p <- function(N1, N2, x1, x2, tsmethod) {
  attr(BinaryRR(N1, N2, 0.05, 'Fisher', alternative = 'two.sided', tsmethod = tsmethod),
       'p.value')[x1 + 1, x2 + 1]
}
two.sided <- data.frame(
  source = c(rep('Fay and Hunsberger (2021)', 3), rep('Mehrotra et al. (2003)', 2)),
  outcome = c(rep('8 / 14 against 1 / 7', 3), rep('8 / 148 against 1 / 132', 2)),
  tsmethod = c('minlike', 'central', 'blaker', 'minlike', 'blaker'),
  published = c(0.159, 0.157, 0.087, 0.0388, 0.0388),
  digits = c(3, 3, 3, 4, 4),
  recomputed = c(fisher_p(14, 7, 8, 1, 'minlike'), fisher_p(14, 7, 8, 1, 'central'),
                 fisher_p(14, 7, 8, 1, 'blaker'), fisher_p(148, 132, 8, 1, 'minlike'),
                 fisher_p(148, 132, 8, 1, 'blaker'))
)
two.sided$agrees <- abs(round(two.sided$recomputed, two.sided$digits) -
                          two.sided$published) < 1e-9
two.sided
#>                      source                 outcome tsmethod published digits
#> 1 Fay and Hunsberger (2021)    8 / 14 against 1 / 7  minlike    0.1590      3
#> 2 Fay and Hunsberger (2021)    8 / 14 against 1 / 7  central    0.1570      3
#> 3 Fay and Hunsberger (2021)    8 / 14 against 1 / 7   blaker    0.0870      3
#> 4    Mehrotra et al. (2003) 8 / 148 against 1 / 132  minlike    0.0388      4
#> 5    Mehrotra et al. (2003) 8 / 148 against 1 / 132   blaker    0.0388      4
#>   recomputed agrees
#> 1 0.15882353   TRUE
#> 2 0.15665635   TRUE
#> 3 0.08730650   TRUE
#> 4 0.03878174   TRUE
#> 5 0.03878174   TRUE
```

In the example of Mehrotra, Chan and Berger (2003) the `minlike` and
`blaker` conventions give the same p-value, whereas the `central`
convention gives 0.0543. Their text defines the two-sided Fisher p-value
by formula (2), the `blaker` convention. The reproduction script shows
that their Tables 1 and 3 agree with the `minlike` convention instead.
The two conventions coincide when the groups are of equal size, and
wherever they give different rounded values in the configurations with
groups of unequal size, the published type I error rates and powers of
the Fisher and Boschloo tests are those of `minlike`. The script
therefore compares these tables with `minlike`.

## Summary

The Fisher exact test reproduces
[`stats::fisher.test`](https://rdrr.io/r/stats/fisher.test.html) to
machine precision, and its `blaker` convention reproduces `exact2x2`.
The unconditional tests reproduce a direct evaluation of their
definition to machine precision, and agree with `Exact` and `exact2x2`
up to the difference in the search over the nuisance parameter. In the
configuration examined, the three exact tests hold the type I error rate
below the nominal level for one-sided and two-sided alternatives alike.
A re-estimation design without re-estimation reproduces the fixed-sample
design, the two calculations of the type I error rate of a re-estimation
design agree, the certified largest type I error rate is confirmed by a
dense grid, and published values are recomputed.

## References

Berger, R. L. and Boos, D. D. (1994). P values maximized over a
confidence set for the nuisance parameter. *Journal of the American
Statistical Association*, 89, 1012-1016.

Blackwelder, W. C. (1982). “Proving the null hypothesis” in clinical
trials. *Controlled Clinical Trials*, 3, 345-353.

Boschloo, R. D. (1970). Raised conditional level of significance for the
2 × 2-table when testing the equality of two probabilities. *Statistica
Neerlandica*, 24, 1-9.

Farrington, C. P. and Manning, G. (1990). Test statistics and sample
size formulae for comparative binomial trials with null hypothesis of
non-zero risk difference or non-unity relative risk. *Statistics in
Medicine*, 9, 1447-1454.

Fay, M. P. and Hunsberger, S. A. (2021). Practical valid inferences for
the two-sample binomial problem. *Statistics Surveys*, 15, 72-110.

Friede, T. and Kieser, M. (2004). Sample size recalculation for binary
data in internal pilot study designs. *Pharmaceutical Statistics*, 3,
269-279.

Friede, T., Mitchell, C. and Mueller-Velten, G. (2007). Blinded sample
size reestimation in non-inferiority trials with binary endpoints.
*Biometrical Journal*, 49, 903-916.

Kieser, M. (2020). *Methods and Applications of Sample Size Calculation
and Recalculation in Clinical Trials*. Springer.

Mehrotra, D. V., Chan, I. S. F. and Berger, R. L. (2003). A cautionary
note on exact unconditional inference for a difference between two
independent binomial proportions. *Biometrics*, 59, 441-450.
