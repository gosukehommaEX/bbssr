## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5,
  fig.align = "center",
  warning = FALSE,
  message = FALSE
)

## ----load---------------------------------------------------------------------
library(bbssr)
has_exact <- requireNamespace('Exact', quietly = TRUE)
has_exact2x2 <- requireNamespace('exact2x2', quietly = TRUE)
has_bench <- requireNamespace('microbenchmark', quietly = TRUE)
c(Exact = has_exact, exact2x2 = has_exact2x2, microbenchmark = has_bench)

## ----fisher-vs-stats----------------------------------------------------------
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

## ----central-vs-definition----------------------------------------------------
got <- attr(BinaryRR(N1, N2, 0.05, 'Fisher', alternative = 'two.sided',
                     tsmethod = 'central'), 'p.value')
want <- outer(0:N1, 0:N2, function(i, j) {
  s <- i + j
  pmin(1, 2 * pmin(stats::phyper(i, N1, N2, s),
                   stats::phyper(i - 1, N1, N2, s, lower.tail = FALSE)))
})
max(abs(got - want))

## ----brute-force--------------------------------------------------------------
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

## ----legacy-------------------------------------------------------------------
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

## ----tie-detail---------------------------------------------------------------
fisher77 <- attr(BinaryRR(7, 7, 0.5, 'Fisher'), 'p.value')
data.frame(
  outcome = c('x1 = 5, x2 = 1', 'x1 = 6, x2 = 2'),
  fisher.p = c(fisher77[6, 2], fisher77[7, 3]),
  legacy.p = c(legacy[6, 2], legacy[7, 3]),
  current.p = c(current[6, 2], current[7, 3])
)

## ----tie-type1----------------------------------------------------------------
max_type1 <- function(reject, N1, N2, n.grid = 401) {
  theta <- seq(0, 1, length.out = n.grid)
  max(vapply(theta, function(t) {
    sum(outer(dbinom(0:N1, N1, t), dbinom(0:N2, N2, t)) * reject)
  }, numeric(1)))
}
c(legacy = max_type1(legacy < 0.025, 7, 7),
  current = max_type1(current < 0.025, 7, 7))

## ----exact-comparison, eval = has_exact---------------------------------------
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

## ----exact2x2-comparison, eval = has_exact2x2---------------------------------
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

## ----type1-table--------------------------------------------------------------
tests <- c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')
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

## ----timing, eval = has_bench-------------------------------------------------
microbenchmark::microbenchmark(
  Chisq = BinaryRR(50, 50, 0.025, 'Chisq'),
  Fisher = BinaryRR(50, 50, 0.025, 'Fisher'),
  `Z-pool` = BinaryRR(50, 50, 0.025, 'Z-pool'),
  Boschloo = BinaryRR(50, 50, 0.025, 'Boschloo'),
  `Boschloo, Berger-Boos` = BinaryRR(50, 50, 0.025, 'Boschloo', bb.gamma = 1e-4),
  times = 10L,
  unit = 'ms'
)

## ----timing-exact, eval = has_exact && has_bench------------------------------
N1 <- 20
N2 <- 20
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

## ----samplesize-check---------------------------------------------------------
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

