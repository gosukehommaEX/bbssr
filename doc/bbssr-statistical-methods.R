## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5,
  fig.align = "center"
)

## ----load---------------------------------------------------------------------
library(bbssr)

## ----central-identity---------------------------------------------------------
N1 <- 9
N2 <- 7
alpha <- 0.02
two <- BinaryRR(N1, N2, 2 * alpha, 'Fisher',
                alternative = 'two.sided', tsmethod = 'central')
upper <- BinaryRR(N1, N2, alpha, 'Fisher')
lower <- t(BinaryRR(N2, N1, alpha, 'Fisher'))
identical(as.vector(two), as.vector(upper | lower))

## ----ts-compare---------------------------------------------------------------
data.frame(
  tsmethod = c('minlike', 'central'),
  rejected = c(
    sum(BinaryRR(N1, N2, 0.05, 'Fisher', alternative = 'two.sided',
                 tsmethod = 'minlike')),
    sum(BinaryRR(N1, N2, 0.05, 'Fisher', alternative = 'two.sided',
                 tsmethod = 'central'))
  )
)

## ----grid---------------------------------------------------------------------
coarse <- attr(BinaryRR(12, 12, 0.025, 'Boschloo', n.grid = 20), 'p.value')
fine <- attr(BinaryRR(12, 12, 0.025, 'Boschloo', n.grid = 2000), 'p.value')
data.frame(all.p.values.increased = all(fine >= coarse - 1e-12),
           largest.increase = max(fine - coarse))

## ----ties---------------------------------------------------------------------
stat <- attr(BinaryRR(7, 7, 0.025, 'Fisher'), 'p.value')
c(cell_5_1 = stat[6, 2], cell_6_2 = stat[7, 3])

## ----ties-equal---------------------------------------------------------------
p <- attr(BinaryRR(7, 7, 0.025, 'Boschloo'), 'p.value')
c(cell_5_1 = p[6, 2], cell_6_2 = p[7, 3])

## ----berger-boos--------------------------------------------------------------
plain <- BinaryRR(15, 15, 0.025, 'Boschloo', n.grid = 200)
bb <- BinaryRR(15, 15, 0.025, 'Boschloo', n.grid = 200, bb.gamma = 1e-4)
p.plain <- attr(plain, 'p.value')
p.bb <- attr(bb, 'p.value')
data.frame(rejected.plain = sum(plain), rejected.berger.boos = sum(bb),
           largest.decrease = max(p.plain - p.bb),
           largest.increase = max(p.bb - p.plain))

## ----nesting------------------------------------------------------------------
fisher <- BinaryRR(15, 15, 0.025, 'Fisher')
boschloo <- BinaryRR(15, 15, 0.025, 'Boschloo', n.grid = 200)
data.frame(rejected.fisher = sum(fisher), rejected.boschloo = sum(boschloo),
           fisher.region.contained = all(as.vector(boschloo)[as.vector(fisher)]))

## ----type1--------------------------------------------------------------------
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

## ----bssr-rules---------------------------------------------------------------
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

