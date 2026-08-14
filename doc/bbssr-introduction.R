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

## ----tests--------------------------------------------------------------------
tests <- c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')

## ----rr-----------------------------------------------------------------------
RR <- BinaryRR(N1 = 10, N2 = 10, alpha = 0.025, Test = 'Boschloo')
RR

## ----rr-plot------------------------------------------------------------------
plot(RR)

## ----rr-compare---------------------------------------------------------------
n_rejected <- vapply(tests, function(tst) {
  sum(BinaryRR(10, 10, 0.025, tst))
}, numeric(1))
n_rejected

## ----power--------------------------------------------------------------------
BinaryPower(p1 = 0.6, p2 = 0.3, N1 = 40, N2 = 40, alpha = 0.025, Test = 'Fisher')

## ----power-curve--------------------------------------------------------------
pw <- BinaryPower(p1 = seq(0.35, 0.75, by = 0.05), p2 = rep(0.3, 9),
                  N1 = 40, N2 = 40, alpha = 0.025, Test = 'Fisher')
plot(pw)

## ----samplesize---------------------------------------------------------------
ss <- BinarySampleSize(p1 = 0.6, p2 = 0.3, r = 1, alpha = 0.025,
                       tar.power = 0.8, Test = 'Fisher')
ss

## ----sawtooth-----------------------------------------------------------------
n2 <- (ss$N2 - 6):(ss$N2 + 5)
data.frame(
  N2 = n2,
  Power = round(vapply(n2, function(n) {
    BinaryPower(0.6, 0.3, n, n, 0.025, 'Fisher')$Power
  }, numeric(1)), 4)
)

## ----bssr-plan----------------------------------------------------------------
BinarySampleSize(p1 = 0.45, p2 = 0.09, r = 1, alpha = 0.025,
                 tar.power = 0.8, Test = 'Z-pool')

## ----bssr-power---------------------------------------------------------------
res <- BinaryPowerBSSR(
  p = seq(0.19, 0.37, by = 0.03),
  Delta.A = 0.36, Delta.T = 0.36,
  N1 = 24, N2 = 24, omega = 0.5, r = 1,
  alpha = 0.025, tar.power = 0.8, Test = 'Z-pool'
)
res

## ----bssr-plot----------------------------------------------------------------
plot(res)

## ----rules--------------------------------------------------------------------
args <- list(
  p = seq(0.19, 0.37, by = 0.06),
  Delta.A = 0.36, Delta.T = 0.36, N1 = 24, N2 = 24, omega = 0.5, r = 1,
  alpha = 0.025, tar.power = 0.8, Test = 'Chisq'
)
unrestricted <- do.call(BinaryPowerBSSR, c(args, list(restricted = FALSE)))
restricted <- do.call(BinaryPowerBSSR, c(args, list(restricted = TRUE)))
data.frame(
  p = unrestricted$p,
  power.unrestricted = round(unrestricted$power.BSSR, 4),
  power.restricted = round(restricted$power.BSSR, 4),
  EN.unrestricted = round(unrestricted$E.N, 1),
  EN.restricted = round(restricted$E.N, 1)
)

## ----bssr-interim-------------------------------------------------------------
BinaryBSSR(n1 = 12, n2 = 12, S = 8, Delta.A = 0.36, r = 1,
           alpha = 0.025, tar.power = 0.8, Test = 'Z-pool')

## ----two-sided----------------------------------------------------------------
pw <- function(alpha, alternative, tsmethod = 'minlike') {
  BinaryPower(0.6, 0.3, 40, 40, alpha, 'Fisher',
              alternative = alternative, tsmethod = tsmethod)$Power
}
round(c(
  one.sided.alpha.0.025 = pw(0.025, 'greater'),
  two.sided.alpha.0.025.minlike = pw(0.025, 'two.sided', 'minlike'),
  two.sided.alpha.0.025.central = pw(0.025, 'two.sided', 'central'),
  two.sided.alpha.0.05.minlike = pw(0.05, 'two.sided', 'minlike'),
  two.sided.alpha.0.05.central = pw(0.05, 'two.sided', 'central')
), 4)

