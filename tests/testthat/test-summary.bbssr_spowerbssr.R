# Reference values from tools/reference/reference_strat.py (blocks_summary_sbssr), which
# enumerates the interim pooled numbers of responders and re-estimates the final group
# sizes independently of the package. For each design: E.N, SD.N, N.q25, N.q50, N.q75 and
# P.N.max with the true effect 0.3 (p1 = (0.55, 0.75), p2 = (0.25, 0.45)), then without
# an effect
ssum <- function(Delta.T, ...) {
  args <- list(p = c(0.4, 0.6), Delta.A = 0.3, Delta.T = Delta.T, N1 = c(10, 10),
               N2 = c(10, 10), omega = 0.5, r = 1, alpha = 0.05, tar.power = 0.8,
               Test = 'CMH', ss.method = 'standard')
  new <- list(...)
  args[names(new)] <- new
  s <- summary(do.call(StratPowerBSSR, args))
  unlist(s[1, c('E.N', 'SD.N', 'N.q25', 'N.q50', 'N.q75', 'P.N.max')], use.names = FALSE)
}

test_that("summary.bbssr_spowerbssr agrees with the reference values", {
  expect_equal(c(ssum(0.3), ssum(0)),
               c(59.740212490731, 7.12815613980113, 56, 60, 64, 0.198521899360162,
                 59.1902574387838, 7.52672540880168, 56, 60, 64, 0.185660316278052),
               tolerance = 1e-10)
  fisher <- function(Delta.T) {
    ssum(Delta.T, Delta.A = 0.4, N1 = c(8, 10), N2 = c(8, 10), alpha = 0.1,
         Test = 'Strat-Fisher', ss.method = 'exact', N.max = 48)
  }
  # P.N.max is the probability of the largest final total of the table of the
  # re-estimation, 48 in this design
  expect_equal(c(fisher(0.3), fisher(0)),
               c(32.6818367739194, 2.4998983667878, 32, 34, 34, 5.90356207233369e-05,
                 32.6534657362278, 2.63854503937714, 32, 34, 34, 0.000107352592482304),
               tolerance = 1e-10)
})

test_that("summary.bbssr_spowerbssr reports the distribution of every scenario", {
  res <- StratPowerBSSR(p = rbind(c(0.4, 0.6), c(0.3, 0.5)), Delta.A = 0.3, Delta.T = 0.3,
                        N1 = c(10, 10), N2 = c(10, 10), omega = 0.5, r = 1, alpha = 0.05,
                        tar.power = 0.8, Test = 'CMH', ss.method = 'standard')
  s <- summary(res, probs = c(0.1, 0.9))
  expect_named(s, c('p.S1', 'p.S2', 'p1.S1', 'p1.S2', 'p2.S1', 'p2.S2', 'E.N', 'SD.N',
                    'N.q10', 'N.q90', 'P.N.max'))
  expect_equal(s$E.N, res$E.N, tolerance = 1e-12)
  expect_identical(s$p1.S2, res$p1.S2)
  d <- attr(res, 'N.dist')
  expect_equal(as.vector(tapply(d$prob, d$scenario, sum)), c(1, 1), tolerance = 1e-12)
  expect_true(all(s$N.q10 <= s$N.q90))
  top <- max(attr(res, 'reestimation')$N)
  expect_equal(s$P.N.max[2], sum(d$prob[d$scenario == 2 & d$N == top]))
  expect_error(summary(res, probs = 0), 'probs')
  bare <- structure(data.frame(p.S1 = 0.4), class = c('bbssr_spowerbssr', 'data.frame'))
  expect_error(summary(bare), 'no distribution')
})

test_that("P.N.max refers to N.max without a table of the re-estimation", {
  run <- function(N.max) {
    set.seed(5)
    StratPowerBSSR(p = c(0.4, 0.6), Delta.A = 0.3, Delta.T = 0.3, N1 = c(6, 6),
                   N2 = c(6, 6), omega = 0.5, r = 1, alpha = 0.05, tar.power = 0.8,
                   Test = 'CMH', ss.method = 'standard', alloc = 'stratum', nseq = 5,
                   N.max = N.max)
  }
  res <- run(40)
  s <- summary(res)
  d <- attr(res, 'N.dist')
  expect_null(attr(res, 'reestimation'))
  expect_equal(s$P.N.max, sum(d$prob[d$N == 40]))
  expect_gt(s$P.N.max, 0)
  expect_equal(s$E.N, res$E.N, tolerance = 1e-12)
  expect_true(is.na(summary(run(NULL))$P.N.max))
})
