sb <- function(...) {
  args <- list(n1 = c(10, 12), n2 = c(10, 12), S = c(9, 16), Delta.A = 0.2, r = 1,
               alpha = 0.025, tar.power = 0.8, Test = 'CMH', N1 = c(20, 24),
               N2 = c(20, 24), ss.method = 'standard')
  new <- list(...)
  args[names(new)] <- new
  do.call(StratBSSR, args)
}

test_that("StratBSSR agrees with the reference values", {
  re <- sb()
  expect_equal(c(re$N1.final, re$N2.final), c(42, 50, 42, 50))
  expect_equal(attr(re, 'Power'), 0.809417719209081, tolerance = 1e-10)
  expect_equal(re$hat.p1, c(9 / 20 + 0.1, 16 / 24 + 0.1))
  expect_identical(re$n1.stage2, re$N1.final - c(10L, 12L))
  expect_identical(attr(re, 'N.final'), as.integer(sum(re$N1.final) + sum(re$N2.final)))
  expect_s3_class(re, 'bbssr_sbssr')
})

test_that("StratBSSR agrees with the re-estimation of StratPowerBSSR", {
  res <- StratPowerBSSR(p = c(0.5, 0.6), Delta.A = 0.2, Delta.T = 0.2, N1 = c(20, 24),
                        N2 = c(20, 24), n.interim = rbind(c(10, 12), c(10, 12)), r = 1,
                        alpha = 0.025, tar.power = 0.8, Test = 'CMH',
                        ss.method = 'standard')
  map <- attr(res, 'reestimation')
  row <- map[map$s1 == 9 & map$s2 == 16, ]
  re <- sb()
  expect_equal(c(row$N1.S1, row$N1.S2, row$N2.S1, row$N2.S2),
               c(re$N1.final, re$N2.final))
})

test_that("StratBSSR checks its arguments", {
  expect_error(sb(S = c(9, 30)), 'S must')
  expect_error(sb(N1 = NULL, N2 = NULL), 'supply frac')
  expect_equal(sb(N1 = NULL, N2 = NULL, frac = c(20, 24))$N1.final, c(42, 50))
  expect_error(sb(N1 = c(20, 24, 10), N2 = c(20, 24, 10)), 'one element per stratum')
})
