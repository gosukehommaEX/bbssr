ss <- function(...) {
  args <- list(p1 = c(0.5, 0.7), p2 = c(0.3, 0.5), frac = c(0.4, 0.6), alpha = 0.025,
               tar.power = 0.8, Test = 'CMH', method = 'standard')
  new <- list(...)
  args[names(new)] <- new
  do.call(StratSampleSize, args)
}

test_that("StratSampleSize agrees with the reference values", {
  expect_equal(c(attr(ss(), 'n'), attr(ss(Test = 'CMH-cc'), 'n')), c(186, 206))
  ex <- function(Test, ...) {
    attr(ss(p1 = c(0.65, 0.8), p2 = c(0.3, 0.45), frac = c(0.5, 0.5), Test = Test,
            method = 'exact', ...), 'n')
  }
  expect_equal(c(ex('CMH'), ex('Strat-Fisher')), c(57, 69))
  expect_equal(ex('CMH', search = 'smallest'), 57)
  expect_equal(attr(ss(), 'n.raw'), 185.997689655091, tolerance = 1e-10)
})

test_that("StratSampleSize returns the group sizes and their exact power", {
  res <- ss(method = 'exact')
  expect_s3_class(res, 'bbssr_ssamplesize')
  expect_named(res, c('p1', 'p2', 'frac', 'r', 'alpha', 'tar.power', 'Test', 'alternative',
                      'N1', 'N2', 'N', 'Power'))
  sz <- strat_sizes(attr(res, 'n'), c(0.4, 0.6), c(1, 1))
  expect_identical(res$N1, as.integer(sz$N1))
  expect_identical(res$N, as.integer(sum(sz$N1) + sum(sz$N2)))
  expect_equal(res$Power, StratPower(c(0.5, 0.7), c(0.3, 0.5), res$N1, res$N2, 0.025,
                                     'CMH')$Power)
  expect_gte(res$Power, 0.8)
  expect_equal(res$frac, c(0.4, 0.6))
})

test_that("StratSampleSize checks its arguments", {
  expect_error(ss(p2 = c(0.5, 0.7)), 'no effect')
  expect_error(ss(p1 = rbind(c(0.5, 0.7), c(0.5, 0.7))), 'one element per stratum')
  expect_error(ss(tar.power = 1), 'tar.power')
  expect_error(ss(frac = c(0.4, 0.3, 0.3)), 'one element per stratum')
  expect_error(ss(p1 = 0.5, p2 = 0.3, frac = 1), 'at least two strata')
})
