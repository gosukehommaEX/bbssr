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

test_that("StratSampleSize of the unstratified tests agrees with the reference values", {
  expect_equal(c(attr(ss(Test = 'Unstrat-Chisq', r = c(1, 2)), 'n'),
                 attr(ss(Test = 'Unstrat-Chisq-cc', r = c(1, 2)), 'n')), c(148, 165))
  ex <- function(Test) {
    attr(ss(p1 = c(0.65, 0.8), p2 = c(0.3, 0.45), frac = c(0.5, 0.5), Test = Test,
            method = 'exact'), 'n')
  }
  expect_equal(c(ex('Unstrat-Chisq'), ex('Unstrat-Fisher')), c(57, 73))
  # The normal approximation is that of the table of totals
  expect_equal(attr(ss(Test = 'Unstrat-Chisq', r = c(1, 2)), 'n.raw'),
               strat_raw_n(c(0.5, 0.7), c(0.3, 0.5), c(0.4, 0.6), c(1, 2), 0.025, 0.8,
                           'greater', FALSE, TRUE))
  expect_error(ss(p1 = c(0.6, 0.2), p2 = c(0.4, 0.4), Test = 'Unstrat-Chisq'), 'averages')
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

test_that("StratSampleSize with a random allocation agrees with the reference values", {
  # Reference values from tools/reference/reference_strat.py (alloc_exact_n): the search
  # 'crossing' over the total sample size with the power of a random allocation
  ra <- function(alloc, ...) {
    StratSampleSize(c(0.6, 0.8), c(0.2, 0.4), frac = c(0.5, 0.5), alpha = 0.05,
                    tar.power = 0.8, Test = 'CMH', alloc = alloc, ...)
  }
  st <- ra('stratum')
  gr <- ra('group')
  expect_equal(c(st$N, gr$N), c(36, 37))
  expect_null(st$N1)
  expect_null(st$N2)
  expect_identical(attr(st, 'n'), st$N)
  expect_equal(st$Power, StratPower(c(0.6, 0.8), c(0.2, 0.4), alpha = 0.05, Test = 'CMH',
                                    alloc = 'stratum', N = st$N,
                                    frac = c(0.5, 0.5))$Power)
  expect_gte(st$Power, 0.8)
  expect_identical(attr(gr, 'n.strata'),
                   as.integer(c(floor(gr$N / 2), gr$N - floor(gr$N / 2))))
  # The normal approximation gives the total sample size itself
  sd <- ra('both', method = 'standard')
  expect_identical(sd$N, as.integer(ceiling(attr(sd, 'n.raw'))))
  expect_identical(attr(sd, 'alloc'), 'both')
  expect_error(ra('stratum', block = 3), 'whole numbers')
  expect_error(StratSampleSize(c(0.6, 0.8, 0.7), c(0.2, 0.4, 0.3), frac = c(1, 1, 1),
                               alpha = 0.05, tar.power = 0.8, Test = 'CMH',
                               alloc = 'group'), 'two strata')
})
