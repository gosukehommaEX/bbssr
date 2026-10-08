test_that("BinaryBSSR returns a bbssr_bssr data frame", {
  res <- BinaryBSSR(n1 = 20, n2 = 20, S = 11, Delta.A = 0.3, r = 1,
                    alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
  expect_s3_class(res, 'bbssr_bssr')
  expect_equal(nrow(res), 1L)
  expect_named(res, c('n1', 'n2', 'n', 'S', 'hat.p', 'hat.p1', 'hat.p2',
                      'N1.re', 'N2.re', 'N.re',
                      'n1.stage2', 'n2.stage2', 'n.stage2',
                      'N1.final', 'N2.final', 'N.final', 'Power'))
})

test_that("the reported sample sizes are internally consistent", {
  res <- BinaryBSSR(n1 = 15, n2 = 15, S = 9, Delta.A = 0.25, r = 1,
                    alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
  expect_equal(res$n, res$n1 + res$n2)
  expect_equal(res$N.re, res$N1.re + res$N2.re)
  expect_equal(res$n.stage2, res$n1.stage2 + res$n2.stage2)
  expect_equal(res$N1.final, res$n1 + res$n1.stage2)
  expect_equal(res$N2.final, res$n2 + res$n2.stage2)
  expect_equal(res$N.final, res$N1.final + res$N2.final)
  expect_gte(res$n1.stage2, 0)
  expect_gte(res$n2.stage2, 0)
})

test_that("the blinded proportions are recovered from the pooled data", {
  res <- BinaryBSSR(n1 = 20, n2 = 20, S = 12, Delta.A = 0.3, r = 1,
                    alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
  expect_equal(res$hat.p, 12 / 40)
  expect_equal(res$hat.p1, 12 / 40 + 0.15)
  expect_equal(res$hat.p2, 12 / 40 - 0.15)
  expect_equal(res$hat.p1 - res$hat.p2, 0.3)
})

test_that("the recovered proportions stay inside the unit interval", {
  res <- BinaryBSSR(n1 = 10, n2 = 10, S = 20, Delta.A = 0.4, r = 1,
                    alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
  expect_equal(res$hat.p1, 1)
  res <- BinaryBSSR(n1 = 10, n2 = 10, S = 0, Delta.A = 0.4, r = 1,
                    alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
  expect_equal(res$hat.p2, 0)
})

test_that("the re-estimated sample size matches BinarySampleSize", {
  res <- BinaryBSSR(n1 = 18, n2 = 18, S = 10, Delta.A = 0.3, r = 1,
                    alpha = 0.025, tar.power = 0.8, Test = 'Fisher')
  ss <- BinarySampleSize(res$hat.p1, res$hat.p2, 1, 0.025, 0.8, 'Fisher')
  expect_equal(res$N1.re, ss$N1)
  expect_equal(res$N2.re, ss$N2)
})

test_that("the restricted rule never shrinks the trial below the plan", {
  res <- BinaryBSSR(n1 = 30, n2 = 30, S = 30, Delta.A = 0.3, r = 1,
                    alpha = 0.025, tar.power = 0.8, Test = 'Chisq',
                    restricted = TRUE, N1 = 90, N2 = 90)
  expect_gte(res$N2.final, 90)
  unres <- BinaryBSSR(n1 = 30, n2 = 30, S = 30, Delta.A = 0.3, r = 1,
                      alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
  expect_lte(unres$N.final, res$N.final)
})

test_that("the allocation ratio is respected in the second stage", {
  res <- BinaryBSSR(n1 = 20, n2 = 10, S = 9, Delta.A = 0.3, r = 2,
                    alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
  expect_equal(res$n1.stage2, as.integer(ceiling(2 * res$n2.stage2)))
})

test_that("BinaryBSSR validates its arguments", {
  expect_error(BinaryBSSR(20, 20, 41, 0.3, 1, 0.025, 0.8, 'Chisq'), 'between 0 and')
  expect_error(BinaryBSSR(20, 20, -1, 0.3, 1, 0.025, 0.8, 'Chisq'), 'between 0 and')
  expect_error(BinaryBSSR(0, 20, 5, 0.3, 1, 0.025, 0.8, 'Chisq'), 'positive integers')
  expect_error(BinaryBSSR(20, 20, 5, 0, 1, 0.025, 0.8, 'Chisq'), 'Delta.A')
  expect_error(BinaryBSSR(20, 20, 5, 0.3, 1, 0.025, 0.8, 'Chisq', restricted = TRUE),
               'must be supplied')
})

test_that("the final sizes honour the allocation ratio for a whole r", {
  for (r in c(1, 2, 3)) {
    res <- BinaryBSSR(n1 = ceiling(r * 12), n2 = 12, S = 9, Delta.A = 0.3, r = r,
                      alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
    expect_equal(res$N1.final, as.integer(ceiling(r * res$N2.final)),
                 info = sprintf('r = %g', r))
  }
})

test_that("an imbalance in the observed interim data is corrected, not carried forward", {
  # Group 1 is two patients short of the two to one target at the interim
  res <- BinaryBSSR(n1 = 18, n2 = 10, S = 9, Delta.A = 0.3, r = 2,
                    alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
  expect_equal(res$N1.final, as.integer(2 * res$N2.final))
  expect_gte(res$n1.stage2, 0)
})

test_that("a fractional allocation ratio is handled", {
  res <- BinaryBSSR(n1 = 18, n2 = 12, S = 9, Delta.A = 0.3, r = 1.5,
                    alpha = 0.025, tar.power = 0.8, Test = 'Chisq')
  expect_equal(res$N1.final, as.integer(ceiling(1.5 * res$N2.final)))
  expect_lt(abs(res$N1.final - 1.5 * res$N2.final), 1)
})

test_that("the lower-tail alternative mirrors the upper-tail one", {
  up <- BinaryBSSR(20, 20, 11, 0.3, 1, 0.025, 0.8, 'Boschloo')
  down <- BinaryBSSR(20, 20, 11, -0.3, 1, 0.025, 0.8, 'Boschloo', alternative = 'less')
  expect_equal(c(down$hat.p1, down$hat.p2), c(up$hat.p2, up$hat.p1))
  expect_equal(c(down$N1.final, down$N2.final), c(up$N1.final, up$N2.final))
  expect_equal(down$Power, up$Power, tolerance = 1e-12)
})

test_that("the normal methods, the rounding rules and the bounds are applied", {
  res <- BinaryBSSR(20, 20, 11, 0.3, 1, 0.025, 0.8, 'Chisq', ss.method = 'standard')
  n2 <- ss_raw_n2(res$hat.p1, res$hat.p2, 1, 0.025, 0.8, 'greater', 'standard', 0, 'RD')
  expect_equal(res$N2.re, ceil_tol(n2))
  capped <- BinaryBSSR(20, 20, 11, 0.3, 1, 0.025, 0.8, 'Chisq', ss.method = 'standard',
                       N.max = 50)
  expect_lte(capped$N.final, 50)
  floor.n <- BinaryBSSR(20, 20, 11, 0.3, 1, 0.025, 0.8, 'Chisq', ss.method = 'standard',
                        N.min = 120)
  expect_gte(floor.n$N.final, 120)
  tot <- BinaryBSSR(20, 20, 11, 0.3, 1, 0.025, 0.8, 'Chisq', ss.method = 'standard',
                    rounding = 'total')
  expect_equal(tot$N.final, max(40, ceil_tol(2 * n2)))
})

test_that("coinciding recovered rates require the planned sample size", {
  expect_error(BinaryBSSR(10, 10, 0, 1.5, 1, 0.025, 0.8, 'Chisq', effect = 'RR'),
               'coincide')
  res <- BinaryBSSR(10, 10, 0, 1.5, 1, 0.025, 0.8, 'Chisq', effect = 'RR', N1 = 30,
                    N2 = 30)
  expect_equal(c(res$N1.final, res$N2.final), c(30L, 30L))
})

test_that("BinaryBSSR passes ref.pvalue to the re-estimation and the final analysis", {
  seen <- ref_pvalue_calls(BinaryBSSR(10, 10, 9, 0.3, 1, 0.025, 0.8, 'Boschloo',
                                      ref.pvalue = TRUE))
  expect_gt(nrow(seen), 1)
  expect_true(all(seen$ref.pvalue))
  res <- BinaryBSSR(10, 10, 9, 0.3, 1, 0.025, 0.8, 'Boschloo', ref.pvalue = TRUE)
  expect_true(attr(res, 'ref.pvalue'))
  seen <- ref_pvalue_calls(BinaryBSSR(10, 10, 9, 0.3, 1, 0.025, 0.8, 'Boschloo'))
  expect_gt(nrow(seen), 0)
  expect_false(any(seen$ref.pvalue))
})

test_that("BinaryBSSR re-estimates under a non-inferiority margin", {
  res <- BinaryBSSR(30, 30, 24, 0, 1, 0.025, 0.8, 'Farrington-Manning',
                    ss.method = 'standard', rounding = 'nearest', margin = 0.2)
  n2 <- ss_raw_n2(0.4, 0.4, 1, 0.025, 0.8, 'greater', 'standard', 0.2, 'RD')
  expect_equal(c(res$N1.final, res$N2.final), rep(max(30, floor(n2 + 0.5)), 2))
  expect_equal(attr(res, 'margin'), 0.2)
  expect_error(BinaryBSSR(30, 30, 24, -0.2, 1, 0.025, 0.8, 'Farrington-Manning',
                          ss.method = 'standard', margin = 0.2), 'exceed -margin')
  seen <- ref_pvalue_calls(BinaryBSSR(30, 30, 24, 0, 1, 0.025, 0.8, 'Farrington-Manning',
                                      margin = 0.2))
  expect_gt(nrow(seen), 1)
  expect_true(all(seen$margin == 0.2))
})

test_that("BinaryBSSR passes tsmethod to the re-estimation", {
  seen <- ref_pvalue_calls(res <- BinaryBSSR(8, 4, 5, 0.3, 2, 0.05, 0.8, 'Fisher',
                                             alternative = 'two.sided',
                                             tsmethod = 'blaker'))
  expect_gt(nrow(seen), 0)
  expect_true(all(seen$tsmethod == 'blaker'))
})

test_that("BinaryBSSR passes search and search.limit to the re-estimation", {
  seen <- search_calls(res <- BinaryBSSR(10, 10, 9, 0.3, 1, 0.025, 0.8, 'Chisq',
                                         search = 'stable', search.limit = c(3, 10)))
  expect_gt(nrow(seen), 0)
  expect_true(all(seen$search == 'stable' & seen$a == 3 & seen$b == 10))
  expect_equal(attr(res, 'search'), 'stable')
  # 9 responders among 20 give p1 = 0.6 and p2 = 0.3, for which the normal approximation
  # gives 42 patients in group 2 and the limit is 3 x 42. Reference value from
  # tools/reference/reference_values.py
  expect_equal(attr(res, 'search.limit'), 126)
  expect_true(is.na(attr(BinaryBSSR(10, 10, 9, 0.3, 1, 0.025, 0.8, 'Chisq'),
                         'search.limit')))
})

test_that("BinaryBSSR re-estimates under a margin on the scale of the risk ratio", {
  run <- function(S, ...) {
    BinaryBSSR(20, 20, S, 1, 1, 0.025, 0.8, 'Farrington-Manning', effect = 'RR',
               ss.method = 'standard', rounding = 'nearest', N.max = 300, margin = 0.8,
               margin.scale = 'RR', ...)
  }
  # 28 responders give the pooled rate 0.7 for both groups, and formula (8) of
  # Farrington and Manning (1990) gives 141 per group, see test-binary-power-bssr.R
  seen <- ref_pvalue_calls(res <- run(28))
  expect_true(all(seen$margin.scale == 'RR'))
  n2 <- ss_raw_n2(0.7, 0.7, 1, 0.025, 0.8, 'greater', 'standard', 0.8, 'RR')
  expect_equal(c(res$N1.final, res$N2.final), rep(floor(n2 + 0.5), 2))
  expect_equal(res$N2.final, 141)
  expect_identical(attr(res, 'margin.scale'), 'RR')
  # Without responders the recovered rates lie on the null boundary, so the planned
  # sample sizes are needed
  expect_error(run(0), 'null boundary')
  expect_equal(run(0, N1 = 100, N2 = 100)$N.final, 200)
})
