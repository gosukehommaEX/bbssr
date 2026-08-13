test_that("BinarySampleSize returns a bbssr_samplesize data frame", {
  res <- BinarySampleSize(p1 = 0.5, p2 = 0.2, r = 1, alpha = 0.025,
                          tar.power = 0.8, Test = 'Chisq')
  expect_s3_class(res, 'bbssr_samplesize')
  expect_equal(nrow(res), 1L)
  expect_named(res, c('p1', 'p2', 'r', 'alpha', 'tar.power', 'Test', 'alternative',
                      'Power', 'N1', 'N2', 'N'))
  expect_type(res$N1, 'integer')
  expect_type(res$N2, 'integer')
  expect_equal(res$N, res$N1 + res$N2)
})

test_that("the returned sample size is the smallest one attaining the target power", {
  for (tst in c('Chisq', 'Fisher')) {
    res <- BinarySampleSize(0.6, 0.25, 1, 0.025, 0.8, tst)
    at <- BinaryPower(0.6, 0.25, res$N1, res$N2, 0.025, tst)$Power
    below <- BinaryPower(0.6, 0.25, ceiling(res$N2 - 1), res$N2 - 1, 0.025, tst)$Power
    expect_gte(at, 0.8)
    expect_lt(below, 0.8)
  }
})

test_that("the allocation ratio is respected", {
  for (r in c(0.5, 1, 2)) {
    res <- BinarySampleSize(0.6, 0.3, r, 0.025, 0.8, 'Chisq')
    expect_equal(res$N1, as.integer(ceiling(r * res$N2)))
  }
})

test_that("a larger target power requires a larger sample size", {
  n <- vapply(c(0.7, 0.8, 0.9), function(tp) {
    BinarySampleSize(0.5, 0.25, 1, 0.025, tp, 'Chisq')$N
  }, numeric(1))
  expect_true(all(diff(n) > 0))
})

test_that("a smaller treatment effect requires a larger sample size", {
  n <- vapply(c(0.7, 0.6, 0.5), function(p1) {
    BinarySampleSize(p1, 0.3, 1, 0.025, 0.8, 'Chisq')$N
  }, numeric(1))
  expect_true(all(diff(n) > 0))
})

test_that("the two-sided test requires at least as many patients as the one-sided test", {
  one <- BinarySampleSize(0.6, 0.3, 1, 0.025, 0.8, 'Fisher')$N
  two <- BinarySampleSize(0.6, 0.3, 1, 0.025, 0.8, 'Fisher',
                          alternative = 'two.sided')$N
  expect_gte(two, one)
})

test_that("the Boschloo test needs no more patients than the Fisher test", {
  fisher <- BinarySampleSize(0.6, 0.2, 1, 0.025, 0.8, 'Fisher')$N
  boschloo <- BinarySampleSize(0.6, 0.2, 1, 0.025, 0.8, 'Boschloo')$N
  expect_lte(boschloo, fisher)
})

test_that("BinarySampleSize validates its arguments", {
  expect_error(BinarySampleSize(0.4, 0.4, 1, 0.025, 0.8, 'Chisq'), 'must differ')
  expect_error(BinarySampleSize(0.5, 0.2, 0, 0.025, 0.8, 'Chisq'), 'positive')
  expect_error(BinarySampleSize(0.5, 0.2, 1, 0.025, 1, 'Chisq'), 'tar.power')
  expect_error(BinarySampleSize(c(0.5, 0.6), 0.2, 1, 0.025, 0.8, 'Chisq'), 'single value')
})
