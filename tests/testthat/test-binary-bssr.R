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
