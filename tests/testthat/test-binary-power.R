all_tests <- c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')

test_that("BinaryPower returns a bbssr_power data frame", {
  res <- BinaryPower(p1 = 0.5, p2 = 0.2, N1 = 8, N2 = 8, alpha = 0.025, Test = 'Chisq')
  expect_s3_class(res, 'bbssr_power')
  expect_s3_class(res, 'data.frame')
  expect_equal(nrow(res), 1L)
  expect_named(res, c('p1', 'p2', 'N1', 'N2', 'alpha', 'Test', 'alternative', 'Power'))
  expect_type(res$Power, 'double')
  expect_true(res$Power >= 0 && res$Power <= 1)
})

test_that("BinaryPower is vectorised over the response probabilities", {
  res <- BinaryPower(p1 = c(0.4, 0.6, 0.8), p2 = rep(0.2, 3),
                     N1 = 10, N2 = 10, alpha = 0.025, Test = 'Chisq')
  expect_equal(nrow(res), 3L)
  expect_true(all(res$Power >= 0 & res$Power <= 1))
  expect_true(all(diff(res$Power) > 0))
})

test_that("BinaryPower agrees with a direct summation over the rejection region", {
  N1 <- 9
  N2 <- 6
  p1 <- 0.7
  p2 <- 0.25
  for (tst in all_tests) {
    RR <- BinaryRR(N1, N2, 0.025, tst, n.grid = 30)
    manual <- sum(outer(stats::dbinom(0:N1, N1, p1), stats::dbinom(0:N2, N2, p2)) *
                    as_plain(RR))
    got <- BinaryPower(p1, p2, N1, N2, 0.025, tst, n.grid = 30)$Power
    expect_equal(got, manual, tolerance = 1e-12, info = tst)
  }
})

test_that("BinaryPower rises with the sample size apart from the discreteness sawtooth", {
  # The exact power is not monotone in the sample size, because the attainable
  # significance level changes with each increment. At N1 = N2 = 8 the power of the
  # chi-squared test is 0.447 and at N1 = N2 = 10 it drops to 0.425, so only the overall
  # trend and a bound on the local decrease are asserted
  powers <- vapply(seq(6, 36, by = 6), function(n) {
    BinaryPower(0.6, 0.2, n, n, 0.025, 'Chisq')$Power
  }, numeric(1))
  expect_gt(powers[length(powers)], powers[1])
  expect_true(all(diff(powers) > -0.05))
})

test_that("BinaryPower is monotone in the effect size", {
  powers <- BinaryPower(seq(0.3, 0.8, by = 0.1), rep(0.2, 6), 15, 15, 0.025, 'Fisher')$Power
  expect_true(all(diff(powers) >= 0))
})

test_that("power under the null does not exceed the significance level", {
  alpha <- 0.05
  for (tst in c('Fisher', 'Z-pool', 'Boschloo')) {
    for (alt in c('greater', 'two.sided')) {
      p <- BinaryPower(0.4, 0.4, 10, 10, alpha, tst,
                       alternative = alt, n.grid = 100)$Power
      expect_lte(p, alpha + 1e-10)
    }
  }
})

test_that("the two-sided test is less powerful than the one-sided test", {
  one <- BinaryPower(0.7, 0.3, 15, 15, 0.025, 'Fisher')$Power
  two <- BinaryPower(0.7, 0.3, 15, 15, 0.025, 'Fisher', alternative = 'two.sided')$Power
  expect_lt(two, one)
})

test_that("BinaryPower validates its arguments", {
  expect_error(BinaryPower(c(0.4, 0.5), 0.2, 10, 10, 0.025, 'Chisq'), 'same length')
  expect_error(BinaryPower(-0.1, 0.2, 10, 10, 0.025, 'Chisq'), 'lie in')
  expect_error(BinaryPower(0.4, 1.2, 10, 10, 0.025, 'Chisq'), 'lie in')
})

test_that("BinaryPower passes ref.pvalue to the p-value computation", {
  seen <- ref_pvalue_calls(BinaryPower(0.6, 0.2, 10, 10, 0.025, 'Z-pool',
                                       ref.pvalue = TRUE))
  expect_equal(seen$ref.pvalue, TRUE)
  # The refined region at 32 patients per group drops the outcomes 18 versus 10 and
  # 22 versus 14, so the power at a common response probability decreases
  grid <- BinaryPower(0.45, 0.45, 32, 32, 0.025, 'Z-pool')$Power
  ref <- BinaryPower(0.45, 0.45, 32, 32, 0.025, 'Z-pool', ref.pvalue = TRUE)$Power
  expect_lt(ref, grid)
})

test_that("BinaryPower distinguishes the three two-sided conventions", {
  # A configuration of Table 3 of Mehrotra, Chan and Berger (2003), with values from
  # tools/reference/reference_values.py
  pw <- vapply(c('blaker', 'minlike', 'central'), function(ts) {
    BinaryPower(0.5, 0.86, 10, 40, 0.05, 'Fisher', alternative = 'two.sided',
                tsmethod = ts)$Power
  }, numeric(1))
  expect_equal(unname(pw), c(0.608054130818687, 0.641394819757925, 0.517056474695217),
               tolerance = 1e-10)
  seen <- ref_pvalue_calls(BinaryPower(0.5, 0.86, 10, 40, 0.05, 'Boschloo',
                                       alternative = 'two.sided', tsmethod = 'blaker'))
  expect_equal(seen$tsmethod, 'blaker')
})

test_that("the power of the Farrington-Manning test reproduces the article", {
  # Example of Farrington and Manning (1990): p1 = 0.4, p2 = 0.05, s0 = 0.2 and 80 patients
  # per group give a true power of 81.3 per cent (p. 1451)
  pw <- BinaryPower(0.4, 0.05, 80, 80, 0.05, 'Farrington-Manning', margin = -0.2)$Power
  expect_equal(round(100 * pw, 1), 81.3)
  expect_equal(pw, 0.81320090915784, tolerance = 1e-8)
  # Table I: p1 = 0.2, p2 = 0.1, s0 = -0.1 with 57 per group (91.34), p1 = 0.5, p2 = 0.1,
  # s0 = 0.2 with 67 and 101 (90.23), and p1 = 0.1, p2 = 0.05, s0 = -0.05 with 168 and 112
  # (90.84). The exact values come from tools/reference/reference_values.py
  fm <- function(p1, p2, N1, N2, margin) {
    BinaryPower(p1, p2, N1, N2, 0.05, 'Farrington-Manning', margin = margin)$Power
  }
  pw <- c(fm(0.2, 0.1, 57, 57, 0.1), fm(0.5, 0.1, 67, 101, -0.2),
          fm(0.1, 0.05, 168, 112, 0.05))
  expect_equal(round(100 * pw, 1), round(c(91.34, 90.23, 90.84), 1))
  expect_equal(pw, c(0.913399732016454, 0.902410125176024, 0.908496094749286),
               tolerance = 1e-8)
})

test_that("BinaryPower passes the margin to the p-value computation", {
  seen <- ref_pvalue_calls(pw <- BinaryPower(0.65, 0.7, 120, 100, 0.025, 'Blackwelder',
                                             margin = 0.15))
  expect_equal(seen$margin, 0.15)
  expect_equal(pw$Power, 0.356600952518614, tolerance = 1e-8)
  expect_equal(attr(pw, 'margin'), 0.15)
})

test_that("BinaryPower computes the power for a margin on the scale of the risk ratio", {
  # Reference values from tools/reference/reference_values.py: the Farrington-Manning
  # test of p1 / p2 <= 0.8 with 150 per group and the Blackwelder test of
  # p1 / p2 >= 1.5 with 120 and 100 patients, both at equal true rates
  seen <- ref_pvalue_calls(
    pw <- c(BinaryPower(0.6, 0.6, 150, 150, 0.025, 'Farrington-Manning', margin = 0.8,
                        margin.scale = 'RR')$Power,
            BinaryPower(0.3, 0.3, 120, 100, 0.025, 'Blackwelder', alternative = 'less',
                        margin = 1.5, margin.scale = 'RR')$Power)
  )
  expect_identical(seen$margin.scale, c('RR', 'RR'))
  expect_equal(pw, c(0.652502871057245, 0.457689359850127), tolerance = 1e-8)
  res <- BinaryPower(0.6, 0.6, 15, 15, 0.025, 'Farrington-Manning', margin = 0.8,
                     margin.scale = 'RR')
  expect_identical(attr(res, 'margin.scale'), 'RR')
})
