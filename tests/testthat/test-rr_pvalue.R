all_tests <- c('Chisq', 'Fisher', 'Fisher-midP', 'Z-pool', 'Boschloo')

test_that("the lower-tail p-values are those of the exchanged groups", {
  for (tst in all_tests) {
    less <- rr_pvalue(6L, 4L, tst, 'less', 'minlike', 30L, 0, FALSE, 0, 'RD')
    swap <- rr_pvalue(4L, 6L, tst, 'greater', 'minlike', 30L, 0, FALSE, 0, 'RD')
    expect_equal(dim(less), c(7L, 5L))
    expect_identical(less, t(swap), info = tst)
  }
})

test_that("the lower-tail Fisher and chi-squared p-values follow their definitions", {
  N1 <- 6
  N2 <- 5
  fisher <- rr_pvalue(N1, N2, 'Fisher', 'less', 'minlike', 100L, 0, FALSE, 0, 'RD')
  manual <- outer(0:N1, 0:N2, function(i, j) stats::phyper(i, N1, N2, i + j))
  expect_equal(fisher, manual, tolerance = 1e-12)
  chisq <- rr_pvalue(N1, N2, 'Chisq', 'less', 'minlike', 100L, 0, FALSE, 0, 'RD')
  expect_equal(chisq, stats::pnorm(zstat(N1, N2)), tolerance = 1e-12)
})

test_that("the lower-tail p-value agrees with stats::fisher.test", {
  N1 <- 5
  N2 <- 7
  p <- rr_pvalue(N1, N2, 'Fisher', 'less', 'minlike', 100L, 0, FALSE, 0, 'RD')
  for (i in 0:N1) {
    for (j in 0:N2) {
      tab <- matrix(c(i, j, N1 - i, N2 - j), nrow = 2)
      expect_equal(p[i + 1, j + 1], stats::fisher.test(tab, alternative = 'less')$p.value,
                   tolerance = 1e-12)
    }
  }
})

test_that("the upper-tail and two-sided p-values are those of version 2.0.0", {
  N1 <- 7
  N2 <- 6
  expect_equal(rr_pvalue(N1, N2, 'Fisher', 'greater', 'minlike', 100L, 0, FALSE, 0, 'RD'),
               fisher_greater_ref(N1, N2))
  stat <- zstat(N1, N2)
  expect_equal(rr_pvalue(N1, N2, 'Z-pool', 'greater', 'minlike', 25L, 0, FALSE, 0, 'RD'),
               unconditional_ref(stat, N1, N2, 25L, decreasing = TRUE), tolerance = 1e-12)
  expect_equal(rr_pvalue(N1, N2, 'Chisq', 'two.sided', 'minlike', 100L, 0, FALSE, 0,
                         'RD'),
               pmin(2 * stats::pnorm(abs(stat), lower.tail = FALSE), 1))
})

test_that("the lower tail of a ratio margin inverts the margin of the exchanged groups", {
  for (tst in c('Blackwelder', 'Farrington-Manning')) {
    se <- if (tst == 'Blackwelder') 'unpooled' else 'restricted'
    less <- rr_pvalue(9L, 12L, tst, 'less', 'minlike', 100L, 0, FALSE, 1.25, 'RR')
    swap <- rr_pvalue(12L, 9L, tst, 'greater', 'minlike', 100L, 0, FALSE, 0.8, 'RR')
    expect_identical(less, t(swap), info = tst)
    # The exchange only changes the sign of the statistic with the ratio 1.25
    z <- zstat_margin(9L, 12L, 1.25, se, 'RR')
    expect_equal(less, stats::pnorm(z), tolerance = 1e-12, info = tst)
    greater <- rr_pvalue(9L, 12L, tst, 'greater', 'minlike', 100L, 0, FALSE, 1.25, 'RR')
    expect_equal(greater, stats::pnorm(z, lower.tail = FALSE), tolerance = 1e-12,
                 info = tst)
  }
})
