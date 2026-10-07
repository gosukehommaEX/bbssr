test_that("the exact method reproduces BinarySampleSize", {
  for (tst in c('Chisq', 'Fisher', 'Z-pool')) {
    n <- sample_size_n(0.6, 0.25, 1, 0.025, 0.8, tst, 'greater', 'minlike', 100L, 0,
                       'exact', 'group', FALSE, 0)
    ss <- BinarySampleSize(0.6, 0.25, 1, 0.025, 0.8, tst)
    expect_identical(as.vector(n), c(ss$N1, ss$N2), info = tst)
    expect_true(is.na(attr(n, 'limit')))
  }
})

test_that("sample_size_n passes the search and returns the limit of the stable search", {
  # Reference values from tools/reference/reference_values.py, see test-ss_exact_search.R
  n <- sample_size_n(0.6, 0.3, 1, 0.025, 0.85, 'Fisher', 'greater', 'minlike', 100L, 0,
                     'exact', 'group', FALSE, 0, 'stable')
  expect_equal(as.vector(n), c(56, 56))
  expect_equal(attr(n, 'limit'), 98)
})

test_that("the three rounding rules follow their definitions", {
  # Unrounded size of group 2 is 23.5466 for r = 3, so the total is 94.19
  args <- list(p1 = 0.95, p2 = 0.75, r = 3, alpha = 0.05, tar.power = 0.8, Test = 'Chisq',
               alternative = 'two.sided', tsmethod = 'minlike', n.grid = 100L,
               bb.gamma = 0, method = 'null.variance')
  rest <- list(ref.pvalue = FALSE, margin = 0)
  grp <- do.call(sample_size_n, c(args, list(rounding = 'group'), rest))
  fk <- do.call(sample_size_n, c(args, list(rounding = 'friede-kieser'), rest))
  tot <- do.call(sample_size_n, c(args, list(rounding = 'total'), rest))
  expect_identical(grp, c(N1 = 72L, N2 = 24L))
  # Friede and Kieser (2004), Table I: 95 patients for theta = 3 and p1 = 0.75
  expect_identical(fk, c(N1 = 71L, N2 = 24L))
  expect_identical(tot, c(N1 = 72L, N2 = 23L))
})

test_that("sample_size_n rounds each group to the nearest whole number under 'nearest'", {
  n <- sample_size_n(0.2, 0.1, 1, 0.05, 0.9, 'Farrington-Manning', 'greater', 'minlike',
                     100L, 0, 'standard', 'nearest', FALSE, 0.1)
  n2 <- ss_raw_n2(0.2, 0.1, 1, 0.05, 0.9, 'greater', 'standard', 0.1)
  expect_identical(unname(n), as.integer(c(floor(n2 + 0.5), floor(n2 + 0.5))))
  # 57 per group in Table I of Farrington and Manning (1990)
  expect_identical(unname(n), c(57L, 57L))
})
