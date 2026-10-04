test_that("the exact method reproduces BinarySampleSize", {
  for (tst in c('Chisq', 'Fisher', 'Z-pool')) {
    n <- sample_size_n(0.6, 0.25, 1, 0.025, 0.8, tst, 'greater', 'minlike', 100L, 0,
                       'exact', 'group')
    ss <- BinarySampleSize(0.6, 0.25, 1, 0.025, 0.8, tst)
    expect_identical(unname(n), c(ss$N1, ss$N2), info = tst)
  }
})

test_that("the three rounding rules follow their definitions", {
  # Unrounded size of group 2 is 23.5466 for r = 3, so the total is 94.19
  args <- list(0.95, 0.75, 3, 0.05, 0.8, 'Chisq', 'two.sided', 'minlike', 100L, 0,
               'null.variance')
  grp <- do.call(sample_size_n, c(args, 'group'))
  fk <- do.call(sample_size_n, c(args, 'friede-kieser'))
  tot <- do.call(sample_size_n, c(args, 'total'))
  expect_identical(grp, c(N1 = 72L, N2 = 24L))
  # Friede and Kieser (2004), Table I: 95 patients for theta = 3 and p1 = 0.75
  expect_identical(fk, c(N1 = 71L, N2 = 24L))
  expect_identical(tot, c(N1 = 72L, N2 = 23L))
})
