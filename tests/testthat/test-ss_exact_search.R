test_that("ss_exact_search reproduces separate searches for every pair", {
  p1 <- c(0.5, 0.6, 0.45, 0.7)
  p2 <- c(0.2, 0.3, 0.25, 0.3)
  for (tst in c('Chisq', 'Fisher', 'Boschloo')) {
    got <- ss_exact_search(p1, p2, 1, 0.025, 0.8, tst, 'greater', 'minlike', 100L, 0,
                           FALSE)
    expect_type(got, 'integer')
    one <- vapply(seq_along(p1), function(k) {
      ss_exact_search(p1[k], p2[k], 1, 0.025, 0.8, tst, 'greater', 'minlike', 100L, 0,
                      FALSE)
    }, integer(1))
    expect_identical(got, one, info = tst)
  }
})

test_that("the size attains the target power and one patient fewer does not", {
  N2 <- ss_exact_search(0.6, 0.25, 2, 0.025, 0.8, 'Chisq', 'greater', 'minlike', 100L, 0,
                        FALSE)
  at <- BinaryPower(0.6, 0.25, ceiling(2 * N2), N2, 0.025, 'Chisq')$Power
  below <- BinaryPower(0.6, 0.25, ceiling(2 * (N2 - 1)), N2 - 1, 0.025, 'Chisq')$Power
  expect_gte(at, 0.8)
  expect_lt(below, 0.8)
  expect_identical(ss_exact_search(numeric(0), numeric(0), 1, 0.025, 0.8, 'Chisq',
                                   'greater', 'minlike', 100L, 0, FALSE), integer(0))
})
