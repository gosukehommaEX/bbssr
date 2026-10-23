test_that("strat_n1 derives the interim sizes from omega or takes n.interim", {
  expect_identical(strat_n1(c(30, 21), c(1, 2), 0.5, NULL),
                   list(n1 = c(15L, 22L), n2 = c(15L, 11L)))
  expect_identical(strat_n1(c(30, 21), c(1, 2), NULL, rbind(c(10, 8), c(9, 7))),
                   list(n1 = c(10L, 8L), n2 = c(9L, 7L)))
})

test_that("strat_n1 rejects invalid interim sizes", {
  expect_error(strat_n1(c(30, 21), 1, 0.5, rbind(c(10, 8), c(9, 7))), 'either')
  expect_error(strat_n1(c(30, 21), 1, NULL, NULL), 'either')
  expect_error(strat_n1(c(30, 21), 1, 1.5, NULL), 'omega')
  expect_error(strat_n1(c(30, 21), 1, NULL, c(10, 8)), 'n.interim')
  expect_error(strat_n1(c(30, 21), 1, NULL, rbind(c(10, 0), c(9, 7))), 'n.interim')
})

test_that("strat_n1 is not raised by a rounding error of the products", {
  # 0.28 * 100 and 0.28 * 50 exceed 28 and 14 by a rounding error, and so does 1.1 * 50
  expect_identical(strat_n1(c(100, 50), c(1, 1), 0.28, NULL),
                   list(n1 = c(28L, 14L), n2 = c(28L, 14L)))
  expect_identical(strat_n1(c(100, 100), c(1, 1.1), 0.5, NULL),
                   list(n1 = c(50L, 55L), n2 = c(50L, 50L)))
})
