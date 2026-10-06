test_that("maximize_label describes how the largest type I error rate was located", {
  x <- structure(list(), maximize = 'certified', interval = c(0.1, 0.9))
  expect_equal(maximize_label(x), 'certified over theta in [0.1, 0.9]')
  expect_equal(maximize_label(structure(list(), maximize = 'refined')),
               'largest values on the grid refined')
  expect_equal(maximize_label(structure(list(), maximize = 'grid')),
               'largest value on the grid')
  expect_equal(maximize_label(list()), 'not recorded')
})
