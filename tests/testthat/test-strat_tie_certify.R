test_that("strat_tie_certify reports the maximum of strat_bern_max", {
  coef <- outer(0:5, 0:4, function(i, k) (i / 5) * (1 - k / 4) + 0.1)
  box <- rbind(c(0.1, 0.2), c(0.7, 0.9))
  m <- strat_tie_certify(coef, box)
  # The polynomial x (1 - y) + 0.1 is largest at x = 0.7, y = 0.2
  expect_equal(m$x, c(0.7, 0.2))
  expect_equal(m$y, 0.7 * 0.8 + 0.1)
  expect_gte(m$bound, m$y)
})
