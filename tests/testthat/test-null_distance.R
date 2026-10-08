test_that("null_distance is the distance from the boundary on the side of the alternative", {
  p1 <- c(0.5, 0.3)
  p2 <- c(0.4, 0.45)
  expect_equal(null_distance(p1, p2, 'greater', 0.1, 'RD'), c(0.2, -0.05))
  expect_equal(null_distance(p1, p2, 'less', 0.1, 'RD'), c(0, 0.25))
  expect_equal(null_distance(p1, p2, 'greater', 0.8, 'RR'), c(0.18, -0.06))
  expect_equal(null_distance(p1, p2, 'less', 1.25, 'RR'), c(0, 0.2625))
  # Without a margin the distance is the difference p1 - p2
  expect_identical(null_distance(p1, p2, 'two.sided', 0, 'RD'), p1 - p2)
})
