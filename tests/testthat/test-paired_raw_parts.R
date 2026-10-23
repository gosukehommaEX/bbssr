test_that("paired_raw_parts reproduces the example of Nam (1997, Section 5)", {
  # p01 = .05, delta0 = -.05, delta1 = 0: p01.bar = .0845, so v0 = 2 p01.bar - .05 - .0025
  s <- paired_raw_parts(0.05, 0.05, 'greater', 0.05)
  expect_equal(round((s$v0 + 0.05 + 0.05^2) / 2, 4), 0.0845)
  expect_equal(s$d, 0.05)
  expect_equal(s$v1, 0.1)
  # Reference value of the limit, from a root search of the expected score equation
  expect_equal((s$v0 + 0.05 + 0.05^2) / 2, 0.08446012188022556, tolerance = 1e-12)
})

test_that("paired_raw_parts gives the variances of Connor (1987) for a zero margin", {
  s <- paired_raw_parts(c(0.45, 0.2), c(0.05, 0.1), 'greater', 0)
  expect_equal(s$v0, c(0.5, 0.3))
  expect_equal(s$v1, c(0.5 - 0.16, 0.3 - 0.01))
  expect_equal(s$d, c(0.4, 0.1))
})

test_that("paired_raw_parts exchanges the probabilities for the alternative 'less'", {
  expect_equal(paired_raw_parts(0.1, 0.25, 'less', 0.05),
               paired_raw_parts(0.25, 0.1, 'greater', 0.05))
})

test_that("paired_raw_parts handles probabilities without discordant pairs", {
  # p10 = p01 = 0 and margin 0.2: p01.bar = 0.2, so v0 = 0.2 - 0.04
  s <- paired_raw_parts(0, 0, 'greater', 0.2)
  expect_equal(c(s$d, s$v0, s$v1), c(0.2, 0.16, 0))
})
