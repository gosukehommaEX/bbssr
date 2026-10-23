test_that("paired_boundary gives the probabilities on the null boundary", {
  b <- paired_boundary(c(0, 0.3, 1), 'greater', 0)
  expect_equal(b$p10, c(0, 0.15, 0.5))
  expect_equal(b$p01, c(0, 0.15, 0.5))
  expect_true(all(b$ok))
  expect_equal(b$t, c(0, 0.3, 1))
  expect_equal(c(b$type, b$m), c(0, 0))
  expect_equal(b$range, c(0, 1))
  # Non-inferiority, p10 - p01 = -0.2
  b <- paired_boundary(c(0.1, 0.2, 0.6, 1), 'greater', 0.2)
  expect_equal(b$ok, c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(b$p10[-1] - b$p01[-1], rep(-0.2, 3))
  expect_equal(b$t[-1], c(0, 0.5, 1))
  expect_equal(c(b$type, b$m), c(1, 0.2))
  expect_equal(b$range, c(0.2, 1))
  # The alternative 'less' with margin 0.2 gives p10 - p01 = 0.2
  b <- paired_boundary(c(0.4, 0.9), 'less', 0.2)
  expect_equal(b$p10 - b$p01, c(0.2, 0.2))
  expect_equal(c(b$type, b$m), c(2, 0.2))
  expect_equal(b$range, c(0.4, 0.9))
  # A negative margin, superiority by 0.1: p10 - p01 = 0.1
  b <- paired_boundary(0.5, 'greater', -0.1)
  expect_equal(c(b$p10, b$p01, b$type), c(0.3, 0.2, 2))
})

test_that("paired_boundary has the linear form used by the Bernstein coefficients", {
  for (case in list(c('greater', 0.15), c('less', 0.15))) {
    b <- paired_boundary(seq(0.15, 1, by = 0.05), case[1], as.numeric(case[2]))
    m <- b$m
    p.u <- if (b$type == 1) b$p10 else b$p01
    expect_equal(p.u, (1 - m) * b$t / 2)
    expect_equal(1 - b$p10 - b$p01, (1 - m) * (1 - b$t))
  }
})

test_that("paired_boundary reports no range without values on the boundary", {
  b <- paired_boundary(c(0.05, 0.1), 'greater', 0.2)
  expect_false(any(b$ok))
  expect_true(all(is.na(b$range)))
})
