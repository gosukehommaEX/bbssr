test_that("fm_restricted reproduces the example of Farrington and Manning (1990)", {
  # p1 = 0.4, p2 = 0.05 and s0 = 0.2 with equal groups give the large sample values
  # 0.2935 and 0.0935 (p. 1451)
  est <- fm_restricted(0.4, 0.05, 1, 0.2)
  expect_equal(round(c(est$p1, est$p2), 4), c(0.2935, 0.0935))
})

test_that("fm_restricted agrees with a direct maximization of the likelihood", {
  # Reference values from tools/reference/reference_values.py, which maximizes the
  # restricted likelihood by bisection on the score. Each row is x1, n1, x2, n2 and s0
  cases <- rbind(c(3, 12, 7, 20, -0.15), c(0, 18, 0, 18, -0.2), c(12, 12, 20, 20, 0.1),
                 c(5, 10, 0, 15, 0.3), c(8, 25, 9, 30, -0.05))
  got <- unlist(lapply(seq_len(nrow(cases)), function(i) {
    x <- cases[i, ]
    est <- fm_restricted(x[1] / x[2], x[3] / x[4], x[4] / x[2], x[5])
    c(est$p1, est$p2)
  }))
  want <- c(0.222389288048392, 0.372389288048392, 0, 0.2, 1, 0.9, 0.3, 0,
            0.283388372295963, 0.333388372295963)
  expect_equal(got, want, tolerance = 1e-9)
  # Large sample values for p1 = p2 = 0.7, theta = 1 / 3 and s0 = -0.1
  est <- fm_restricted(0.7, 0.7, 1 / 3, -0.1)
  expect_equal(c(est$p1, est$p2), c(0.670595327339649, 0.770595327339649),
               tolerance = 1e-9)
})

test_that("fm_restricted gives the pooled proportion when s0 is 0", {
  est <- fm_restricted(c(0.2, 0.5), c(0.4, 0.1), 2, 0)
  pooled <- (c(0.2, 0.5) + 2 * c(0.4, 0.1)) / 3
  expect_equal(est$p1, pooled, tolerance = 1e-12)
  expect_equal(est$p2, pooled, tolerance = 1e-12)
})

test_that("fm_restricted stays in the admissible range over a whole outcome grid", {
  p1 <- rep((0:12) / 12, times = 16)
  p2 <- rep((0:15) / 15, each = 13)
  for (s0 in c(-0.3, -0.05, 0.1, 0.4)) {
    est <- fm_restricted(p1, p2, 15 / 12, s0)
    expect_true(all(est$p1 >= max(0, s0) & est$p1 <= min(1, 1 + s0)), info = s0)
    expect_equal(est$p1 - est$p2, rep(s0, length(p1)), tolerance = 1e-12)
  }
})
