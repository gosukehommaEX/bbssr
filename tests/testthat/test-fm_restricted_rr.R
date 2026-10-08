test_that("fm_restricted_rr reproduces the second example of Farrington and Manning (1990)", {
  # p1 = p2 = 0.01, R0 = 1.5 and equal groups give the large sample values 0.012 and
  # 0.008 (p. 1453)
  est <- fm_restricted_rr(0.01, 0.01, 1, 1.5)
  expect_equal(round(c(est$p1, est$p2), 3), c(0.012, 0.008))
})

test_that("fm_restricted_rr agrees with a direct maximization of the likelihood", {
  # Reference values from tools/reference/reference_values.py, which maximizes the
  # restricted likelihood by bisection on the score. Each row is x1, n1, x2, n2 and R0
  cases <- rbind(c(3, 12, 7, 20, 0.5), c(0, 18, 0, 18, 0.8), c(12, 12, 20, 20, 0.9),
                 c(5, 10, 0, 15, 2), c(8, 25, 9, 30, 1.5), c(10, 10, 5, 10, 2))
  got <- unlist(lapply(seq_len(nrow(cases)), function(i) {
    x <- cases[i, ]
    est <- fm_restricted_rr(x[1] / x[2], x[3] / x[4], x[4] / x[2], x[5])
    c(est$p1, est$p2)
  }))
  want <- c(0.189028034951513, 0.378056069903026, 0, 0, 0.9, 1, 0.310102051443364,
            0.155051025721682, 0.372064896664658, 0.248043264443105, 1, 0.5)
  expect_equal(got, want, tolerance = 1e-9)
  # Large sample values for p1 = 0.65, p2 = 0.7, theta = 1 / 2 and R0 = 0.8
  est <- fm_restricted_rr(0.65, 0.7, 0.5, 0.8)
  expect_equal(c(est$p1, est$p2), c(0.604729180758055, 0.755911475947568),
               tolerance = 1e-9)
})

test_that("fm_restricted_rr gives the pooled proportion when R0 is 1", {
  est <- fm_restricted_rr(c(0.2, 0.5), c(0.4, 0.1), 2, 1)
  pooled <- (c(0.2, 0.5) + 2 * c(0.4, 0.1)) / 3
  expect_equal(est$p1, pooled, tolerance = 1e-12)
  expect_equal(est$p2, pooled, tolerance = 1e-12)
})

test_that("fm_restricted_rr stays in the admissible range over a whole outcome grid", {
  p1 <- rep((0:12) / 12, times = 16)
  p2 <- rep((0:15) / 15, each = 13)
  for (R0 in c(0.3, 0.8, 1.25, 4)) {
    est <- fm_restricted_rr(p1, p2, 15 / 12, R0)
    expect_true(all(est$p1 >= 0 & est$p1 <= min(1, R0)), info = R0)
    expect_true(all(est$p2 >= 0 & est$p2 <= 1), info = R0)
    expect_equal(est$p1, R0 * est$p2, tolerance = 1e-12)
  }
  # Proportions already on the restriction are returned unchanged
  est <- fm_restricted_rr(c(0.3, 0.4), c(0.6, 0.8), 1.5, 0.5)
  expect_equal(est$p1, c(0.3, 0.4), tolerance = 1e-12)
  expect_equal(est$p2, c(0.6, 0.8), tolerance = 1e-12)
})
