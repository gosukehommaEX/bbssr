test_that("paired_score_z reproduces the examples of Nam (1997) and Tango (1998)", {
  # Nam (1997, Section 5): liver scans, n = 344, x10 = 32, x01 = 27, delta0 = -0.05,
  # z = 2.795 and one-sided p = .0026
  z <- paired_score_z(344, 0.05)[33, 28]
  expect_equal(round(z, 3), 2.795)
  expect_equal(round(stats::pnorm(z, lower.tail = FALSE), 4), 0.0026)
  # Tango (1998, Example 1): n = 44, b = 0, c = 1, Delta = 0.1, Z = 1.709, p = 0.044
  z <- paired_score_z(44, 0.1)[1, 2]
  expect_equal(round(z, 3), 1.709)
  expect_equal(round(stats::pnorm(z, lower.tail = FALSE), 3), 0.044)
  # Tango (1998, Example 2): n = 32, b = 9, c = 3, Delta = 0, Z = 1.73, two-sided p = 0.08
  z <- paired_score_z(32, 0)[10, 4]
  expect_equal(round(z, 2), 1.73)
  expect_equal(round(2 * stats::pnorm(z, lower.tail = FALSE), 2), 0.08)
})

test_that("paired_score_z reproduces Table I of Tango (1998)", {
  # b = 0 and c = 0, 1, 2 with Delta = 0.1, where the estimate lies on the boundary
  tab <- rbind(c(1.83, 1.22, 0.61), c(2.36, 1.89, 1.41), c(2.98, 2.61, 2.24))
  n <- c(30, 50, 80)
  for (i in 1:3) {
    expect_equal(round(paired_score_z(n[i], 0.1)[1, 1:3], 2), tab[i, ], info = n[i])
  }
})

test_that("paired_score_z agrees with the restricted maximum likelihood estimate", {
  # Reference values with the estimate found by a root search of the score equation
  Z <- paired_score_z(20, 0.1)
  expect_equal(c(Z[1, 1], Z[1, 4], Z[5, 1], Z[8, 3], Z[2, 20]),
               c(1.49071198499986, -0.74535599249993, 2.58198889747161, 2.26501206586678,
                 -3.59573259980396), tolerance = 1e-10)
  Z <- paired_score_z(20, -0.15)
  expect_equal(c(Z[1, 1], Z[1, 4], Z[5, 1], Z[8, 3], Z[20, 2]),
               c(-1.87867287325545, -2.65684465662029, 0.62622429108515, 0.701497013833065,
                 3.39248454773925), tolerance = 1e-10)
})

test_that("paired_score_z is the McNemar statistic for a zero margin", {
  Z <- paired_score_z(10, 0)
  expect_equal(Z[8, 3], 5 / 3)
  expect_equal(Z[1, 1], 0)
  expect_true(all(is.na(Z[outer(0:10, 0:10, '+') > 10])))
  expect_false(anyNA(Z[outer(0:10, 0:10, '+') <= 10]))
  # A small margin approaches the McNemar statistic
  expect_equal(paired_score_z(10, 1e-9)[8, 3], 5 / 3, tolerance = 1e-6)
})

test_that("paired_score_z gives sqrt(n Delta / (1 - Delta)) without discordant pairs", {
  # Tango (1998), formula (29)
  expect_equal(paired_score_z(30, 0.1)[1, 1], sqrt(30 * 0.1 / 0.9))
})
