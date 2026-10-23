test_that("strat_raw_n agrees with the reference values of Woolson et al. and Nam", {
  expect_equal(c(strat_raw_n(c(0.5, 0.7), c(0.3, 0.5), c(0.4, 0.6), c(1, 1), 0.025, 0.8,
                             'greater', FALSE),
                 strat_raw_n(c(0.5, 0.7), c(0.3, 0.5), c(0.4, 0.6), c(1, 1), 0.025, 0.8,
                             'greater', TRUE),
                 strat_raw_n(c(0.3, 0.5, 0.2), c(0.45, 0.6, 0.35), c(1, 2, 1) / 4,
                             c(2, 2, 2), 0.05, 0.9, 'less', FALSE),
                 strat_raw_n(c(0.6, 0.4), c(0.4, 0.3), c(0.5, 0.5), c(0.5, 2), 0.05, 0.8,
                             'two.sided', TRUE)),
               c(185.997689655091, 205.511097928377, 566.558428740696, 402.041885184174),
               tolerance = 1e-10)
})

test_that("strat_raw_n is vectorized over scenarios and NaN without an effect", {
  p1 <- rbind(c(0.5, 0.7), c(0.3, 0.5), c(0.6, 0.4))
  p2 <- rbind(c(0.3, 0.5), c(0.3, 0.5), c(0.4, 0.6))
  n <- strat_raw_n(p1, p2, c(0.5, 0.5), c(1, 1), 0.025, 0.8, 'greater', FALSE)
  expect_length(n, 3)
  expect_true(is.finite(n[1]))
  # No effect, and effects that cancel over the strata
  expect_true(is.nan(n[2]) && is.nan(n[3]))
  expect_true(is.nan(strat_raw_n(c(0.5, 0.7), c(0.3, 0.5), c(0.5, 0.5), c(1, 1), 0.025,
                                 0.8, 'less', FALSE)))
  # With one stratum the formula is that of two independent binomial proportions
  one <- strat_raw_n(c(0.5, 0.5), c(0.3, 0.3), c(0.5, 0.5), c(1, 1), 0.025, 0.8,
                     'greater', FALSE)
  z <- stats::qnorm(c(0.975, 0.8))
  ref <- 2 * (z[1] * sqrt(2 * 0.4 * 0.6) + z[2] * sqrt(0.5 * 0.5 + 0.3 * 0.7))^2 / 0.2^2
  expect_equal(one, ref)
})

test_that("strat_raw_n reproduces published sizes of Nam (1992) and Woolson et al. (1986)", {
  # Nam (1992, Section 2.1): the Iowa example with four strata, the odds ratio 3 and the
  # one-sided level 0.05; N = 192 and N* = 171 for the power 0.9, 115 and 95 for 0.7
  t <- c(0.1, 0.4, 0.35, 0.15)
  p0 <- c(0.75, 0.7, 0.65, 0.6)
  p1 <- p0 * 3 / (1 - p0 + p0 * 3)
  n <- c(strat_raw_n(p1, p0, t, rep(1, 4), 0.05, 0.9, 'greater', TRUE),
         strat_raw_n(p1, p0, t, rep(1, 4), 0.05, 0.9, 'greater', FALSE),
         strat_raw_n(p1, p0, t, rep(1, 4), 0.05, 0.7, 'greater', TRUE),
         strat_raw_n(p1, p0, t, rep(1, 4), 0.05, 0.7, 'greater', FALSE))
  expect_equal(round(n), c(192, 171, 115, 95))
  # Woolson et al. (1986, Table 1): two equal strata with the control probabilities
  # (pi21, 1 - pi21), the fractions (s1, 1 - s1) of group 1 and the quantiles 1.645 and
  # 0.841 of the article, rounded to the nearest integer
  wo <- function(pi21, theta, s1) {
    q <- c(pi21, 1 - pi21)
    s <- c(s1, 1 - s1)
    strat_raw_n(q * theta / (1 - q + q * theta), q, c(0.5, 0.5), s / (1 - s),
                stats::pnorm(-1.645), stats::pnorm(0.841), 'greater', FALSE)
  }
  n <- c(wo(0.1, 4, 0.1), wo(0.1, 2, 0.1), wo(0.2, 4, 0.1), wo(0.4, 2, 0.4))
  expect_equal(round(n), c(318, 1442, 205, 231))
})
