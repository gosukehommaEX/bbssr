re.args <- list(Delta.A = 0.3, effect = 'RD', frac = c(0.5, 0.5), r = c(1, 1),
                ss.method = 'standard', ss.Test = 'CMH', ss.alpha = 0.05, tar.power = 0.8,
                N1.plan = c(10L, 10L), N2.plan = c(10L, 10L), restricted = FALSE,
                N.min = NULL, N.max = NULL, search = 'crossing', search.limit = c(2, 50))

test_that("strat_map lists every interim outcome with the first stratum varying fastest", {
  m <- strat_map(c(5L, 4L), c(5L, 3L), 0.05, 'CMH', 'greater', 'minlike', re.args)
  expect_identical(nrow(m$map), 11L * 8L)
  expect_identical(m$map$s1[1:12], c(0:10, 0L))
  expect_identical(m$map$s2[1:12], c(rep(0L, 11), 1L))
  expect_identical(m$map$N, as.integer(m$map$N1.S1 + m$map$N1.S2 + m$map$N2.S1 +
                                         m$map$N2.S2))
  expect_length(m$id, nrow(m$map))
  # Each identifier points to the design of its row, with its critical values
  for (i in c(1, 17, 50, 88)) {
    d <- m$designs[[m$id[i] + 1]]
    expect_identical(d$N1, c(m$map$N1.S1[i], m$map$N1.S2[i]))
    expect_identical(d$N2, c(m$map$N2.S1[i], m$map$N2.S2[i]))
    expect_length(d$lo, prod(d$N1 + d$N2 + 1))
  }
  expect_identical(length(m$designs), length(unique(m$id)))
})

test_that("strat_map gives the rows of strat_reestimate", {
  m <- strat_map(c(5L, 4L), c(5L, 3L), 0.05, 'CMH', 'greater', 'minlike', re.args)
  re <- do.call(strat_reestimate, c(list(s = cbind(c(3, 7), c(5, 0)), n1 = c(5L, 4L),
                                         n2 = c(5L, 3L), alternative = 'greater',
                                         tsmethod = 'minlike'), re.args))
  rows <- c(3 + 1 + 11 * 5, 7 + 1)
  expect_identical(cbind(m$map$N1.S1[rows], m$map$N1.S2[rows]), re$N1)
  expect_equal(m$map$n.raw[rows], re$n.raw)
})
