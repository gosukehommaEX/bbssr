test_that("bssr_map agrees with BinaryBSSR for every interim total", {
  map <- bssr_map(0.3, 12, 12, 0.5, NULL, 1, 0.025, 0.8, 'Chisq', FALSE, 'greater',
                  'minlike', 100L, 0, 'RD', 'exact', 'Chisq', 0.025, 'group', NULL, NULL,
                  FALSE, 0, margin.scale = 'RD')
  expect_equal(attr(map, 'n11'), 6L)
  expect_equal(attr(map, 'n12'), 6L)
  expect_equal(map$s, 0:12)
  for (s in c(0, 5, 12)) {
    b <- BinaryBSSR(6, 6, s, 0.3, 1, 0.025, 0.8, 'Chisq')
    row <- map[map$s == s, ]
    expect_equal(c(row$N1, row$N2), c(b$N1.final, b$N2.final), info = s)
    expect_equal(c(row$hat.p1, row$hat.p2), c(b$hat.p1, b$hat.p2), info = s)
  }
})

test_that("n.interim and omega describe the same interim analysis", {
  args <- list(Delta.A = 0.3, N1 = 20, N2 = 10, r = 2, alpha = 0.025, tar.power = 0.8,
               Test = 'Chisq', restricted = FALSE, alternative = 'greater',
               tsmethod = 'minlike', n.grid = 100L, bb.gamma = 0, effect = 'RD',
               ss.method = 'standard', ss.Test = 'Chisq', ss.alpha = 0.025,
               rounding = 'group', N.min = NULL, N.max = NULL, ref.pvalue = FALSE,
               margin = 0, margin.scale = 'RD')
  a <- do.call(bssr_map, c(args, list(omega = 0.5, n.interim = NULL)))
  b <- do.call(bssr_map, c(args, list(omega = NULL, n.interim = c(10, 5))))
  expect_identical(a, b)
  expect_error(do.call(bssr_map, c(args, list(omega = 0.5, n.interim = c(10, 5)))),
               'not both')
  expect_error(do.call(bssr_map, c(args, list(omega = NULL, n.interim = NULL))),
               'either omega')
  expect_error(do.call(bssr_map, c(args, list(omega = NULL, n.interim = c(10, 5.5)))),
               'n.interim')
})

test_that("bssr_map rejects an exact re-estimation with another rounding rule", {
  expect_error(bssr_map(0.3, 12, 12, 0.5, NULL, 1, 0.025, 0.8, 'Chisq', FALSE, 'greater',
                        'minlike', 100L, 0, 'RD', 'exact', 'Chisq', 0.025, 'total', NULL,
                        NULL, FALSE, 0, margin.scale = 'RD'), 'rounding')
})

test_that("the interim sizes are not raised by a rounding error of the products", {
  args <- list(Delta.A = 0.3, n.interim = NULL, alpha = 0.025, tar.power = 0.8,
               Test = 'Chisq', restricted = FALSE, alternative = 'greater',
               tsmethod = 'minlike', n.grid = 100L, bb.gamma = 0, effect = 'RD',
               ss.method = 'standard', ss.Test = 'Chisq', ss.alpha = 0.025,
               rounding = 'group', N.min = NULL, N.max = NULL, ref.pvalue = FALSE,
               margin = 0, margin.scale = 'RD')
  # 0.55 * 100 exceeds 55, and 1.1 * 50 exceeds 55, by a rounding error
  a <- do.call(bssr_map, c(args, list(N1 = 100, N2 = 100, omega = 0.55, r = 1)))
  expect_identical(c(attr(a, 'n11'), attr(a, 'n12')), c(55L, 55L))
  b <- do.call(bssr_map, c(args, list(N1 = 55, N2 = 50, omega = 1, r = 1.1)))
  expect_identical(c(attr(b, 'n11'), attr(b, 'n12')), c(55L, 50L))
  args$N.max <- Inf
  expect_error(do.call(bssr_map, c(args, list(N1 = 20, N2 = 20, omega = 0.5, r = 1))),
               'N.max must be a single finite number')
})
