test_that("bssr_setup lists every interim cell with its final sample size", {
  map <- bssr_map(0.3, 10, 10, NULL, c(4, 3), 1, 0.025, 0.8, 'Chisq', FALSE, 'greater',
                  'minlike', 100L, 0, 'RD', 'standard', 'Chisq', 0.025, 'group', NULL,
                  NULL)
  st <- bssr_setup(map)
  expect_equal(length(st$x11), 5 * 4)
  expect_equal(st$x11[1:6], c(0:4, 0))
  expect_equal(st$x12[1:6], c(rep(0, 5), 1))
  # The final sizes reached from each cell are those of its pooled count
  expect_equal(st$N1[st$rr.id + 1L], map$N1[st$x11 + st$x12 + 1L])
  expect_equal(st$N2[st$rr.id + 1L], map$N2[st$x11 + st$x12 + 1L])
  expect_false(anyDuplicated(paste(st$N1, st$N2)) > 0)
  expect_equal(st$n21, st$N1 - 4L)
  expect_equal(st$n22, st$N2 - 3L)
})
