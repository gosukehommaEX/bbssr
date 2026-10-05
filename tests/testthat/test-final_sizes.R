re_of <- function(n.raw, N1.re, N2.re) data.frame(n.raw = n.raw, N1.re = N1.re, N2.re = N2.re)

test_that("the group rule keeps the interim size and the allocation ratio", {
  re <- re_of(NA, c(10, 40, 90), c(5, 20, 45))
  fin <- final_sizes(re, 16, 8, 2, 'group', FALSE, NULL, NULL, NULL, NULL)
  expect_identical(fin$N2, c(8L, 20L, 45L))
  expect_identical(fin$N1, c(16L, 40L, 90L))
})

test_that("the group rule applies the restricted rule and the bounds", {
  re <- re_of(NA, c(10, 40, 90), c(5, 20, 45))
  fin <- final_sizes(re, 16, 8, 2, 'group', TRUE, 60, 30, NULL, NULL)
  expect_identical(fin$N2, c(30L, 30L, 45L))
  fin <- final_sizes(re, 16, 8, 2, 'group', FALSE, NULL, NULL, 36, 100)
  # ceiling(36 / 3) = 12 from below, and 33 + 66 = 99 is the largest total within 100
  expect_identical(fin$N2, c(12L, 20L, 33L))
  expect_identical(fin$N1, c(24L, 40L, 66L))
  expect_true(all(fin$N1 + fin$N2 <= 100))
  expect_error(final_sizes(re, 16, 8, 2, 'group', FALSE, NULL, NULL, NULL, 20), 'N.max')
})

test_that("the group rule respects N.max for a fractional allocation ratio", {
  re <- re_of(NA, 1000, 1000)
  for (N.max in 50:60) {
    fin <- final_sizes(re, 3, 2, 1.5, 'group', FALSE, NULL, NULL, NULL, N.max)
    expect_lte(fin$N1 + fin$N2, N.max)
    # One more patient in group 2 would exceed the bound
    expect_gt(fin$N2 + 1 + ceiling(1.5 * (fin$N2 + 1)), N.max)
  }
})

test_that("the Friede-Kieser rule rounds the second stage up group by group", {
  # Interim 30 + 10, unrounded totals 94.19 and 250, cap 190
  re <- re_of(c(94.19, 250, 20), NA, NA)
  fin <- final_sizes(re, 30, 10, 3, 'friede-kieser', FALSE, NULL, NULL, NULL, 190)
  # 94.19: n2 = 55, 14 + 42 more; 250 -> 190: n2 = 150, 38 + 113 more; 20: no stage 2
  expect_identical(fin$N2, c(24L, 48L, 10L))
  expect_identical(fin$N1, c(72L, 143L, 30L))
})

test_that("the total rule rounds the total and splits it", {
  re <- re_of(c(312.4, 100), NA, NA)
  fin <- final_sizes(re, 79, 79, 1, 'total', FALSE, NULL, NULL, NULL, NULL)
  expect_identical(fin$N1 + fin$N2, c(313L, 158L))
  expect_identical(fin$N2, c(156L, 79L))
})

test_that("final_sizes rounds each group to the nearest whole number under 'nearest'", {
  re <- re_of(c(50, 101.2, 30), c(34L, 68L, 20L), c(17L, 34L, 10L))
  fin <- final_sizes(re, 30, 15, 2, 'nearest', FALSE, NULL, NULL, NULL, NULL)
  expect_equal(fin$N1, c(33L, 67L, 30L))
  expect_equal(fin$N2, c(17L, 34L, 15L))
  fin <- final_sizes(re, 30, 15, 2, 'nearest', FALSE, NULL, NULL, NULL, 60)
  expect_equal(fin$N1, c(33L, 40L, 30L))
  expect_equal(fin$N2, c(17L, 20L, 15L))
})
