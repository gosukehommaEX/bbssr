test_that("print methods return their argument invisibly and produce output", {
  objs <- list(
    rr = BinaryRR(5, 5, 0.025, 'Chisq'),
    power = BinaryPower(0.5, 0.2, 8, 8, 0.025, 'Chisq'),
    samplesize = BinarySampleSize(0.5, 0.2, 1, 0.025, 0.8, 'Chisq'),
    powerbssr = BinaryPowerBSSR(p = 0.45,
                                Delta.A = 0.3, Delta.T = 0.3, N1 = 6, N2 = 6,
                                omega = 0.5, r = 1, alpha = 0.025,
                                tar.power = 0.8, Test = 'Chisq'),
    bssr = BinaryBSSR(20, 20, 11, 0.3, 1, 0.025, 0.8, 'Chisq')
  )
  for (nm in names(objs)) {
    out <- utils::capture.output(res <- print(objs[[nm]]))
    expect_gt(length(out), 3, label = nm)
    expect_identical(res, objs[[nm]], label = nm)
  }
})

test_that("the rejection region map is shown for small grids and hidden for large ones", {
  small <- utils::capture.output(print(BinaryRR(5, 5, 0.025, 'Chisq')))
  expect_true(any(grepl('x1=0', small, fixed = TRUE)))
  large <- utils::capture.output(print(BinaryRR(40, 40, 0.025, 'Chisq')))
  expect_true(any(grepl('suppressed', large, fixed = TRUE)))
  forced <- utils::capture.output(print(BinaryRR(40, 40, 0.025, 'Chisq'), show.map = TRUE))
  expect_true(any(grepl('x1=0', forced, fixed = TRUE)))
})

test_that("plot methods return ggplot objects", {
  plots <- list(
    plot(BinaryRR(6, 6, 0.025, 'Chisq')),
    plot(BinaryPower(seq(0.3, 0.7, by = 0.1), rep(0.2, 5), 10, 10, 0.025, 'Chisq')),
    plot(BinarySampleSize(0.5, 0.2, 1, 0.025, 0.8, 'Chisq'), N2.range = 20:24),
    plot(BinaryPowerBSSR(p = c(0.4, 0.45, 0.5),
                         Delta.A = 0.3, Delta.T = 0.3, N1 = 6, N2 = 6,
                         omega = 0.5, r = 1, alpha = 0.025,
                         tar.power = 0.8, Test = 'Chisq'))
  )
  for (p in plots) {
    expect_s3_class(p, 'ggplot')
    expect_no_error(ggplot2::ggplot_build(p))
  }
})

test_that("the printed output reports the design settings", {
  out <- utils::capture.output(
    print(BinaryRR(6, 6, 0.05, 'Boschloo', alternative = 'two.sided',
                   n.grid = 40, bb.gamma = 0.001))
  )
  expect_true(any(grepl('Boschloo', out)))
  expect_true(any(grepl('two.sided', out)))
  expect_true(any(grepl('minlike', out)))
  expect_true(any(grepl('Berger-Boos', out)))
})

test_that("the rejection region plot uses integer breaks and matches the printed layout", {
  p <- plot(BinaryRR(10, 10, 0.025, 'Chisq'))
  built <- ggplot2::ggplot_build(p)
  x.breaks <- built$layout$panel_params[[1]]$x$get_breaks()
  y.breaks <- built$layout$panel_params[[1]]$y$get_breaks()
  x.breaks <- x.breaks[is.finite(x.breaks)]
  y.breaks <- y.breaks[is.finite(y.breaks)]
  expect_equal(x.breaks, round(x.breaks))
  expect_equal(y.breaks, round(y.breaks))
  # The vertical axis runs downwards, so no responders in group 1 sits at the top.
  # The layer keeps the row order of expand.grid(x1 = 0:10, x2 = 0:10), so row 1 holds
  # the outcome x1 = 0 and row 11 the outcome x1 = 10
  tiles <- built$data[[1]]
  expect_gt(tiles$y[1], tiles$y[11])
})

test_that("the sample size plot uses integer breaks", {
  p <- plot(BinarySampleSize(0.5, 0.2, 1, 0.025, 0.8, 'Chisq'), N2.range = 20:30)
  built <- ggplot2::ggplot_build(p)
  br <- built$layout$panel_params[[1]]$x$get_breaks()
  br <- br[is.finite(br)]
  expect_equal(br, round(br))
})
