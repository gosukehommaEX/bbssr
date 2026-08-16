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

test_that("the plotting methods keep their previous defaults", {
  pw <- plot(BinaryPower(seq(0.3, 0.7, by = 0.1), rep(0.2, 5), 10, 10, 0.025, 'Chisq'))
  expect_identical(pw$labels$x, 'Response probability of group 1')
  expect_identical(pw$labels$y, 'Power')
  expect_identical(pw$theme$text$size, 11)
  yr <- ggplot2::ggplot_build(pw)$layout$panel_params[[1]]$y.range
  expect_lt(yr[1], 0)
  expect_gt(yr[2], 1)
  rr <- plot(BinaryRR(6, 6, 0.025, 'Chisq'))
  expect_null(rr$labels$fill)
  ss <- plot(BinarySampleSize(0.5, 0.2, 1, 0.025, 0.8, 'Chisq'), N2.range = 20:24)
  expect_length(ggplot2::ggplot_build(ss)$data, 4)
  expect_true(grepl('target power', ss$labels$subtitle, fixed = TRUE))
})

test_that("the labels and the base font size can be replaced", {
  pw <- BinaryPower(seq(0.3, 0.7, by = 0.1), rep(0.2, 5), 10, 10, 0.025, 'Chisq')
  p <- plot(pw, main = 'Custom title', sub = NA, xlab = 'x axis', ylab = 'y axis',
            ylim = c(0.2, 0.9), colours = 'darkgreen', base_size = 14)
  expect_identical(p$labels$title, 'Custom title')
  expect_null(p$labels$subtitle)
  expect_identical(p$labels$x, 'x axis')
  expect_identical(p$labels$y, 'y axis')
  expect_identical(p$theme$text$size, 14)
  built <- ggplot2::ggplot_build(p)
  expect_true(all(built$data[[1]]$colour == 'darkgreen'))
  yr <- built$layout$panel_params[[1]]$y.range
  expect_true(yr[1] > 0.1 && yr[1] < 0.2)
  expect_true(yr[2] > 0.9 && yr[2] < 1)
})

test_that("the markers can be dropped so that only the curve is drawn", {
  pw <- plot(BinaryPower(seq(0.3, 0.7, by = 0.1), rep(0.2, 5), 10, 10, 0.025, 'Chisq'),
             show.points = FALSE)
  expect_length(ggplot2::ggplot_build(pw)$data, 1)
  expect_true(inherits(pw$layers[[1]]$geom, 'GeomLine'))
  res <- BinaryPowerBSSR(p = c(0.4, 0.45, 0.5),
                         Delta.A = 0.3, Delta.T = 0.3, N1 = 6, N2 = 6,
                         omega = 0.5, r = 1, alpha = 0.025,
                         tar.power = 0.8, Test = 'Chisq')
  bssr <- plot(res, show.points = FALSE)
  expect_length(ggplot2::ggplot_build(bssr)$data, 2)
  ss <- plot(BinarySampleSize(0.5, 0.2, 1, 0.025, 0.8, 'Chisq'),
             N2.range = 20:24, show.points = FALSE)
  expect_length(ggplot2::ggplot_build(ss)$data, 3)
})

test_that("the legend of the rejection region plot can be titled and relabelled", {
  p <- plot(BinaryRR(6, 6, 0.025, 'Chisq'), legend.title = 'Decision',
            legend.labels = c('keep', 'drop'), colours = c('white', 'black'))
  expect_identical(p$labels$fill, 'Decision')
  # A scale is trained by the build, and calling get_labels() on the plot object queries
  # it before that has happened, so the labels are read from the scale itself
  expect_identical(p$scales$get_scales('fill')$labels, c('keep', 'drop'))
  expect_true(all(unique(ggplot2::ggplot_build(p)$data[[1]]$fill) %in% c('white', 'black')))
})

test_that("the reference line of the BSSR plot can be moved and dropped", {
  res <- BinaryPowerBSSR(p = c(0.4, 0.45, 0.5),
                         Delta.A = 0.3, Delta.T = 0, N1 = 6, N2 = 6,
                         omega = 0.5, r = 1, alpha = 0.025,
                         tar.power = 0.8, Test = 'Chisq')
  moved <- plot(res, ref.line = 0.025, ylab = 'Type I error rate',
                legend.labels = c('BSSR', 'Fixed'))
  built <- ggplot2::ggplot_build(moved)
  expect_length(built$data, 3)
  expect_equal(built$data[[3]]$yintercept, 0.025)
  expect_identical(moved$labels$y, 'Type I error rate')
  expect_true(grepl('dashed line at 0.025', moved$labels$subtitle, fixed = TRUE))
  expect_identical(moved$scales$get_scales('colour')$labels, c('BSSR', 'Fixed'))
  dropped <- plot(res, ref.line = NA)
  expect_length(ggplot2::ggplot_build(dropped)$data, 2)
  expect_false(grepl('dashed', dropped$labels$subtitle, fixed = TRUE))
})

test_that("both reference lines of the sample size plot can be moved and dropped", {
  ss <- BinarySampleSize(0.5, 0.2, 1, 0.025, 0.8, 'Chisq')
  moved <- plot(ss, N2.range = 20:24, ref.line = 0.9, ref.line.N2 = ss$N2[1] + 1)
  built <- ggplot2::ggplot_build(moved)
  expect_equal(built$data[[3]]$yintercept, 0.9)
  expect_equal(built$data[[4]]$xintercept, ss$N2[1] + 1)
  expect_true(grepl('dashed line at 0.9', moved$labels$subtitle, fixed = TRUE))
  bare <- plot(ss, N2.range = 20:24, ref.line = NA, ref.line.N2 = NA)
  expect_length(ggplot2::ggplot_build(bare)$data, 2)
  expect_identical(bare$labels$subtitle, sprintf('Selected N2 = %d', ss$N2[1]))
})
