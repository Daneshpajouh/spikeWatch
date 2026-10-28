test_that("Shewhart limits use the moving range", {
    x <- c(1, 2, 1, 2, 1, 6)
    ch <- driftChart(x, baseline = 5, method = "shewhart")
    s <- 1 / 1.128
    expect_equal(attr(ch, "sigma"), s)
    expect_equal(ch$upper[1], 1.4 + 3 * s)
    expect_equal(firstSignal(ch), 6)
    expect_false(any(ch$signal[1:5]))
})

test_that("EWMA statistic and first limit match the formulas", {
    x <- c(1, 2, 1, 2, 1, 3)
    ch <- driftChart(x, baseline = 5, method = "ewma", lambda = 0.2)
    expect_equal(ch$stat[1], 0.2 * 1 + 0.8 * 1.4)
    s <- 1 / 1.128
    expect_equal(ch$upper[1] - 1.4, 3 * s * sqrt(0.2 / 1.8 * (1 - 0.8^2)))
})

test_that("CUSUM accumulates in sigma units", {
    x <- c(0, 1, 0, 1, 0, 3, 3)
    ch <- driftChart(x, baseline = 5, method = "cusum", sigma = "sd")
    s <- stats::sd(c(0, 1, 0, 1, 0))
    z <- (3 - 0.4) / s - 0.5
    expect_equal(ch$cusumUp[7] - ch$cusumUp[5], 2 * z)
})

test_that("missing values carry the statistic and never signal", {
    x <- c(1, 2, 1, 2, 1, NA, 9)
    ch <- driftChart(x, baseline = 5, method = "ewma")
    expect_equal(ch$stat[6], ch$stat[5])
    expect_false(ch$signal[6])
    expect_true(ch$signal[7])
})

test_that("driftChart refuses a baseline it cannot use", {
    expect_error(driftChart(c(1, 1, 1, 1), baseline = 3), "no spread")
    expect_error(driftChart(1:5, baseline = 2), "at least three")
    expect_error(driftChart(1:5, baseline = 9), "outside")
    expect_error(driftChart(data.frame(a = 1:5)), "metric")
})
