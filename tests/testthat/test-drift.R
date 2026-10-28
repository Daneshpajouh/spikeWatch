test_that("Shewhart limits use the moving range", {
    x <- c(1, 2, 1, 2, 1, 6)
    expect_warning(ch <- driftChart(x, baseline = 5, method = "shewhart"),
                   "fewer than 20")
    s <- 1 / 1.128
    expect_equal(attr(ch, "sigma"), s)
    expect_equal(ch$upper[1], 1.4 + 3 * s)
    expect_equal(firstSignal(ch), 6)
    expect_true(is.na(firstSignal(ch, after = 6)))
    expect_false(any(ch$signal[1:5]))
})

test_that("EWMA statistic and first limit match the formulas", {
    x <- c(1, 2, 1, 2, 1, 3, 2)
    ch <- suppressWarnings(
        driftChart(x, baseline = 5, method = "ewma", lambda = 0.2))
    expect_true(all(is.na(ch$stat[1:5])))
    expect_equal(ch$center[1:5], rep(1.4, 5))
    expect_equal(ch$stat[6], 0.2 * 3 + 0.8 * 1.4)
    expect_equal(ch$stat[7], 0.2 * 2 + 0.8 * ch$stat[6])
    s <- 1 / 1.128
    expect_equal(ch$upper[6] - 1.4, 3 * s * sqrt(0.2 / 1.8 * (1 - 0.8^2)))
})

test_that("EWMA starts after the last baseline run when runs are skipped", {
    # Baseline 1, 2, 4, 5 has mean 0.5; run 3 is outside it but before
    # its end, so it neither enters the statistic nor signals.
    x <- c(0, 1, 100, 0, 1, 50)
    ch <- suppressWarnings(driftChart(x, baseline = c(1, 2, 4, 5),
                                      method = "ewma", lambda = 0.2))
    expect_true(all(is.na(ch$stat[1:5])))
    expect_false(any(ch$signal[1:5]))
    expect_equal(ch$stat[6], 0.2 * 50 + 0.8 * 0.5)
    # Moving ranges stop at the gap: |1 - 0| twice, so sigma is 1 / 1.128,
    # and run 6 takes the first-step half-width, 3 sigma sqrt(0.2 / 1.8 * 0.36).
    expect_equal(ch$upper[6], 0.5 + 0.6 / 1.128)
})

test_that("CUSUM accumulates in sigma units from the first monitoring run", {
    x <- c(0, 1, 0, 1, 0, 3, 3)
    ch <- suppressWarnings(
        driftChart(x, baseline = 5, method = "cusum", sigma = "sd"))
    s <- stats::sd(c(0, 1, 0, 1, 0))
    z <- (3 - 0.4) / s - 0.5
    expect_true(all(is.na(ch$cusumUp[1:5])))
    expect_equal(ch$cusumUp[6], z)
    expect_equal(ch$cusumUp[7], 2 * z)
})

test_that("monitoring does not inherit what the baseline built up", {
    x <- c(rep(-1, 10), rep(1, 10), 0, 0)
    for (m in c("cusum", "ewma")) {
        ch <- driftChart(x, baseline = 20, method = m)
        expect_true(all(is.na(ch$stat[1:20])))
        expect_false(any(ch$signal))
        expect_equal(ch$stat[21], 0)
    }
})

test_that("a baseline covering every run leaves nothing to monitor", {
    x <- c(1, 2, 1, 2, 1)
    ch <- suppressWarnings(driftChart(x, baseline = 5, method = "cusum"))
    expect_equal(nrow(ch), 5L)
    expect_true(all(is.na(ch$stat)))
    expect_false(any(ch$signal))
    expect_true(is.na(firstSignal(ch)))
})

test_that("missing values carry the statistic and never signal", {
    x <- c(1, 2, 1, 2, 1, NA, 9)
    ch <- suppressWarnings(driftChart(x, baseline = 5, method = "ewma"))
    expect_equal(ch$stat[6], 1.4)
    expect_false(ch$signal[6])
    expect_true(ch$signal[7])
})

test_that("the default baseline is 20 runs and a short one warns", {
    set.seed(3)
    x <- rnorm(25)
    expect_silent(ch <- driftChart(x, method = "cusum"))
    expect_equal(sum(ch$baseline), 20L)
    expect_warning(driftChart(x, baseline = 19), "fewer than 20")
})

test_that("driftChart refuses a baseline it cannot use", {
    expect_error(suppressWarnings(driftChart(c(1, 1, 1, 1), baseline = 3)),
                 "no spread")
    expect_error(driftChart(1:5, baseline = 2), "at least three")
    expect_error(driftChart(1:5, baseline = 9), "outside")
    expect_error(driftChart(data.frame(a = 1:5)), "metric")
})

test_that("moving ranges never pair runs across a gap", {
    # only the ranges 10 to 12 and 12 to 10 are between neighbours
    x <- c(10, 12, NA, 12, 10, 11)
    ch <- suppressWarnings(driftChart(x, baseline = 5, method = "shewhart"))
    expect_equal(attr(ch, "sigma"), 2 / 1.128)
    y <- c(10, 12, 99, 12, 10, 11)
    ch2 <- suppressWarnings(driftChart(y, baseline = c(1, 2, 4, 5),
                                       method = "shewhart"))
    expect_equal(attr(ch2, "sigma"), 2 / 1.128)
    expect_error(suppressWarnings(
        driftChart(c(1, NA, 2, NA, 3, 4), baseline = 5)), "adjacent")
    # the standard deviation uses every observed baseline value
    ch3 <- suppressWarnings(driftChart(x, baseline = 5, sigma = "sd"))
    expect_equal(attr(ch3, "sigma"), stats::sd(c(10, 12, 12, 10)))
})

test_that("the lower CUSUM catches a fall", {
    ch <- suppressWarnings(driftChart(c(-1, 0, 1, -6), baseline = 3,
                                      method = "cusum", sigma = "sd"))
    expect_equal(ch$cusumLow[4], 5.5)
    expect_equal(ch$cusumUp[4], 0)
    expect_true(ch$signal[4])
})

test_that("firstSignal refuses a chart without a baseline", {
    ch <- suppressWarnings(driftChart(c(1, 2, 1, 2, 1, 6), baseline = 5,
                                      method = "shewhart"))
    ch$baseline <- FALSE
    expect_error(firstSignal(ch), "no baseline")
})

test_that("driftChart takes a named limit with attributes", {
    set.seed(3)
    x <- c(rnorm(20), rnorm(10, 2))
    h <- structure(c(h = 5), alpha = 0.05, method = "cusum")
    L <- structure(c(L = 3), alpha = 0.05, method = "ewma")
    a <- driftChart(x, baseline = 20, method = "cusum", h = h)
    b <- driftChart(x, baseline = 20, method = "cusum", h = 5)
    expect_equal(a$signal, b$signal)
    expect_equal(nrow(a), 30)
    e <- driftChart(x, baseline = 20, method = "ewma", L = L)
    expect_equal(e$upper, driftChart(x, baseline = 20, L = 3)$upper)
})
