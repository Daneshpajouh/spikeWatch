test_that("calibrated CUSUM limits hold the false-alarm target on normal data", {
    set.seed(5)
    fa <- vapply(seq_len(300), function(i) {
        x <- rnorm(30)
        h <- calibrateChart(x, baseline = 20, method = "cusum",
                            horizon = 10, B = 400)
        ch <- driftChart(x, baseline = 20, method = "cusum", h = h)
        !is.na(firstSignal(ch))
    }, logical(1))
    expect_lt(mean(fa), 0.10)
    expect_gt(mean(fa), 0.01)
})

test_that("calibrated limits are wider than the defaults with a short baseline", {
    set.seed(2)
    x <- rnorm(20)
    h <- calibrateChart(x, baseline = 20, method = "cusum", horizon = 20,
                        B = 2000)
    expect_named(h, "h")
    expect_gt(unname(h), 5)
    L <- calibrateChart(x, baseline = 20, method = "ewma", horizon = 20,
                        B = 2000)
    expect_named(L, "L")
    expect_gt(unname(L), 3)
    expect_equal(attr(L, "horizon"), 20)
})

test_that("the bootstrap maxima match driftChart's statistics", {
    u <- matrix(c(0.5, 1.5, -2, 0.2, 3, -1), nrow = 1)
    expect_warning(ch <- driftChart(c(-1, 1, -1, 1, as.vector(u)),
                                    baseline = 4, method = "cusum",
                                    h = 100, sigma = "sd"), "fewer than 20")
    s <- attr(ch, "sigma")
    expect_equal(spikeWatch:::.cusumMax(u / s, 0.5),
                 max(ch$stat, na.rm = TRUE))
    expect_warning(ce <- driftChart(c(-1, 1, -1, 1, as.vector(u)),
                                    baseline = 4, method = "ewma", L = 1,
                                    sigma = "sd"), "fewer than 20")
    mon <- 5:10
    expect_equal(spikeWatch:::.ewmaMax(u / s, 0.2),
                 max(abs(ce$stat[mon]) / (ce$upper[mon] - ce$center[mon])))
})

test_that("calibrateChart refuses bad inputs", {
    expect_error(calibrateChart(rnorm(30), alpha = 1), "alpha")
    expect_error(calibrateChart(rnorm(30), B = 10), "B")
    expect_error(calibrateChart(c(1, NA, NA, 2), baseline = 4), "three")
    expect_error(calibrateChart(data.frame(a = 1:30)), "metric")
})
