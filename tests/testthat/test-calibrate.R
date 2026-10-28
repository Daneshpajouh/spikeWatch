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

test_that("calibrateChart matches a plain loop on the same draws", {
    x <- c(0.3, -1.2, 0.8, 0.1, -0.4, 1.5, -0.9, 0.2, 0.6, -0.1)
    n <- length(x)
    B <- 400
    H <- 6
    ref <- function(method, alpha) {
        set.seed(11)
        xs <- matrix(x[sample.int(n, n * B, replace = TRUE)], B)
        fut <- matrix(x[sample.int(n, B * H, replace = TRUE)], B)
        m <- numeric(B)
        for (b in seq_len(B)) {
            mu <- mean(xs[b, ])
            s <- mean(abs(diff(xs[b, ]))) / 1.128
            cu <- cl <- z <- 0
            for (t in seq_len(H)) {
                u <- (fut[b, t] - mu) / s
                if (method == "cusum") {
                    cu <- max(0, cu + u - 0.5)
                    cl <- max(0, cl - u - 0.5)
                    m[b] <- max(m[b], cu, cl)
                } else {
                    z <- 0.2 * u + 0.8 * z
                    w <- sqrt(0.2 / 1.8 * (1 - 0.8^(2 * t)))
                    m[b] <- max(m[b], abs(z) / w)
                }
            }
        }
        sort(m)[ceiling((1 - alpha) * B)]
    }
    for (method in c("cusum", "ewma")) for (alpha in c(0.05, 0.2)) {
        set.seed(11)
        got <- suppressWarnings(calibrateChart(x, baseline = n,
            method = method, alpha = alpha, horizon = H, B = B))
        expect_equal(unname(c(got)), ref(method, alpha))
    }
})

test_that("bootstrap baselines with no spread count as alarms", {
    set.seed(3)
    expect_warning(h <- calibrateChart(c(1, 1, 2, rnorm(17)), baseline = 3,
                                       B = 400, horizon = 5), "no spread")
    expect_true(is.infinite(h))
    expect_warning(e <- calibrateChart(c(1, 1, 2, rnorm(17)), baseline = 3,
                                       B = 400, horizon = 5, method = "ewma"),
                   "no spread")
    expect_true(is.infinite(e))
})
