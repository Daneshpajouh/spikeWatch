test_that("simulated truth equals the design and drift labels are right", {
    d <- erccDesign()
    set.seed(3)
    s <- simulateSpikeRuns(d, nRuns = 8, mix = "Mix2", driftAt = 5,
                           driftType = "offset", driftSize = -1)
    e <- SummarizedExperiment::assay(s, "expected")
    m2 <- d[d$mix == "Mix2", ]
    expect_equal(unname(e[m2$feature, 1]), m2$amount)
    expect_equal(s$drifted, rep(c(FALSE, TRUE), each = 4))
    expect_true(all(SummarizedExperiment::assay(s, "counts") >= 0))
})

test_that("a planted slope drift is found soon after it starts", {
    set.seed(11)
    s <- simulateSpikeRuns(erccDesign(), nRuns = 40, driftAt = 31,
                           driftType = "slope", driftSize = -0.08)
    m <- runMetrics(s, minAmount = 1)
    ch <- driftChart(m, "slope", baseline = 20, method = "cusum")
    hit <- firstSignal(ch)
    expect_false(is.na(hit))
    expect_gte(hit, 31)
    expect_lte(hit, 36)
})

test_that("a ramp reaches full size after rampRuns", {
    set.seed(1)
    s <- simulateSpikeRuns(erccDesign(), nRuns = 12, driftAt = 5,
                           driftShape = "ramp", rampRuns = 4,
                           driftType = "offset", driftSize = 2)
    expect_equal(s$drifted, c(rep(FALSE, 4), rep(TRUE, 8)))
})

test_that("departures are off by default and change only what they name", {
    d <- erccDesign()
    set.seed(9)
    a <- simulateSpikeRuns(d, nRuns = 6)
    set.seed(9)
    b <- simulateSpikeRuns(d, nRuns = 6, pivot = "median",
                           counts = "poisson", saturation = 0,
                           spikeBiasSD = 0, slopeSD = 0, runCor = 0)
    expect_identical(SummarizedExperiment::assay(a, "counts"),
                     SummarizedExperiment::assay(b, "counts"))
    expect_error(simulateSpikeRuns(d, runCor = 1), "runCor")
    expect_error(simulateSpikeRuns(d, saturation = 1), "saturation")
})

test_that("a top pivot raises the spikes below the top under compression", {
    d <- erccDesign()
    run <- function(pivot) {
        set.seed(2)
        s <- simulateSpikeRuns(d, nRuns = 2, driftAt = 2, driftSize = -0.3,
                               noise = 0, depthSD = 0, pivot = pivot,
                               scale = 1000)
        cnt <- SummarizedExperiment::assay(s, "counts")
        e <- SummarizedExperiment::assay(s, "expected")[, 1]
        list(r = log2((cnt[, 2] + 1) / (cnt[, 1] + 1)), e = e)
    }
    top <- run("top")
    mid <- run("median")
    upper <- top$e > stats::median(top$e) & top$e < max(top$e)
    expect_gt(stats::median(top$r[upper]), 0.2)
    expect_lt(stats::median(mid$r[upper]), -0.2)
    expect_lt(abs(stats::median(top$r[top$e == max(top$e)])), 0.1)
})

test_that("saturation compresses above log2(scale) plus the log2 quantile", {
    # Amounts 1 and 100, slope 1: log2 means are 0 and log2(100) above
    # log2(scale). The ceiling is 0.9 * log2(100) above it, and half the
    # excess is kept. The quantile of the amounts themselves would give
    # 6.57 for the top spike instead of 6.31.
    d <- spikeDesign(c("lo", "hi"), c(1, 100))
    set.seed(11)
    s <- simulateSpikeRuns(d, nRuns = 1, noise = 0, depthSD = 0,
                           scale = 1e4, saturation = 0.5)
    y <- log2(SummarizedExperiment::assay(s, "counts")[, 1] / 1e4)
    top <- 0.9 * log2(100)
    expect_lt(abs(y[["hi"]] - (top + 0.5 * (log2(100) - top))), 0.02)
    expect_lt(abs(y[["lo"]]), 0.05)
})

test_that("saturation flattens the top and nbinom spreads the counts", {
    d <- erccDesign()
    set.seed(6)
    flat <- simulateSpikeRuns(d, nRuns = 1, noise = 0, depthSD = 0,
                              saturation = 0.5, scale = 100)
    set.seed(6)
    full <- simulateSpikeRuns(d, nRuns = 1, noise = 0, depthSD = 0,
                              scale = 100)
    cf <- SummarizedExperiment::assay(flat, "counts")[, 1]
    cu <- SummarizedExperiment::assay(full, "counts")[, 1]
    e <- SummarizedExperiment::assay(full, "expected")[, 1]
    hi <- e == max(e)
    expect_true(all(cf[hi] < cu[hi]))
    expect_equal(cf[e < stats::quantile(e, 0.5)],
                 cu[e < stats::quantile(e, 0.5)], tolerance = 0.5)
    set.seed(7)
    p <- simulateSpikeRuns(d, nRuns = 30, noise = 0, depthSD = 0,
                           scale = 100)
    set.seed(7)
    n <- simulateSpikeRuns(d, nRuns = 30, noise = 0, depthSD = 0,
                           scale = 100, counts = "nbinom", size = 5)
    v <- function(s) stats::median(apply(
        SummarizedExperiment::assay(s, "counts")[hi, , drop = FALSE], 1,
        stats::var) / rowMeans(
        SummarizedExperiment::assay(s, "counts")[hi, , drop = FALSE]))
    expect_gt(v(n), 10 * v(p))
})

test_that("spike biases persist across runs and runCor links runs", {
    d <- erccDesign()
    set.seed(12)
    s <- simulateSpikeRuns(d, nRuns = 2, noise = 0, depthSD = 0,
                           spikeBiasSD = 1, scale = 1000)
    cnt <- SummarizedExperiment::assay(s, "counts")
    e <- SummarizedExperiment::assay(s, "expected")[, 1]
    ok <- e > stats::median(e)
    b1 <- log2(cnt[ok, 1] / (1000 * e[ok]))
    b2 <- log2(cnt[ok, 2] / (1000 * e[ok]))
    expect_gt(stats::cor(b1, b2), 0.95)
    expect_gt(stats::sd(b1), 0.5)
    set.seed(13)
    r <- simulateSpikeRuns(d, nRuns = 400, runCor = 0.8)
    expect_gt(stats::cor(r$depth[-1], r$depth[-400]), 0.65)
    expect_lt(abs(stats::sd(r$depth) - 0.5), 0.12)
})

test_that("a ramp plants the change linearly up to full size", {
    set.seed(1)
    s <- simulateSpikeRuns(erccDesign(), nRuns = 12, driftAt = 5,
                           driftShape = "ramp", rampRuns = 4,
                           driftType = "offset", driftSize = 2,
                           depthSD = 0)
    expect_equal(s$depth, c(0, 0, 0, 0, 0.5, 1, 1.5, 2, 2, 2, 2, 2))
    set.seed(1)
    st <- simulateSpikeRuns(erccDesign(), nRuns = 6, driftAt = 3,
                            driftType = "offset", driftSize = 2,
                            depthSD = 0)
    expect_equal(st$depth, c(0, 0, 2, 2, 2, 2))
})
