test_that("simulated truth equals the design and drift labels are right", {
    d <- erccDesign()
    s <- simulateSpikeRuns(d, nRuns = 8, mix = "Mix2", driftAt = 5,
                           driftType = "offset", driftSize = -1, seed = 3)
    e <- SummarizedExperiment::assay(s, "expected")
    m2 <- d[d$mix == "Mix2", ]
    expect_equal(unname(e[m2$feature, 1]), m2$amount)
    expect_equal(s$drifted, rep(c(FALSE, TRUE), each = 4))
    expect_true(all(SummarizedExperiment::assay(s, "counts") >= 0))
})

test_that("a planted slope drift is found soon after it starts", {
    s <- simulateSpikeRuns(erccDesign(), nRuns = 40, driftAt = 31,
                           driftType = "slope", driftSize = -0.08,
                           seed = 11)
    m <- runMetrics(s, minAmount = 1)
    ch <- driftChart(m, "slope", baseline = 20, method = "cusum")
    hit <- firstSignal(ch)
    expect_false(is.na(hit))
    expect_gte(hit, 31)
    expect_lte(hit, 36)
})

test_that("a ramp reaches full size after rampRuns", {
    s <- simulateSpikeRuns(erccDesign(), nRuns = 12, driftAt = 5,
                           driftShape = "ramp", rampRuns = 4,
                           driftType = "offset", driftSize = 2, seed = 1)
    expect_equal(s$drifted, c(rep(FALSE, 4), rep(TRUE, 8)))
})
