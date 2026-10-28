powerLawSE <- function(o, e) {
    se <- SummarizedExperiment::SummarizedExperiment(
        list(counts = matrix(o, ncol = 1,
                             dimnames = list(paste0("s", seq_along(o)),
                                             "r1"))))
    addTruth(se, spikeDesign(paste0("s", seq_along(e)), e))
}

test_that("runMetrics recovers an exact power law", {
    e <- c(1, 2, 4, 8)
    se <- powerLawSE(4 * e^0.8, e)
    m <- runMetrics(se)
    expect_equal(m$slope, 0.8)
    expect_equal(m$intercept, 2)
    expect_equal(m$r2, 1)
    expect_equal(m$residSD, 0, tolerance = 1e-12)
    # median of 2 - 0.2 * c(0, 1, 2, 3)
    expect_equal(m$offset, 1.7)
    expect_equal(m$lod, 1)
    expect_true(m$lodCensored)
    expect_equal(m$nDetected, 4)
})

test_that("runMetrics leaves the slope empty when amounts do not vary", {
    se <- powerLawSE(c(10, 12, 8), c(5, 5, 5))
    m <- runMetrics(se)
    expect_true(is.na(m$slope))
    expect_equal(m$offset, log2(10 / 5))
})

test_that("minAmount keeps low spikes out of the fit but not the limit", {
    e <- c(1, 2, 4, 8, 16)
    o <- c(9, 2 * e[2:5])
    se <- powerLawSE(o, e)
    expect_false(isTRUE(all.equal(runMetrics(se)$slope, 1)))
    m <- runMetrics(se, minAmount = 2)
    expect_equal(m$slope, 1)
    expect_equal(m$offset, 1)
    expect_equal(m$lod, 1)
    # a spike exactly at minAmount is in the fit
    se2 <- powerLawSE(c(9, 4, 8, 16), c(1, 2, 4, 8))
    m2 <- runMetrics(se2, minAmount = 2)
    expect_equal(m2$slope, 1)
    expect_equal(m2$offset, 1)
})

test_that("the detection limit sits between missed and found spikes", {
    e <- 2^(0:9)
    o <- c(0, 0, 0, 1, 0, 6, 12, 30, 50, 100)
    m <- runMetrics(powerLawSE(o, e), lodProb = 0.95)
    expect_false(m$lodCensored)
    expect_gte(m$lod, 8)
    expect_lte(m$lod, 512)
    # perfectly separated detections fall back to the empirical limit
    o2 <- c(0, 0, 0, 5, 9, 20, 40, 80, 160, 320)
    m2 <- runMetrics(powerLawSE(o2, e))
    expect_equal(m2$lod, 8)
    expect_true(is.na(runMetrics(powerLawSE(rep(0, 10), e))$lod))
})

test_that("ratioAccuracy matches hand-computed errors", {
    # background rows read twice as high in x, so centring on them must
    # remove a shift of 1 on the log2 scale
    m <- matrix(c(80, 20, 20, 20, 20, 20,
                  10, 10, 20, 10, 10, 10), nrow = 6,
                dimnames = list(c("u", "s", "d", "b1", "b2", "b3"),
                                c("x", "y")))
    se <- SummarizedExperiment::SummarizedExperiment(
        list(counts = m),
        colData = S4Vectors::DataFrame(mix = c("A", "B"),
                                       row.names = c("x", "y")))
    d <- spikeDesign(rep(c("u", "s", "d"), 2),
                     amount = c(4, 1, 1, 1, 1, 2),
                     mix = rep(c("A", "B"), each = 3))
    se <- addTruth(se, d, mix = "mix")
    r <- ratioAccuracy(se, "x", "y", center = "none")
    expect_equal(r$features$error, c(1, 1, 1))
    expect_equal(ratioAccuracy(se, 1, 2)$features$error, c(0, 0, 0))
    bg <- rownames(se) %in% c("b1", "b2", "b3")
    r2 <- ratioAccuracy(se, 1, 2, center = bg)
    expect_equal(r2$features$measuredLog2Ratio, c(2, 0, -1))
    expect_equal(r2$summary$expectedLog2Ratio, c(2, 0, -1))
    expect_error(ratioAccuracy(se, 1, 2, center = "median"), "center")
})
