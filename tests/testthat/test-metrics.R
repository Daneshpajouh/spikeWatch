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
    # every spike seen: the limit is not estimated
    expect_true(is.na(m$lod))
    expect_equal(m$lodStatus, "all_detected")
    expect_equal(m$allDetectedFrom, 1)
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
    expect_false(isTRUE(all.equal(runMetrics(se, minAmount = 0)$slope, 1)))
    m <- runMetrics(se, minAmount = 2)
    expect_equal(m$slope, 1)
    expect_equal(m$offset, 1)
    expect_equal(m$allDetectedFrom, 1)
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
    expect_equal(m$lodStatus, "fit")
    expect_gte(m$lod, 8)
    expect_lte(m$lod, 512)
    expect_equal(m$allDetectedFrom, 32)
    # perfectly separated detections give no estimate, only the empirical
    # figure in its own column
    o2 <- c(0, 0, 0, 5, 9, 20, 40, 80, 160, 320)
    m2 <- runMetrics(powerLawSE(o2, e))
    expect_true(is.na(m2$lod))
    expect_equal(m2$lodStatus, "separated")
    expect_equal(m2$allDetectedFrom, 8)
    m3 <- runMetrics(powerLawSE(rep(0, 10), e), minAmount = 0)
    expect_true(is.na(m3$lod))
    expect_equal(m3$lodStatus, "none_detected")
    expect_true(is.na(m3$allDetectedFrom))
})

test_that("a limit outside the tested amounts is refused, not clamped", {
    # the fitted 95 percent point lies near 27, above the largest amount 8
    m <- runMetrics(powerLawSE(c(0, 5, 0, 5), c(1, 2, 4, 8)))
    expect_true(is.na(m$lod))
    expect_equal(m$lodStatus, "outside_range")
    expect_equal(m$allDetectedFrom, 8)
})

test_that("detection that falls with amount is a failed fit", {
    m <- runMetrics(powerLawSE(c(5, 0, 5, 0, 0), c(1, 2, 4, 8, 16)),
                    minAmount = 0)
    expect_true(is.na(m$lod))
    expect_equal(m$lodStatus, "fit_failed")
    expect_true(is.na(m$allDetectedFrom))
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

twoRunSE <- function(o1, o2, e) {
    ids <- paste0("s", seq_along(e))
    se <- SummarizedExperiment::SummarizedExperiment(
        list(counts = matrix(c(o1, o2), ncol = 2,
                             dimnames = list(ids, c("r1", "r2")))))
    addTruth(se, spikeDesign(ids, e))
}

test_that("minAmount auto takes the lowest amount every run detects", {
    e <- c(1, 2, 4, 8, 16)
    # run 1 misses amount 1, run 2 misses amount 2
    se <- twoRunSE(c(0, 4, 8, 16, 32), c(2, 0, 8, 16, 32), e)
    m <- runMetrics(se)
    expect_equal(attr(m, "minAmount"), 4)
    expect_equal(m$slope, c(1, 1))
    # only the named baseline run chooses the range
    m1 <- runMetrics(se, baseline = 1)
    expect_equal(attr(m1, "minAmount"), 2)
    # run 2 alone misses amount 2, so its range starts at 4
    expect_equal(attr(runMetrics(se, baseline = "r2"), "minAmount"), 4)
    expect_equal(attr(runMetrics(se, baseline = c(TRUE, FALSE)),
                      "minAmount"), 2)
    expect_equal(attr(runMetrics(se, baseline = "r1"), "minAmount"), 2)
    # a single number counts leading runs: both runs here
    expect_equal(attr(runMetrics(se, baseline = 2), "minAmount"), 4)
    # a number is kept as given
    expect_equal(attr(runMetrics(se, minAmount = 0), "minAmount"), 0)
})

test_that("minAmount auto gives NA fit figures when no amount qualifies", {
    e <- c(1, 2, 4)
    se <- twoRunSE(c(2, 4, 0), c(2, 4, 8), e)
    expect_warning(m <- runMetrics(se), "no expected amount")
    expect_true(is.na(attr(m, "minAmount")))
    expect_true(all(is.na(m$slope)) && all(is.na(m$offset)))
    expect_equal(m$nDetected, c(2, 3))
    expect_equal(attr(runMetrics(se, baseline = "r2"), "minAmount"), 1)
    expect_error(runMetrics(se, baseline = "r9"), "not in the data")
    expect_error(runMetrics(se, baseline = TRUE), "one TRUE or FALSE")
    expect_error(runMetrics(se, minAmount = -1), "minAmount")
    expect_error(runMetrics(se, minAmount = "low"), "minAmount")
    expect_error(runMetrics(se, baseline = 3), "outside")
})
