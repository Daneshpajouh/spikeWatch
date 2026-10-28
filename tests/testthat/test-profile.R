test_that("the profile slope is 1 in control and falls under compression", {
    set.seed(8)
    s <- simulateSpikeRuns(erccDesign(), nRuns = 30, driftAt = 21,
                           driftType = "slope", driftSize = -0.15)
    p <- profileMetrics(s, baseline = 20)
    expect_equal(nrow(p), 30)
    expect_lt(abs(mean(p$profileSlope[1:20]) - 1), 0.01)
    expect_lt(mean(p$profileSlope[21:30]), 0.92)
    expect_true(all(p$nSpikes == length(attr(p, "reference"))))
})

test_that("spikes and reference come from the baseline runs only", {
    m <- matrix(c(10, 100, 1000, 7,
                  12, 90, 1100, 9,
                  0, 80, 900, 8,
                  1, 1, 1, 1), nrow = 4,
                dimnames = list(c("a", "b", "c", "g"), paste0("r", 1:4)))
    se <- SummarizedExperiment::SummarizedExperiment(list(counts = m))
    SummarizedExperiment::rowData(se)$spike <- c(TRUE, TRUE, TRUE, FALSE)
    m4 <- rbind(m, d = c(50, 60, 70, 80))
    se4 <- SummarizedExperiment::SummarizedExperiment(list(counts = m4))
    SummarizedExperiment::rowData(se4)$spike <- c(TRUE, TRUE, TRUE, FALSE,
                                                  TRUE)
    p <- profileMetrics(se4, baseline = 2)
    expect_named(attr(p, "reference"), c("a", "b", "c", "d"))
    expect_equal(unname(attr(p, "reference")["a"]),
                 mean(log2(c(10, 12) + 0.5)))
    p3 <- profileMetrics(se4, baseline = 3)
    expect_named(attr(p3, "reference"), c("b", "c", "d"))
    expect_error(profileMetrics(se, baseline = 3), "three")
})

test_that("the z map flags the spike that moved", {
    set.seed(4)
    s <- simulateSpikeRuns(erccDesign(), nRuns = 22)
    cnt <- SummarizedExperiment::assay(s, "counts")
    top <- names(which.max(rowSums(cnt)))
    cnt[top, 22] <- cnt[top, 22] * 8
    SummarizedExperiment::assay(s, "counts") <- cnt
    z <- attr(profileMetrics(s, baseline = 20), "z")
    expect_equal(names(which.max(abs(z[, 22]))), top)
    expect_gt(abs(z[top, 22]), 5)
})

test_that("non-spike rows are left out", {
    set.seed(1)
    s <- simulateSpikeRuns(erccDesign(), nRuns = 5)
    p <- profileMetrics(s, baseline = 5)
    rd <- SummarizedExperiment::rowData(s)
    expect_true(all(names(attr(p, "reference")) %in%
                    rownames(s)[rd$spike]))
})

test_that("the z map matches a hand-worked example", {
    # log2 values: a (1, 5, 1), b (1, 3, 9), c (4, 2, 4); reference over
    # runs 1 and 2 is (3, 2, 3); minus each run's median shift (-1, 1, 1)
    # gives a (-1, 1, -3), b (0, 0, 6), c (2, -2, 0); baseline SDs are
    # sqrt(2), 0 (so NA) and 2 * sqrt(2).
    m <- matrix(c(2, 2, 16, 32, 8, 4, 2, 512, 16), nrow = 3,
                dimnames = list(c("a", "b", "c"), paste0("r", 1:3)))
    se <- SummarizedExperiment::SummarizedExperiment(list(counts = m))
    z <- attr(profileMetrics(se, baseline = 2, pseudo = 0), "z")
    expect_equal(unname(z[, 3]), c(-3 / sqrt(2), NA, 0))
    expect_equal(unname(z[, 1]), c(-1 / sqrt(2), NA, 2 / (2 * sqrt(2))))
})
