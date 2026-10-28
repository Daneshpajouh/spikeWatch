test_that("plots build", {
    set.seed(1)
    s <- simulateSpikeRuns(erccDesign(), nRuns = 12)
    expect_s3_class(plotDoseResponse(s, 1), "ggplot")
    ch <- driftChart(runMetrics(s), "slope", baseline = 8)
    p <- plotDrift(ch, "slope")
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
})
