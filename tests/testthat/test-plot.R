test_that("plots build and draw without warnings", {
    set.seed(1)
    s <- simulateSpikeRuns(erccDesign(), nRuns = 30)
    expect_s3_class(plotDoseResponse(s, 1), "ggplot")
    for (m in c("cusum", "ewma", "shewhart")) {
        ch <- driftChart(runMetrics(s), "slope", baseline = 20, method = m)
        p <- plotDrift(ch, "slope")
        expect_s3_class(p, "ggplot")
        expect_silent(ggplot2::ggplotGrob(p))
    }
})
