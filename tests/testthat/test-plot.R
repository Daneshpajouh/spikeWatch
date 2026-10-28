test_that("plots build and draw without warnings", {
    set.seed(1)
    s <- simulateSpikeRuns(erccDesign(), nRuns = 30)
    expect_s3_class(plotDoseResponse(s, 1), "ggplot")
    m <- runMetrics(s)
    pd <- plotDoseResponse(s, 1, minAmount = attr(m, "minAmount"))
    expect_silent(ggplot2::ggplotGrob(pd))
    for (m in c("cusum", "ewma", "shewhart")) {
        ch <- driftChart(runMetrics(s), "slope", baseline = 20, method = m)
        p <- plotDrift(ch, "slope")
        expect_s3_class(p, "ggplot")
        expect_silent(ggplot2::ggplotGrob(p))
    }
})
