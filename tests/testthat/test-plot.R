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

test_that("the dose-response line is the runMetrics fit on the log axes", {
    set.seed(3)
    s <- simulateSpikeRuns(erccDesign(), nRuns = 2)
    m <- runMetrics(s, baseline = 2)
    p <- plotDoseResponse(s, 1, minAmount = attr(m, "minAmount"))
    line <- ggplot2::layer_data(p, 2)
    # Scales transform before the smoother runs, so the line is fitted on
    # log10 axes: its slope there is the log2 slope of runMetrics.
    ends <- line[order(line$x), ][c(1, nrow(line)), ]
    slope <- diff(ends$y) / diff(ends$x)
    expect_equal(slope, m$slope[1], tolerance = 1e-8)
})
