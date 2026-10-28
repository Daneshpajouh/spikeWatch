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

test_that("a run with nothing detected stops, and a named run is the title", {
    set.seed(1)
    s <- simulateSpikeRuns(erccDesign(), nRuns = 2)
    expect_error(plotDoseResponse(s, 1, detect = Inf), "no spike is above detect")
    name <- colnames(s)[2]
    expect_equal(plotDoseResponse(s, name)$labels$title, name)
    expect_equal(plotDoseResponse(s, 2)$labels$title, name)
    unnamed <- s
    colnames(unnamed) <- NULL
    expect_equal(plotDoseResponse(unnamed, 2)$labels$title, "2")
})

test_that("by default the line uses the auto range runMetrics chooses", {
    set.seed(3)
    s <- simulateSpikeRuns(erccDesign(), nRuns = 4)
    slopeOf <- function(p) {
        line <- ggplot2::layer_data(p, 2)
        ends <- line[order(line$x), ][c(1, nrow(line)), ]
        diff(ends$y) / diff(ends$x)
    }
    m <- runMetrics(s)
    p <- plotDoseResponse(s, 1)
    expect_equal(ggplot2::layer_data(p, 3)$xintercept, log10(attr(m, "minAmount")))
    expect_equal(slopeOf(p), m$slope[1], tolerance = 1e-8)
    m2 <- runMetrics(s, baseline = 2)
    expect_equal(slopeOf(plotDoseResponse(s, 1, baseline = 2)), m2$slope[1], tolerance = 1e-8)
    expect_error(plotDoseResponse(s, 1, minAmount = -1), "must be")
    expect_error(plotDoseResponse(s, 1, minAmount = "all"), "must be")
})

test_that("an auto range that no amount meets warns and draws no line", {
    set.seed(3)
    s <- simulateSpikeRuns(erccDesign(), nRuns = 2)
    i <- which.max(SummarizedExperiment::assay(s, "expected")[, 2])
    a <- SummarizedExperiment::assay(s, 1)
    a[i, 2] <- 0
    SummarizedExperiment::assay(s, 1) <- a
    expect_warning(p <- plotDoseResponse(s, 1), "no fit line is drawn")
    expect_equal(nrow(p$layers[[2]]$data), 0L)
    expect_warning(runMetrics(s), "fit figures are NA")
})
