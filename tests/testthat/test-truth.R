test_that("addTruth fills design amount times dose, NA elsewhere", {
    m <- matrix(1, 4, 2, dimnames = list(c("s1", "s2", "s3", "g"),
                                         c("x", "y")))
    d <- spikeDesign(c("s1", "s2", "s3", "s1", "s2", "s3"),
                     amount = c(1, 10, 100, 2, 10, 50),
                     mix = rep(c("A", "B"), each = 3),
                     group = c("up", "same", "down"))
    se <- addTruth(tinySE(m), d, mix = "mix", dose = "dose")
    e <- SummarizedExperiment::assay(se, "expected")
    expect_equal(unname(e[, "x"]), c(1, 10, 100, NA))
    expect_equal(unname(e[, "y"]), c(4, 20, 100, NA))
    rd <- SummarizedExperiment::rowData(se)
    expect_equal(rd$spike, c(TRUE, TRUE, TRUE, FALSE))
    expect_equal(rd$spikeGroup, c("up", "same", "down", NA))
})

test_that("addTruth refuses inputs it cannot match", {
    m <- matrix(1, 2, 2, dimnames = list(c("s1", "g"), c("x", "y")))
    d2 <- spikeDesign(c("s1", "s1"), c(1, 2), mix = c("A", "B"))
    expect_error(addTruth(tinySE(m), d2), "2 mixes")
    expect_error(addTruth(tinySE(m), spikeDesign("zz", 1)), "no design")
    expect_error(addTruth(tinySE(m), d2, mix = "nope"), "no column")
    d3 <- spikeDesign(c("s1"), 1, mix = "C")
    expect_error(addTruth(tinySE(m), d3, mix = "mix"), "not in the design")
    expect_error(addTruth(m, d3), "SummarizedExperiment")
})

test_that("addTruth warns about design features it cannot place", {
    m <- matrix(1, 2, 2, dimnames = list(c("s1", "g"), c("x", "y")))
    d <- spikeDesign(c("s1", "s2", "s3"), c(1, 2, 4))
    expect_warning(se <- addTruth(tinySE(m), d), "2 design feature")
    e <- SummarizedExperiment::assay(se, "expected")
    expect_equal(unname(e["s1", ]), c(1, 1))
})
