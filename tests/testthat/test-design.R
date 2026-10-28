test_that("spikeDesign checks its input", {
    expect_error(spikeDesign(character(), numeric()), "at least one")
    expect_error(spikeDesign("a", -1), "positive")
    expect_error(spikeDesign(c("a", "b"), 1), "one value per feature")
    expect_error(spikeDesign(c("a", "a"), c(1, 2)), "once per mix")
    d <- spikeDesign(c("a", "a"), c(1, 2), mix = c("M1", "M2"),
                     unit = "fmol")
    expect_equal(nrow(d), 2)
    expect_equal(attr(d, "unit"), "fmol")
})

test_that("erccDesign has 92 spikes with the published subgroup ratios", {
    d <- erccDesign()
    expect_equal(length(unique(d$feature)), 92)
    expect_setequal(unique(d$mix), c("Mix1", "Mix2"))
    m1 <- d[d$mix == "Mix1", ]
    m2 <- d[d$mix == "Mix2", ]
    m2 <- m2[match(m1$feature, m2$feature), ]
    ratio <- tapply(m1$amount / m2$amount, m1$group,
                    function(r) unique(round(r, 4)))
    expect_equal(as.numeric(ratio[c("4:1", "1:1", "1:1.5", "1:2")]),
                 c(4, 1, 0.6667, 0.5))
    expect_equal(as.vector(table(m1$group)), c(23, 23, 23, 23))
})
