# V0: why runMetrics() takes minAmount. 200 simulated runs, no drift,
# ERCC Mix 1. Compares the slope fitted on every detected spike with the
# slope fitted on spikes chosen by known amount.
suppressMessages(library(spikeWatch))
s <- simulateSpikeRuns(erccDesign(), nRuns = 200, seed = 11)
out <- do.call(rbind, lapply(c(0, 1, 10), function(a) {
    m <- runMetrics(s, minAmount = a)
    data.frame(minAmount = a, meanSlope = mean(m$slope),
               sdSlope = stats::sd(m$slope),
               corWithDepth = stats::cor(m$slope, s$depth))
}))
dir.create("validation/out", showWarnings = FALSE)
write.csv(out, "validation/out/v0_min_amount.csv", row.names = FALSE)
print(out, digits = 3, row.names = FALSE)
