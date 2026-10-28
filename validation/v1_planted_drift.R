# V1: planted drift in simulated series. Rules in PLAN.md.
suppressMessages(library(spikeWatch))
design <- erccDesign()
scen <- list(
    null = list(type = "slope", size = 0, at = NA),
    compression_small = list(type = "slope", size = -0.05, at = 31),
    compression = list(type = "slope", size = -0.10, at = 31),
    sensitivity_loss = list(type = "offset", size = -1, at = 31),
    precision_loss = list(type = "noise", size = 0.25, at = 31))

noTruth <- function(se) {
    cnt <- SummarizedExperiment::assay(se, "counts")
    spike <- SummarizedExperiment::rowData(se)$spike
    b1 <- unname(apply(cnt[spike, ], 2, function(v) log2(sum(v))))
    b2 <- apply(cnt[spike, ], 2, function(v) stats::sd(log2(v[v > 0])))
    list(B1_total = b1, B2_spread = b2)
}

firstAfter <- function(x) {
    ch <- driftChart(x, baseline = 20, method = "cusum")
    s <- which(ch$signal)
    list(fa = any(s >= 21 & s <= 30), faAll = any(s >= 21),
         first = if (any(s >= 31)) min(s[s >= 31]) else NA)
}

rows <- list()
for (nm in names(scen)) {
    sc <- scen[[nm]]
    for (seed in 1:200) {
        set.seed(seed)
        se <- simulateSpikeRuns(design, nRuns = 40, mix = "Mix1",
                                driftAt = sc$at, driftType = sc$type,
                                driftSize = sc$size)
        m <- runMetrics(se, minAmount = 1)
        nt <- noTruth(se)
        series <- list(slope = m$slope, residSD = m$residSD,
                       offset = m$offset, B1_total = nt$B1_total,
                       B2_spread = nt$B2_spread)
        res <- lapply(series, firstAfter)
        u <- list(fa = res$slope$fa || res$residSD$fa,
                  faAll = res$slope$faAll || res$residSD$faAll,
                  first = suppressWarnings(min(c(res$slope$first,
                                                 res$residSD$first),
                                               na.rm = TRUE)))
        if (!is.finite(u$first)) u$first <- NA
        res$union_slope_residSD <- u
        for (d in names(res))
            rows[[length(rows) + 1]] <- data.frame(
                scenario = nm, seed = seed, detector = d,
                falseAlarm = res[[d]]$fa, falseAlarmAll = res[[d]]$faAll,
                firstSignal = res[[d]]$first)
    }
}
raw <- do.call(rbind, rows)
dir.create("validation/out", showWarnings = FALSE)
write.csv(raw, "validation/out/v1_raw.csv", row.names = FALSE)

summ <- do.call(rbind, lapply(split(raw, list(raw$scenario, raw$detector)),
    function(d) {
        clean <- !d$falseAlarm
        det <- clean & !is.na(d$firstSignal)
        data.frame(scenario = d$scenario[1], detector = d$detector[1],
                   n = nrow(d),
                   falseAlarm21to30 = mean(d$falseAlarm),
                   falseAlarm21to40 = mean(d$falseAlarmAll),
                   detectedBy40 = sum(det) / sum(clean),
                   medianDelay = if (any(det))
                       stats::median(d$firstSignal[det] - 31) else NA)
    }))
summ <- summ[order(summ$scenario, summ$detector), ]
write.csv(summ, "validation/out/v1_summary.csv", row.names = FALSE)
print(summ, row.names = FALSE)
cat("R", R.version.string, "spikeWatch", as.character(packageVersion("spikeWatch")), "\n")
