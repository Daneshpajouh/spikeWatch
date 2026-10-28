# Plan v2, step 1: simulate every series and keep its per-run figures.
# Rules in PLAN_v2.md. Usage: Rscript validation/pv2_1_series.R OUTDIR [SET]
# SET is "main" (calibration, evaluation null and drift seeds) or
# "confirm" (the confirmation seeds, run once after the main results).
suppressMessages(library(spikeWatch))
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
set <- if (length(args) > 1) args[2] else "main"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
design <- erccDesign()

families <- list(
    F0 = list(),
    F1 = list(pivot = "top"),
    F2 = list(counts = "nbinom", size = 10),
    F3 = list(saturation = 0.5),
    F4 = list(spikeBiasSD = 0.5),
    F5 = list(runCor = 0.6, slopeSD = 0.02))
drifts <- list(
    null = list(type = "slope", size = 0, at = NA, shape = "step"),
    compression_small = list(type = "slope", size = -0.05, at = 21,
                             shape = "step"),
    compression = list(type = "slope", size = -0.10, at = 21,
                       shape = "step"),
    compression_ramp = list(type = "slope", size = -0.10, at = 21,
                            shape = "ramp"),
    sensitivity_loss = list(type = "offset", size = -1, at = 21,
                            shape = "step"),
    precision_loss = list(type = "noise", size = 0.25, at = 21,
                          shape = "step"))

jobs <- if (set == "main") list(
    list(role = "calib", drift = "null", seeds = 100001:101000),
    list(role = "eval", drift = "null", seeds = 200001:201000),
    list(role = "eval", drift = "compression_small", seeds = 300001:300400),
    list(role = "eval", drift = "compression", seeds = 300001:300400),
    list(role = "eval", drift = "compression_ramp", seeds = 300001:300400),
    list(role = "eval", drift = "sensitivity_loss", seeds = 300001:300400),
    list(role = "eval", drift = "precision_loss", seeds = 300001:300400)) else
    list(list(role = "confirm", drift = "compression", seeds = 900001:900400),
         list(role = "confirm", drift = "sensitivity_loss",
              seeds = 900001:900400))
famNames <- if (set == "main") names(families) else c("F0", "F5")

figures <- function(se) {
    m <- suppressWarnings(runMetrics(se, minAmount = "auto", baseline = 20))
    cnt <- SummarizedExperiment::assay(se, "counts")
    spike <- SummarizedExperiment::rowData(se)$spike
    sc <- cnt[spike, , drop = FALSE]
    b1 <- unname(apply(sc, 2, function(v) log2(sum(v))))
    b2 <- unname(apply(sc, 2, function(v) stats::sd(log2(v[v > 0]))))
    p <- profileMetrics(se, baseline = 20)
    cbind(slope = m$slope, residSD = m$residSD, offset = m$offset,
          B1 = b1, B2 = b2, B3 = p$profileSlope)
}

runFamily <- function(fam) {
    out <- list()
    for (jb in jobs) {
        dr <- drifts[[jb$drift]]
        arr <- array(NA_real_, c(length(jb$seeds), 40, 6),
                     dimnames = list(jb$seeds, NULL,
                                     c("slope", "residSD", "offset", "B1",
                                       "B2", "B3")))
        for (i in seq_along(jb$seeds)) {
            set.seed(jb$seeds[i])
            se <- do.call(simulateSpikeRuns, c(list(
                design, nRuns = 40, mix = "Mix1", driftAt = dr$at,
                driftType = dr$type, driftSize = dr$size,
                driftShape = dr$shape, rampRuns = 10), families[[fam]]))
            arr[i, , ] <- figures(se)
        }
        out[[paste(jb$role, jb$drift, sep = ":")]] <- arr
    }
    saveRDS(out, file.path(outdir, sprintf("series_%s_%s.rds", set, fam)))
    fam
}

t0 <- Sys.time()
done <- parallel::mclapply(famNames, runFamily, mc.cores = 3)
print(unlist(done))
cat("elapsed", format(Sys.time() - t0), "\n")
print(sessionInfo())
