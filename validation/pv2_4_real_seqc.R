# Plan v2 rerun of V2 at the v2 code: minAmount "auto" over every run.
# V2: real RNA-seq, SEQC lab 5 ERCC counts from erccdashboard. PLAN.md.
suppressMessages({
    library(spikeWatch)
    library(SummarizedExperiment)
    library(erccdashboard)
})
data(SEQC.Example)
data(ERCC)
cnt <- UHRR.HBRR.countDat
rownames(cnt) <- cnt$Feature
cnt <- as.matrix(cnt[grepl("^ERCC-", cnt$Feature), -1])
mix <- ifelse(grepl("^UHRR", colnames(cnt)), "Mix1", "Mix2")
cd <- S4Vectors::DataFrame(mix = mix, row.names = colnames(cnt))
design <- erccDesign()
len <- ERCCDef$Length[match(rownames(cnt), ERCCDef$Feature)]
stopifnot(!anyNA(len))
se <- SummarizedExperiment(list(counts = cnt, cpk = cnt / len * 1000),
                           colData = cd)
se <- addTruth(se, design, mix = "mix")
cat("ERCC rows:", nrow(se), " in design:", sum(rowData(se)$spike), "\n\n")

m <- runMetrics(se, assay = "cpk", minAmount = "auto")
cat("minAmount chosen:", attr(m, "minAmount"), "\n")
dir.create("validation/out", showWarnings = FALSE)
write.csv(m, "validation/out-v2/real_seqc_run_metrics.csv", row.names = FALSE)
print(m, digits = 3, row.names = FALSE)

r <- ratioAccuracy(se, mix == "Mix1", mix == "Mix2", assay = "counts",
                   center = "spikes")
write.csv(r$features, "validation/out-v2/real_seqc_ratio_features.csv",
          row.names = FALSE)
write.csv(r$summary, "validation/out-v2/real_seqc_ratio_summary.csv",
          row.names = FALSE)
cat("\nRatio accuracy, UHRR (Mix1) over HBRR (Mix2), centred on spikes:\n")
print(r$summary, digits = 3, row.names = FALSE)
