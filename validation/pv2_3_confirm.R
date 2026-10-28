# Plan v2, confirmation: C1 and C2 for F0 and F5 on seeds 900001 to
# 900400, with thresholds from the calibration seeds. Run once.
# Usage: Rscript validation/pv2_3_confirm.R SERIESDIR OUTDIR
suppressMessages(library(spikeWatch))
args <- commandArgs(trailingOnly = TRUE)
sdir <- args[1]
outdir <- args[2]
src <- readLines(file.path(dirname(sub("--file=", "", grep("--file=",
    commandArgs(FALSE), value = TRUE))), "pv2_2_score.R"))
keep <- seq_len(grep("^scoreFamily <- ", src) - 1)
eval(parse(text = src[keep]))
thr <- utils::read.csv(file.path(outdir, "thresholds.csv"))
set.seed(20261099)
rows <- list()
for (fam in c("F0", "F5")) {
    ser <- readRDS(file.path(sdir, sprintf("series_confirm_%s.rds", fam)))
    for (key in names(ser)) {
        st <- statsOf(ser[[key]])
        for (pn in c("C1", "C2")) for (chart in charts) for (H in Hs) {
            a <- pairs[[pn]][1]
            b <- pairs[[pn]][2]
            th <- function(d) thr$threshold[thr$family == fam &
                thr$detector == d & thr$chart == chart & thr$H == H]
            ha <- rowMaxH(detStat(st, a, chart), H) > th(a)
            hb <- rowMaxH(detStat(st, b, chart), H) > th(b)
            d <- as.numeric(ha) - as.numeric(hb)
            n <- length(d)
            idx <- matrix(sample.int(n, n * 2000, replace = TRUE), 2000)
            q <- stats::quantile(rowMeans(matrix(d[idx], 2000)),
                                 c(0.025, 0.975), names = FALSE)
            rows[[length(rows) + 1]] <- data.frame(
                family = fam, pair = pn, chart = chart, H = H,
                scenario = sub("confirm:", "", key), n = n,
                detect_truth = mean(ha), detect_without = mean(hb),
                diff = mean(d), ci_low = q[1], ci_high = q[2],
                verdict = if (q[1] > 0) "gain" else if (q[2] < 0) "loss"
                    else "no clear difference")
        }
    }
}
utils::write.csv(do.call(rbind, rows), file.path(outdir, "confirm.csv"),
                 row.names = FALSE)
print(sessionInfo())
