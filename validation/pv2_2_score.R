# Plan v2, step 2: calibrate, score and compare every detector.
# Rules in PLAN_v2.md. Usage: Rscript validation/pv2_2_score.R SERIESDIR OUTDIR
suppressMessages(library(spikeWatch))
args <- commandArgs(trailingOnly = TRUE)
sdir <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
fams <- paste0("F", 0:5)
figs <- c("slope", "residSD", "offset", "B1", "B2", "B3", "B4")
dets <- list(slope = "slope", residSD = "residSD", offset = "offset",
             B1 = "B1", B2 = "B2", B3 = "B3",
             UT = c("slope", "residSD"), UF = c("B2", "B3"),
             B4 = "B4", UF2 = c("B3", "B4"))
charts <- c("cusum", "ewma")
Hs <- c(10, 20)
defaults <- c(cusum = 5, ewma = 3)
drifts <- c("compression_small", "compression", "compression_ramp",
            "sensitivity_loss", "precision_loss")
pairs <- list(C1 = c("UT", "UF"), C2 = c("slope", "B3"),
              C3 = c("residSD", "B2"), C4 = c("offset", "B1"),
              C5 = c("residSD", "B4"), C6 = c("UT", "UF2"))

# Standardized statistic over the 20 monitored runs: CUSUM sums in sigma
# units, EWMA distance over its exact standard deviation. A missing run
# cannot signal, so its statistic is NA.
monitorStat <- function(x, chart) {
    ch <- tryCatch(suppressWarnings(
        driftChart(x, baseline = 20, method = chart, L = 1)),
        error = function(e) NULL)
    if (is.null(ch))
        return(rep(NA_real_, 20))
    st <- if (chart == "cusum") ch$stat else
        abs(ch$stat - ch$center) / (ch$upper - ch$center)
    st[is.na(x)] <- NA_real_
    st[21:40]
}

# stats[[fig]][[chart]] is a series by run matrix.
statsOf <- function(arr) {
    lapply(stats::setNames(figs, figs), function(f)
        lapply(stats::setNames(charts, charts), function(ch)
            t(apply(arr[, , f], 1, monitorStat, chart = ch))))
}

detStat <- function(st, det, chart) {
    m <- st[[dets[[det]][1]]][[chart]]
    if (length(dets[[det]]) > 1)
        m <- pmax(m, st[[dets[[det]][2]]][[chart]], na.rm = TRUE)
    m
}

rowMaxH <- function(m, H) {
    v <- apply(m[, seq_len(H), drop = FALSE], 1, function(r)
        if (all(is.na(r))) -Inf else max(r, na.rm = TRUE))
    v
}

firstOver <- function(m, H, thr) {
    apply(m[, seq_len(H), drop = FALSE], 1, function(r) {
        w <- which(!is.na(r) & r > thr)
        if (length(w)) w[1] else H + 1
    })
}

cp <- function(k, n) {
    if (n == 0) return(c(NA, NA))
    stats::binom.test(k, n)$conf.int[1:2]
}

scoreFamily <- function(fam) {
    ser <- readRDS(file.path(sdir, sprintf("series_main_%s.rds", fam)))
    st <- lapply(ser, statsOf)
    thr <- list()
    rows <- list()
    for (det in names(dets)) for (chart in charts) for (H in Hs) {
        mc <- rowMaxH(detStat(st[["calib:null"]], det, chart), H)
        cal <- sort(mc)[950]
        thr[[length(thr) + 1]] <- data.frame(family = fam, detector = det,
            chart = chart, H = H, threshold = cal)
        for (lim in c("calibrated", "default")) {
            th <- if (lim == "calibrated") cal else defaults[[chart]]
            for (key in c("eval:null", paste0("eval:", drifts))) {
                m <- detStat(st[[key]], det, chart)
                hit <- rowMaxH(m, H) > th
                fo <- firstOver(m, H, th)
                ci <- cp(sum(hit), length(hit))
                rows[[length(rows) + 1]] <- data.frame(
                    family = fam, detector = det, chart = chart, H = H,
                    limits = lim, threshold = th,
                    scenario = sub("eval:", "", key), n = length(hit),
                    signalled = sum(hit), rate = mean(hit),
                    ci_low = ci[1], ci_high = ci[2],
                    delay_mean_censored = mean(fo),
                    delay_median_detected = if (any(hit))
                        stats::median(fo[hit]) else NA)
            }
        }
    }
    thrs <- do.call(rbind, thr)
    res <- do.call(rbind, rows)

    # Paired bootstrap of detection differences at calibrated limits.
    set.seed(20261028 + match(fam, fams))
    cmp <- list()
    for (pn in names(pairs)) for (chart in charts) for (H in Hs) {
        a <- pairs[[pn]][1]
        b <- pairs[[pn]][2]
        ta <- thrs$threshold[thrs$detector == a & thrs$chart == chart &
                             thrs$H == H]
        tb <- thrs$threshold[thrs$detector == b & thrs$chart == chart &
                             thrs$H == H]
        fa <- res[res$detector == a & res$chart == chart & res$H == H &
                  res$limits == "calibrated" & res$scenario == "null", ]
        fb <- res[res$detector == b & res$chart == chart & res$H == H &
                  res$limits == "calibrated" & res$scenario == "null", ]
        for (dr in drifts) {
            key <- paste0("eval:", dr)
            ha <- rowMaxH(detStat(st[[key]], a, chart), H) > ta
            hb <- rowMaxH(detStat(st[[key]], b, chart), H) > tb
            d <- as.numeric(ha) - as.numeric(hb)
            n <- length(d)
            idx <- matrix(sample.int(n, n * 2000, replace = TRUE), 2000)
            bs <- rowMeans(matrix(d[idx], 2000))
            q <- stats::quantile(bs, c(0.025, 0.975), names = FALSE)
            verdict <- if (q[1] > 0) "gain" else if (q[2] < 0) "loss" else
                "no clear difference"
            cmp[[length(cmp) + 1]] <- data.frame(
                family = fam, pair = pn, truth = a, without = b,
                chart = chart, H = H, scenario = dr, n = n,
                detect_truth = mean(ha), detect_without = mean(hb),
                diff = mean(d), ci_low = q[1], ci_high = q[2],
                verdict = verdict,
                fap_truth = fa$rate, fap_without = fb$rate,
                fap_flag = (fa$ci_low > 0.05 | fa$ci_high < 0.05) |
                    (fb$ci_low > 0.05 | fb$ci_high < 0.05))
        }
    }
    cmps <- do.call(rbind, cmp)

    # calibrateChart() on each series' own baseline; it is not in this
    # release (addendum 2), so the block runs only where it exists.
    cc <- list()
    if (exists("calibrateChart", mode = "function"))
    for (key in c("eval:null", paste0("eval:", drifts))) {
        arr <- ser[[key]]
        seeds <- as.integer(dimnames(arr)[[1]])
        hit <- array(NA, c(length(seeds), 4, 2, 2),
                     dimnames = list(NULL, c("slope", "residSD", "offset",
                                             "B3"), charts, Hs))
        for (i in seq_along(seeds)) {
            set.seed(seeds[i] + 5000000)
            for (f in dimnames(hit)[[2]]) for (chart in charts)
                for (H in Hs) {
                    x <- arr[i, , f]
                    th <- tryCatch(calibrateChart(x, baseline = 20,
                        method = chart, alpha = 0.05, horizon = H, B = 500),
                        error = function(e) NA_real_)
                    if (is.na(th)) next
                    s <- st[[key]][[f]][[chart]][i, seq_len(H)]
                    hit[i, f, chart, as.character(H)] <-
                        any(!is.na(s) & s > th)
                }
        }
        for (f in dimnames(hit)[[2]]) for (chart in charts) for (H in Hs) {
            h <- hit[, f, chart, as.character(H)]
            ok <- !is.na(h)
            ci <- cp(sum(h[ok]), sum(ok))
            cc[[length(cc) + 1]] <- data.frame(
                family = fam, figure = f, chart = chart, H = H,
                scenario = sub("eval:", "", key), n = length(h),
                failed = sum(!ok), signalled = sum(h[ok]),
                rate = mean(h[ok]), ci_low = ci[1], ci_high = ci[2])
        }
    }
    list(thresholds = thrs, scores = res, comparisons = cmps,
         calibrate = do.call(rbind, cc))
}

t0 <- Sys.time()
out <- parallel::mclapply(fams, scoreFamily, mc.cores = 3)
bad <- vapply(out, inherits, logical(1), "try-error")
if (any(bad)) {
    print(out[bad])
    stop("scoring failed for ", paste(fams[bad], collapse = ", "))
}
for (part in c("thresholds", "scores", "comparisons", "calibrate"))
    if (!is.null(out[[1]][[part]])) utils::write.csv(do.call(rbind, lapply(out, `[[`, part)),
                     file.path(outdir, paste0(part, ".csv")),
                     row.names = FALSE)

# Headline rule: a pair, drift and chart with the same verdict in at least
# five of six families at both horizons.
cmps <- do.call(rbind, lapply(out, `[[`, "comparisons"))
hl <- do.call(rbind, lapply(split(cmps, list(cmps$pair, cmps$scenario,
                                             cmps$chart), drop = TRUE),
    function(d) {
        both <- tapply(d$verdict, d$family, function(v)
            if (all(v == "gain")) "gain" else if (all(v == "loss")) "loss"
            else "mixed")
        data.frame(pair = d$pair[1], scenario = d$scenario[1],
                   chart = d$chart[1],
                   families_gain_both_H = sum(both == "gain"),
                   families_loss_both_H = sum(both == "loss"),
                   headline = if (sum(both == "gain") >= 5) "gain" else
                       if (sum(both == "loss") >= 5) "loss" else "none")
    }))
utils::write.csv(hl, file.path(outdir, "headline.csv"), row.names = FALSE)
cat("elapsed", format(Sys.time() - t0), "\n")
print(sessionInfo())
