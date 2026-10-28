#' Accuracy figures for each run against the known spike-in truth
#'
#' For each sample, compares the measured spike-in values with their
#' expected amounts on the log2 scale.
#'
#' * `slope`, `intercept`, `r2`: least-squares fit of log2 measured on
#'   log2 expected, over detected spikes. A slope of 1 means the measured
#'   values follow the known amounts in proportion; below 1 means
#'   compression. Needs at least three detected spikes and more than one
#'   distinct expected amount, otherwise `NA`. Only spikes with an expected
#'   amount of at least `minAmount` enter the fit (see below).
#' * `residSD`: residual standard error of the fit (divisor n - 2), a
#'   precision figure.
#' * `offset`: median of log2 measured minus log2 expected, over the same
#'   spikes as the fit. It moves with sequencing depth or instrument
#'   response, and is defined even when all spikes in a run share one
#'   amount.
#' * `lod`: the expected amount at which a spike is detected with
#'   probability `lodProb`, from a logistic fit of detection on log2
#'   expected amount. It is given only when the data can support it, and
#'   `lodStatus` says why when it is `NA`:
#'   `"fit"` (the fit converged with a positive slope and the point lies
#'   inside the tested amounts), `"all_detected"` (every spike was seen, so
#'   the limit is at or below the smallest amount but not estimated),
#'   `"none_detected"`, `"separated"` (detected and missed spikes do not
#'   overlap, so the fit has no finite estimate), `"fit_failed"` (no
#'   convergence, or a slope that is not positive) and `"outside_range"`
#'   (the fitted point lies outside the tested amounts; it is never
#'   clamped).
#' * `allDetectedFrom`: the smallest expected amount at and above which
#'   every spike was detected, an empirical figure reported for every run
#'   and never a substitute for `lod`.
#'
#' @param se A `SummarizedExperiment` with an expected-amount assay from
#'   [addTruth()].
#' @param assay Name or index of the measured assay.
#' @param expected Name of the expected-amount assay.
#' @param detect A spike counts as detected when its measured value is
#'   above this threshold and not `NA`.
#' @param lodProb Detection probability that defines `lod`.
#' @param minAmount Smallest expected amount used for `slope`,
#'   `intercept`, `r2`, `residSD` and `offset`. Spikes near the detection
#'   limit are seen only when their count happens to be high, which pulls
#'   the slope down and ties it to depth. Choosing spikes by their known
#'   amount, not by what was measured, avoids that. `"auto"` (the default)
#'   takes the smallest expected amount at and above which every spike is
#'   detected in every run named by `baseline`; a number fixes it. When no
#'   amount qualifies, the fit figures are `NA` with a warning rather than
#'   fitted on a biased range. `lod` always uses all spikes.
#' @param baseline Runs used to choose `minAmount = "auto"`: a single
#'   number is the count of leading runs; otherwise run indices, run names
#'   or a logical vector. `NULL` (the default) uses every run, which suits
#'   one experiment scored once, such as a dilution series. For runs that
#'   will be charted, name the good baseline runs: with `NULL`, a later
#'   run that misses a low spike raises the range, and so changes the fit
#'   figures, of every run.
#'
#' @return A `data.frame` with one row per sample and columns `run`,
#'   `nSpikes` (spikes with an expected amount above 0), `nDetected` (of
#'   those, the ones detected), `slope`, `intercept`, `r2`, `residSD`,
#'   `offset`, `lod`, `lodStatus` and `allDetectedFrom`, as described
#'   above. Its `"minAmount"` attribute holds the amount used, to pass on
#'   to [plotDoseResponse()].
#'
#' @examples
#' set.seed(1)
#' runs <- simulateSpikeRuns(erccDesign(), nRuns = 3)
#' runMetrics(runs)
#' @export
runMetrics <- function(se, assay = 1L, expected = "expected", detect = 0,
                       lodProb = 0.95, minAmount = "auto",
                       baseline = NULL) {
    obs <- as.matrix(SummarizedExperiment::assay(se, assay))
    expd <- as.matrix(SummarizedExperiment::assay(se, expected))
    if (!identical(dim(obs), dim(expd)))
        stop("the measured and expected assays differ in shape")
    if (identical(minAmount, "auto")) {
        use <- .runIndex(baseline, colnames(se), ncol(obs))
        minAmount <- .autoMinAmount(obs[, use, drop = FALSE],
                                    expd[, use, drop = FALSE], detect)
        if (is.na(minAmount))
            warning("no expected amount is detected in every baseline run; ",
                    "the fit figures are NA. Set minAmount or name good ",
                    "runs in 'baseline'")
    } else if (!is.numeric(minAmount) || length(minAmount) != 1L ||
               is.na(minAmount) || minAmount < 0) {
        stop("'minAmount' must be \"auto\" or one amount of at least 0")
    }
    fitFrom <- if (is.na(minAmount)) Inf else minAmount
    rows <- lapply(seq_len(ncol(se)), function(j)
        .oneRun(obs[, j], expd[, j], detect, lodProb, fitFrom))
    out <- do.call(rbind, rows)
    res <- data.frame(run = colnames(se) %||% as.character(seq_len(ncol(se))),
                      out, row.names = NULL, stringsAsFactors = FALSE)
    attr(res, "minAmount") <- minAmount
    res
}

# Runs named by a baseline argument, as indices.
.runIndex <- function(baseline, names, n) {
    if (is.null(baseline))
        return(seq_len(n))
    if (is.logical(baseline)) {
        if (length(baseline) != n || anyNA(baseline))
            stop("a logical baseline needs one TRUE or FALSE per run")
        return(which(baseline))
    }
    if (is.character(baseline)) {
        idx <- match(baseline, names)
        if (anyNA(idx))
            stop("baseline names a run that is not in the data")
        return(idx)
    }
    idx <- if (length(baseline) == 1L) seq_len(baseline) else baseline
    if (any(idx < 1L | idx > n))
        stop("baseline indices fall outside the runs")
    idx
}

# Smallest expected amount at and above which every spike is detected in
# every run given; NA when even the largest amount is missed somewhere.
.autoMinAmount <- function(o, e, detect) {
    ok <- !is.na(e) & e > 0
    det <- !is.na(o) & o > detect
    for (a in sort(unique(e[ok]))) {
        sel <- ok & e >= a
        if (all(det[sel]))
            return(a)
    }
    NA_real_
}

.oneRun <- function(o, e, detect, lodProb, minAmount) {
    keep <- !is.na(e) & e > 0
    o <- o[keep]
    e <- e[keep]
    det <- !is.na(o) & o > detect
    res <- data.frame(nSpikes = length(e), nDetected = sum(det),
                      slope = NA_real_, intercept = NA_real_, r2 = NA_real_,
                      residSD = NA_real_, offset = NA_real_, lod = NA_real_,
                      lodStatus = NA_character_, allDetectedFrom = NA_real_)
    if (!length(e))
        return(res)
    use <- det & e >= minAmount
    if (any(use)) {
        lo <- log2(o[use])
        le <- log2(e[use])
        res$offset <- stats::median(lo - le)
        if (sum(use) >= 3L && length(unique(le)) > 1L) {
            fit <- stats::lm(lo ~ le)
            cf <- stats::coef(fit)
            r <- stats::residuals(fit)
            res$intercept <- unname(cf[1])
            res$slope <- unname(cf[2])
            ss <- sum((lo - mean(lo))^2)
            # A run whose fitted spikes all read the same has no variance to explain.
            res$r2 <- if (ss > 0) 1 - sum(r^2) / ss else NA_real_
            res$residSD <- sqrt(sum(r^2) / (length(r) - 2L))
        }
    }
    lodFit <- .lod(e, det, lodProb)
    res$lod <- lodFit$lod
    res$lodStatus <- lodFit$status
    res$allDetectedFrom <- .allDetectedFrom(e, det)
    res
}

.lod <- function(e, det, p) {
    none <- function(status) list(lod = NA_real_, status = status)
    if (!any(det))
        return(none("none_detected"))
    if (all(det))
        return(none("all_detected"))
    if (max(e[!det]) < min(e[det]))
        return(none("separated"))
    x <- log2(e)
    fit <- tryCatch(suppressWarnings(
        stats::glm(det ~ x, family = stats::binomial())),
        error = function(err) NULL)
    if (is.null(fit) || !fit$converged)
        return(none("fit_failed"))
    cf <- stats::coef(fit)
    if (anyNA(cf) || cf[2] <= 0)
        return(none("fit_failed"))
    lx <- (stats::qlogis(p) - cf[1]) / cf[2]
    if (lx < min(x) || lx > max(x))
        return(none("outside_range"))
    list(lod = unname(2^lx), status = "fit")
}

# Smallest amount at and above which every spike was detected.
.allDetectedFrom <- function(e, det) {
    if (!any(det))
        return(NA_real_)
    if (all(det))
        return(min(e))
    above <- e[e > max(e[!det])]
    if (length(above)) min(above) else NA_real_
}

#' Fold-change accuracy between two groups of runs
#'
#' Compares the measured log2 ratio of each spike between two groups of
#' samples with its expected log2 ratio. Ratios use the mean of each group.
#' A spike enters only if it is detected in every sample of both groups.
#'
#' @param se A `SummarizedExperiment` with an expected-amount assay.
#' @param groupA,groupB Column indices, names or logical vectors selecting
#'   the two groups. The ratio is A over B.
#' @param assay Name or index of the measured assay.
#' @param expected Name of the expected-amount assay.
#' @param center How to remove the overall scale between the groups.
#'   `"spikes"` subtracts the median error over all spikes; `"none"` keeps
#'   the raw ratios; a logical vector over the rows of `se` subtracts the
#'   median measured log2 ratio of those rows, for example a background
#'   proteome whose true ratio is 1.
#' @param detect Detection threshold, as in [runMetrics()].
#'
#' @return A list of two `data.frame`s. `features` has one row per spike
#'   with columns `feature`, `expectedLog2Ratio`, `measuredLog2Ratio`
#'   (after centring) and `error` (measured minus expected). `summary` has
#'   one row per expected ratio with columns `expectedLog2Ratio`,
#'   `nSpikes`, `medianError` and `rmse` (root mean square error).
#'
#' @examples
#' ercc <- erccDesign()
#' set.seed(1)
#' a <- simulateSpikeRuns(ercc, nRuns = 3, mix = "Mix1")
#' set.seed(2)
#' b <- simulateSpikeRuns(ercc, nRuns = 3, mix = "Mix2")
#' colnames(b) <- paste0("b", seq_len(3))
#' both <- cbind(a, b)
#' ratioAccuracy(both, 1:3, 4:6)$summary
#' @export
ratioAccuracy <- function(se, groupA, groupB, assay = 1L,
                          expected = "expected", center = "spikes",
                          detect = 0) {
    obs <- as.matrix(SummarizedExperiment::assay(se, assay))
    expd <- as.matrix(SummarizedExperiment::assay(se, expected))
    a <- .cols(se, groupA)
    b <- .cols(se, groupB)
    if (!length(a) || !length(b))
        stop("both groups need at least one sample")
    if (length(intersect(a, b)))
        stop("'groupA' and 'groupB' share a sample")
    ea <- rowMeans(expd[, a, drop = FALSE])
    eb <- rowMeans(expd[, b, drop = FALSE])
    oa <- obs[, a, drop = FALSE]
    ob <- obs[, b, drop = FALSE]
    seen <- function(m) rowSums(is.na(m) | m <= detect) == 0
    measuredLR <- log2(rowMeans(oa)) - log2(rowMeans(ob))
    spike <- !is.na(ea) & !is.na(eb) & ea > 0 & eb > 0 & seen(oa) & seen(ob)
    if (!any(spike))
        stop("no spike is detected in every sample of both groups")
    expectedLR <- log2(ea) - log2(eb)
    shift <- 0
    if (is.logical(center)) {
        if (length(center) != nrow(se))
            stop("a logical 'center' needs one value per row")
        ref <- center & seen(oa) & seen(ob)
        if (!any(ref))
            stop("no reference row is detected in every sample")
        shift <- stats::median(measuredLR[ref])
    } else if (identical(center, "spikes")) {
        shift <- stats::median(measuredLR[spike] - expectedLR[spike])
    } else if (!identical(center, "none")) {
        stop("'center' must be \"spikes\", \"none\" or a logical vector")
    }
    f <- data.frame(feature = rownames(se)[spike],
                    expectedLog2Ratio = expectedLR[spike],
                    measuredLog2Ratio = measuredLR[spike] - shift,
                    stringsAsFactors = FALSE)
    f$error <- f$measuredLog2Ratio - f$expectedLog2Ratio
    key <- round(f$expectedLog2Ratio, 3)
    s <- do.call(rbind, lapply(split(f$error, key), function(err)
        data.frame(nSpikes = length(err), medianError = stats::median(err),
                   rmse = sqrt(mean(err^2)))))
    s <- data.frame(expectedLog2Ratio = as.numeric(rownames(s)), s,
                    row.names = NULL)
    list(features = f, summary = s[order(-s$expectedLog2Ratio), ])
}

.cols <- function(se, g) {
    idx <- seq_len(ncol(se))
    names(idx) <- colnames(se)
    unname(idx[g])
}
