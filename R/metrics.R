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
#' * `residSD`: standard deviation of the fit residuals, a precision figure.
#' * `offset`: median of log2 measured minus log2 expected, over the same
#'   spikes as the fit. It moves with sequencing depth or instrument
#'   response, and is defined even when all spikes in a run share one
#'   amount.
#' * `lod`: the expected amount at which a spike is detected with
#'   probability `lodProb`, from a logistic fit of detection on log2
#'   expected amount. When every spike is detected, `lod` is the smallest
#'   expected amount and `lodCensored` is `TRUE`: the true limit is at or
#'   below it. When the fit fails or its slope is not positive, `lod` is
#'   the smallest amount above which every spike was detected. The same
#'   rule applies when detected and missed spikes do not overlap, where the
#'   logistic fit has no finite estimate.
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
#'   amount, not by what was measured, avoids that. Set it to an amount
#'   every good run detects; `lod` always uses all spikes.
#'
#' @return A `data.frame` with one row per sample.
#'
#' @examples
#' runs <- simulateSpikeRuns(erccDesign(), nRuns = 3, seed = 1)
#' runMetrics(runs)
#' @export
runMetrics <- function(se, assay = 1L, expected = "expected", detect = 0,
                       lodProb = 0.95, minAmount = 0) {
    obs <- as.matrix(SummarizedExperiment::assay(se, assay))
    expd <- as.matrix(SummarizedExperiment::assay(se, expected))
    if (!identical(dim(obs), dim(expd)))
        stop("the measured and expected assays differ in shape")
    rows <- lapply(seq_len(ncol(se)), function(j)
        .oneRun(obs[, j], expd[, j], detect, lodProb, minAmount))
    out <- do.call(rbind, rows)
    data.frame(run = colnames(se) %||% as.character(seq_len(ncol(se))),
               out, row.names = NULL, stringsAsFactors = FALSE)
}

.oneRun <- function(o, e, detect, lodProb, minAmount) {
    keep <- !is.na(e) & e > 0
    o <- o[keep]
    e <- e[keep]
    det <- !is.na(o) & o > detect
    res <- data.frame(nSpikes = length(e), nDetected = sum(det),
                      slope = NA_real_, intercept = NA_real_, r2 = NA_real_,
                      residSD = NA_real_, offset = NA_real_, lod = NA_real_,
                      lodCensored = NA)
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
            res$r2 <- 1 - sum(r^2) / sum((lo - mean(lo))^2)
            res$residSD <- stats::sd(r)
        }
    }
    lodFit <- .lod(e, det, lodProb)
    res$lod <- lodFit$lod
    res$lodCensored <- lodFit$censored
    res
}

.lod <- function(e, det, p) {
    if (!any(det))
        return(list(lod = NA_real_, censored = NA))
    if (all(det))
        return(list(lod = min(e), censored = TRUE))
    # Complete separation: the logistic fit has no finite estimate, so use
    # the smallest amount above which every spike was detected.
    missed <- max(e[!det])
    if (missed < min(e[det]))
        return(list(lod = min(e[det]), censored = FALSE))
    x <- log2(e)
    fit <- tryCatch(suppressWarnings(
        stats::glm(det ~ x, family = stats::binomial())),
        error = function(err) NULL)
    if (!is.null(fit) && fit$converged) {
        cf <- stats::coef(fit)
        if (!anyNA(cf) && cf[2] > 0) {
            lx <- (stats::qlogis(p) - cf[1]) / cf[2]
            lx <- min(max(lx, min(x)), max(x))
            return(list(lod = unname(2^lx), censored = FALSE))
        }
    }
    above <- e[e > missed]
    list(lod = if (length(above)) min(above) else NA_real_, censored = FALSE)
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
#' @return A list with `features`, one row per spike (expected and
#'   measured log2 ratio, error), and `summary`, one row per expected
#'   ratio (number of spikes, median error, root mean square error).
#'
#' @examples
#' ercc <- erccDesign()
#' a <- simulateSpikeRuns(ercc, nRuns = 3, mix = "Mix1", seed = 1)
#' b <- simulateSpikeRuns(ercc, nRuns = 3, mix = "Mix2", seed = 2)
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
    ea <- rowMeans(expd[, a, drop = FALSE])
    eb <- rowMeans(expd[, b, drop = FALSE])
    oa <- obs[, a, drop = FALSE]
    ob <- obs[, b, drop = FALSE]
    seen <- function(m) rowSums(is.na(m) | m <= detect) == 0
    measuredLR <- log2(rowMeans(oa)) - log2(rowMeans(ob))
    spike <- !is.na(ea) & !is.na(eb) & ea > 0 & eb > 0 & seen(oa) & seen(ob)
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
