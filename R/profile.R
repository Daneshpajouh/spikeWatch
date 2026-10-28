#' Truth-free figures from the spike profile of good runs
#'
#' What a lab can chart without the known amounts: each run's spike counts
#' are compared with the average profile of the same spikes over a
#' baseline of good runs. It is the strongest chart we know how to build
#' without truth, and the fair baseline for the truth-based figures of
#' [runMetrics()].
#'
#' Values are taken as log2 (value + `pseudo`). The spikes used are those
#' detected in every baseline run, and the reference profile is their mean
#' over the baseline runs; both are fixed from the baseline only.
#'
#' * `profileSlope`, `profileIntercept`: least-squares fit of each run's
#'   values on the reference profile. Compression gives a slope below 1.
#' * `profileResidSD`: residual standard error of that fit (divisor
#'   n - 2).
#' * `nSpikes`: spikes in the fit.
#'
#' The `"z"` attribute is a spike by run matrix: each value minus its
#' reference and minus the run's median shift, divided by that spike's
#' standard deviation of the same quantity over the baseline runs. It
#' shows which spikes moved in which runs.
#'
#' @param se A `SummarizedExperiment`. Rows with `rowData` column `spike`
#'   `TRUE` are used; without that column every row is.
#' @param assay Name or index of the measured assay.
#' @param baseline Good runs: a single number is the count of leading
#'   runs; otherwise run indices, run names or a logical vector.
#' @param detect A spike counts as detected when its value is above this
#'   threshold and not `NA`.
#' @param pseudo Added before the log2.
#'
#' @return A `data.frame` with one row per run and attributes `"z"`,
#'   `"reference"` (the profile, named by spike) and `"baseline"` (the run
#'   indices).
#'
#' @examples
#' set.seed(1)
#' runs <- simulateSpikeRuns(erccDesign(), nRuns = 25, driftAt = 21,
#'                           driftType = "slope", driftSize = -0.1)
#' p <- profileMetrics(runs, baseline = 20)
#' round(p$profileSlope, 3)
#' @export
profileMetrics <- function(se, assay = 1L, baseline = 20L, detect = 0,
                           pseudo = 0.5) {
    obs <- as.matrix(SummarizedExperiment::assay(se, assay))
    rd <- rowData(se)
    if ("spike" %in% names(rd))
        obs <- obs[!is.na(rd$spike) & rd$spike, , drop = FALSE]
    base <- .runIndex(baseline, colnames(se), ncol(obs))
    seen <- !is.na(obs[, base, drop = FALSE]) &
        obs[, base, drop = FALSE] > detect
    use <- rowSums(seen) == length(base)
    if (sum(use) < 3L)
        stop("fewer than three spikes are detected in every baseline run")
    y <- log2(obs[use, , drop = FALSE] + pseudo)
    ref <- rowMeans(y[, base, drop = FALSE])
    fits <- lapply(seq_len(ncol(y)), function(j) .profileFit(y[, j], ref))
    out <- data.frame(
        run = colnames(se) %||% as.character(seq_len(ncol(y))),
        profileSlope = vapply(fits, `[[`, 0, "slope"),
        profileIntercept = vapply(fits, `[[`, 0, "intercept"),
        profileResidSD = vapply(fits, `[[`, 0, "residSD"),
        nSpikes = vapply(fits, `[[`, 0, "n"),
        stringsAsFactors = FALSE)
    dev <- y - ref
    dev <- sweep(dev, 2, apply(dev, 2, stats::median, na.rm = TRUE))
    s <- apply(dev[, base, drop = FALSE], 1, stats::sd)
    s[!is.finite(s) | s <= 0] <- NA_real_
    attr(out, "z") <- dev / s
    attr(out, "reference") <- ref
    attr(out, "baseline") <- base
    out
}

.profileFit <- function(v, ref) {
    ok <- !is.na(v)
    n <- sum(ok)
    if (n < 3L)
        return(list(slope = NA_real_, intercept = NA_real_,
                    residSD = NA_real_, n = n))
    f <- stats::lm.fit(cbind(1, ref[ok]), v[ok])
    list(slope = unname(f$coefficients[2]),
         intercept = unname(f$coefficients[1]),
         residSD = sqrt(sum(f$residuals^2) / (n - 2)), n = n)
}
