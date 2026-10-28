#' Simulate a series of runs with known spike-in truth
#'
#' Makes count data for a series of runs of one spike-in mix, with an
#' optional drift planted at a known run. It exists to test charts and
#' thresholds before trusting them on real data: the truth of every run,
#' and whether it drifted, is in the output.
#'
#' For run `j` and spike `i`, the log2 mean count is
#' `log2(scale) + depth_j + slope_j * log2(amount_i / pivot) +
#' log2(pivot) + noise`, where `depth_j` varies from run to run as normal
#' with standard deviation `depthSD` (log2 units), `pivot` is the median
#' design amount and `noise` is normal with standard deviation `noise`.
#' Counts are Poisson around that mean, so low spikes drop out the way
#' they do in real runs. Depth changes from run to run are not drift.
#'
#' Drift types, applied from run `driftAt` on:
#' * `"slope"`: `slope_j` changes by `driftSize` (compression when
#'   negative). The slope turns about the median amount, so middle spikes
#'   keep their level while high spikes fall and low spikes rise. Total
#'   spike counts still fall, because the highest spikes dominate them.
#' * `"offset"`: `depth_j` changes by `driftSize`, a loss or gain of
#'   sensitivity across all spikes.
#' * `"noise"`: the noise standard deviation grows by `driftSize`.
#'
#' With `driftShape = "ramp"` the change grows linearly to `driftSize`
#' over `rampRuns` runs instead of stepping.
#'
#' @param design A design from [spikeDesign()] or [erccDesign()].
#' @param nRuns Number of runs.
#' @param mix The design mix used in every run. Defaults to the first.
#' @param scale Mean count per unit amount at depth 0.
#' @param slope Dose-response slope before any drift.
#' @param noise Standard deviation of the log2 noise before any drift.
#' @param depthSD Standard deviation of the run depth, in log2 units.
#' @param driftAt Run at which drift starts, or `NA` for none.
#' @param driftType One of `"slope"`, `"offset"` or `"noise"`.
#' @param driftSize Size of the change.
#' @param driftShape `"step"` or `"ramp"`.
#' @param rampRuns Runs over which a ramp reaches full size.
#' @param seed Optional random seed.
#'
#' @return A `SummarizedExperiment` with assays `counts` and `expected`,
#'   and `colData` columns `mix`, `depth` (log2) and `drifted`.
#'
#' @examples
#' runs <- simulateSpikeRuns(erccDesign(), nRuns = 30, driftAt = 21,
#'                           driftType = "slope", driftSize = -0.1,
#'                           seed = 42)
#' m <- runMetrics(runs)
#' firstSignal(driftChart(m, "slope", baseline = 15))
#' @export
simulateSpikeRuns <- function(design, nRuns = 30L, mix = NULL, scale = 10,
                              slope = 1, noise = 0.25, depthSD = 0.5,
                              driftAt = NA, driftType = c("slope", "offset",
                                                          "noise"),
                              driftSize = 0, driftShape = c("step", "ramp"),
                              rampRuns = 5L, seed = NULL) {
    driftType <- match.arg(driftType)
    driftShape <- match.arg(driftShape)
    if (!is.null(seed))
        set.seed(seed)
    mix <- mix %||% design$mix[1]
    d <- design[design$mix == mix, , drop = FALSE]
    if (!nrow(d))
        stop("the design has no mix '", mix, "'")
    j <- seq_len(nRuns)
    step <- if (is.na(driftAt)) numeric(nRuns) else if (
        driftShape == "step") as.numeric(j >= driftAt) else
        pmin(pmax((j - driftAt + 1) / rampRuns, 0), 1)
    change <- driftSize * step
    depth <- stats::rnorm(nRuns, 0, depthSD)
    slopes <- rep(slope, nRuns)
    sds <- rep(noise, nRuns)
    if (driftType == "slope") slopes <- slopes + change
    if (driftType == "offset") depth <- depth + change
    if (driftType == "noise") sds <- sds + change
    pivot <- log2(stats::median(d$amount))
    la <- log2(d$amount)
    lmu <- outer(la - pivot, slopes) + pivot + log2(scale)
    lmu <- sweep(lmu, 2, depth, "+")
    lmu <- lmu + sweep(matrix(stats::rnorm(length(lmu)), nrow(lmu)), 2,
                       sds, "*")
    counts <- matrix(stats::rpois(length(lmu), 2^lmu), nrow(lmu),
                     dimnames = list(d$feature, sprintf("run%02d", j)))
    cd <- S4Vectors::DataFrame(mix = rep(mix, nRuns), depth = depth,
                               drifted = change != 0,
                               row.names = colnames(counts))
    se <- SummarizedExperiment::SummarizedExperiment(
        list(counts = counts), colData = cd)
    addTruth(se, d, mix = "mix")
}
