#' Simulate a series of runs with known spike-in truth
#'
#' Makes count data for a series of runs of one spike-in mix, with an
#' optional drift planted at a known run. Call [set.seed()] first for a
#' repeatable series. It exists to test charts and
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
#' Departures from that model, for testing a chart away from the model
#' the fit assumes. Each is off by default, and the defaults give the same
#' counts as earlier versions for the same seed.
#' * `pivot = "top"`: the slope turns about the highest design amount, so
#'   the top spike keeps its level and compression raises every lower
#'   spike towards it, as a rising background does. Total spike counts
#'   then barely move.
#' * `counts = "nbinom"`: negative binomial counts with size `size` in
#'   place of Poisson.
#' * `saturation`: above a fixed ceiling, the log2 mean is compressed by
#'   this fraction (0.5 halves the rise). The ceiling is `log2(scale)`
#'   plus the `saturationAt` quantile of the log2 design amounts: the
#'   log2 mean of that amount at depth 0 when the slope is 1. It does not
#'   move with the slope, as a detector's ceiling would not.
#' * `spikeBiasSD`: each spike gets a fixed log2 bias, normal with this
#'   standard deviation, drawn once and kept in every run.
#' * `slopeSD`: the slope varies from run to run as normal with this
#'   standard deviation, before any drift.
#' * `runCor`: run depth and the slope wobble follow AR(1) with this
#'   correlation between neighbouring runs, at the standard deviations
#'   given.
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
#' @param pivot `"median"` or `"top"`: the amount the slope turns about.
#' @param counts `"poisson"` or `"nbinom"`.
#' @param size Negative binomial size; smaller is more overdispersed.
#' @param saturation Fraction by which the log2 mean is compressed above
#'   the saturation level; 0 for none.
#' @param saturationAt Quantile of the log2 design amounts. The ceiling
#'   is `log2(scale)` plus that quantile.
#' @param spikeBiasSD Standard deviation of the fixed per-spike log2 bias.
#' @param slopeSD Standard deviation of the run-to-run slope.
#' @param runCor Correlation between neighbouring runs of depth and slope.
#'
#' @return A `SummarizedExperiment` with assays `counts` and `expected`,
#'   and `colData` columns `mix`, `depth` (log2) and `drifted`.
#'
#' @examples
#' set.seed(42)
#' runs <- simulateSpikeRuns(erccDesign(), nRuns = 30, driftAt = 21,
#'                           driftType = "slope", driftSize = -0.1)
#' m <- runMetrics(runs)
#' firstSignal(driftChart(m, "slope", baseline = 15))
#' @export
simulateSpikeRuns <- function(design, nRuns = 30L, mix = NULL, scale = 10,
                              slope = 1, noise = 0.25, depthSD = 0.5,
                              driftAt = NA, driftType = c("slope", "offset",
                                                          "noise"),
                              driftSize = 0, driftShape = c("step", "ramp"),
                              rampRuns = 5L, pivot = c("median", "top"),
                              counts = c("poisson", "nbinom"), size = 10,
                              saturation = 0, saturationAt = 0.9,
                              spikeBiasSD = 0, slopeSD = 0, runCor = 0) {
    driftType <- match.arg(driftType)
    driftShape <- match.arg(driftShape)
    pivot <- match.arg(pivot)
    counts <- match.arg(counts)
    if (runCor < 0 || runCor >= 1)
        stop("'runCor' must be in [0, 1)")
    if (saturation < 0 || saturation >= 1)
        stop("'saturation' must be in [0, 1)")
    mix <- mix %||% design$mix[1]
    d <- design[design$mix == mix, , drop = FALSE]
    if (!nrow(d))
        stop("the design has no mix '", mix, "'")
    j <- seq_len(nRuns)
    step <- if (is.na(driftAt)) numeric(nRuns) else if (
        driftShape == "step") as.numeric(j >= driftAt) else
        pmin(pmax((j - driftAt + 1) / rampRuns, 0), 1)
    change <- driftSize * step
    depth <- .ar1(stats::rnorm(nRuns), runCor) * depthSD
    slopes <- rep(slope, nRuns)
    if (slopeSD > 0)
        slopes <- slopes + .ar1(stats::rnorm(nRuns), runCor) * slopeSD
    sds <- rep(noise, nRuns)
    if (driftType == "slope") slopes <- slopes + change
    if (driftType == "offset") depth <- depth + change
    if (driftType == "noise") sds <- sds + change
    la <- log2(d$amount)
    piv <- if (pivot == "top") max(la) else log2(stats::median(d$amount))
    lmu <- outer(la - piv, slopes) + piv + log2(scale)
    if (spikeBiasSD > 0)
        lmu <- lmu + stats::rnorm(nrow(d), 0, spikeBiasSD)
    lmu <- sweep(lmu, 2, depth, "+")
    if (saturation > 0) {
        top <- log2(scale) + stats::quantile(la, saturationAt, names = FALSE)
        over <- lmu > top
        lmu[over] <- top + (1 - saturation) * (lmu[over] - top)
    }
    lmu <- lmu + sweep(matrix(stats::rnorm(length(lmu)), nrow(lmu)), 2,
                       sds, "*")
    y <- if (counts == "nbinom")
        stats::rnbinom(length(lmu), size = size, mu = 2^lmu) else
        stats::rpois(length(lmu), 2^lmu)
    y <- matrix(y, nrow(lmu),
                dimnames = list(d$feature, sprintf("run%02d", j)))
    cd <- S4Vectors::DataFrame(mix = rep(mix, nRuns), depth = depth,
                               drifted = change != 0,
                               row.names = colnames(y))
    se <- SummarizedExperiment::SummarizedExperiment(
        list(counts = y), colData = cd)
    addTruth(se, d, mix = "mix")
}

# Stationary AR(1) with unit variance from independent standard normals;
# with rho = 0 it returns e unchanged.
.ar1 <- function(e, rho) {
    if (rho == 0 || length(e) < 2L)
        return(e)
    out <- e
    for (i in 2:length(e))
        out[i] <- rho * out[i - 1L] + sqrt(1 - rho^2) * e[i]
    out
}
