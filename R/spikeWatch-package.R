#' spikeWatch: run accuracy and drift from spike-in truth
#'
#' Start with [spikeDesign()] or [erccDesign()], attach the truth to an
#' experiment with [addTruth()], compute per-run figures with
#' [runMetrics()] and chart them across runs with [driftChart()].
#' [ratioAccuracy()] checks fold changes between two groups, and
#' [simulateSpikeRuns()] makes series with planted drift for testing.
#'
#' @keywords internal
"_PACKAGE"
