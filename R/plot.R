#' Plot measured against expected spike-in amounts for one run
#'
#' @param se A `SummarizedExperiment` with an expected-amount assay.
#' @param run Column name or index of the run to plot.
#' @param assay Name or index of the measured assay.
#' @param expected Name of the expected-amount assay.
#' @param detect Detection threshold, as in [runMetrics()].
#'
#' @return A `ggplot` object. Undetected spikes are drawn at the bottom
#'   edge as open points.
#'
#' @examples
#' runs <- simulateSpikeRuns(erccDesign(), nRuns = 2, seed = 1)
#' plotDoseResponse(runs, 1)
#' @export
plotDoseResponse <- function(se, run = 1L, assay = 1L,
                             expected = "expected", detect = 0) {
    o <- SummarizedExperiment::assay(se, assay)[, run]
    e <- SummarizedExperiment::assay(se, expected)[, run]
    keep <- !is.na(e) & e > 0
    df <- data.frame(expected = e[keep], measured = o[keep],
                     group = SummarizedExperiment::rowData(se)$spikeGroup[keep])
    df$detected <- !is.na(df$measured) & df$measured > detect
    floor <- min(df$measured[df$detected], 1, na.rm = TRUE)
    df$measured[!df$detected] <- floor
    df$group[is.na(df$group)] <- "spike"
    ggplot2::ggplot(df, ggplot2::aes(.data$expected, .data$measured)) +
        ggplot2::geom_point(ggplot2::aes(shape = .data$detected,
                                         colour = .data$group)) +
        ggplot2::geom_smooth(data = df[df$detected, ], method = "lm",
                             formula = y ~ x, se = FALSE, colour = "grey30",
                             linewidth = 0.4) +
        ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
        ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
        ggplot2::labs(x = "Expected amount", y = "Measured",
                      title = colnames(se)[run] %||% as.character(run),
                      shape = "Detected", colour = "Group") +
        ggplot2::theme_bw()
}

#' Plot a control chart
#'
#' @param chart A chart from [driftChart()].
#' @param title Plot title.
#'
#' @return A `ggplot` object: the chart statistic by run, its limits, the
#'   baseline runs shaded and signalling runs in red.
#'
#' @examples
#' ch <- driftChart(c(1, 1.01, 0.99, 1, 1.02, 0.9, 0.85), baseline = 5,
#'                  method = "ewma")
#' plotDrift(ch)
#' @importFrom rlang .data
#' @export
plotDrift <- function(chart, title = NULL) {
    df <- chart
    df$i <- seq_len(nrow(df))
    nb <- max(which(df$baseline))
    ggplot2::ggplot(df, ggplot2::aes(.data$i, .data$stat)) +
        ggplot2::annotate("rect", xmin = 0.5, xmax = nb + 0.5, ymin = -Inf,
                          ymax = Inf, alpha = 0.1) +
        ggplot2::geom_line(colour = "grey40") +
        ggplot2::geom_point(ggplot2::aes(colour = .data$signal)) +
        ggplot2::geom_line(ggplot2::aes(y = .data$upper), linetype = 2) +
        ggplot2::geom_line(ggplot2::aes(y = .data$lower), linetype = 2) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = .data$center),
                            colour = "grey60") +
        ggplot2::scale_colour_manual(values = c(`FALSE` = "black",
                                                `TRUE` = "red")) +
        ggplot2::labs(x = "Run", y = attr(chart, "method"),
                      colour = "Signal", title = title) +
        ggplot2::theme_bw()
}
