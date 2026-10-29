#' Plot measured against expected spike-in amounts for one run
#'
#' @param se A `SummarizedExperiment` with an expected-amount assay.
#' @param run Column name or index of the run to plot.
#' @param assay Name or index of the measured assay.
#' @param expected Name of the expected-amount assay.
#' @param detect Detection threshold, as in [runMetrics()].
#' @param minAmount The fit range, as in [runMetrics()]: the line is
#'   fitted to detected spikes at or above it, and a dashed line marks it.
#'   `"auto"` (the default) chooses it from the runs named by `baseline`
#'   exactly as [runMetrics()] does, so with the same `detect` and
#'   `baseline` the line is the fit behind that run's `slope`. A number
#'   fixes it (0 fits every detected spike); `NA` draws no line.
#' @param baseline Runs used to choose `minAmount = "auto"`, as in
#'   [runMetrics()]. `NULL` (the default) uses every run.
#'
#' @return A `ggplot` object. Undetected spikes are drawn at the bottom
#'   edge as open points. A run with no spike above `detect` is an error.
#'
#' @examples
#' set.seed(1)
#' runs <- simulateSpikeRuns(erccDesign(), nRuns = 2)
#' plotDoseResponse(runs, 1)
#' @export
plotDoseResponse <- function(se, run = 1L, assay = 1L,
                             expected = "expected", detect = 0,
                             minAmount = "auto", baseline = NULL) {
    noAuto <- FALSE
    if (identical(minAmount, "auto")) {
        obs <- as.matrix(SummarizedExperiment::assay(se, assay))
        expd <- as.matrix(SummarizedExperiment::assay(se, expected))
        use <- .runIndex(baseline, colnames(se), ncol(obs))
        minAmount <- .autoMinAmount(obs[, use, drop = FALSE],
                                    expd[, use, drop = FALSE], detect)
        noAuto <- is.na(minAmount)
    } else if (length(minAmount) != 1L || !(is.na(minAmount) ||
               (is.numeric(minAmount) && minAmount >= 0))) {
        stop("'minAmount' must be \"auto\", NA or one amount of at least 0")
    }
    o <- SummarizedExperiment::assay(se, assay)[, run]
    e <- SummarizedExperiment::assay(se, expected)[, run]
    keep <- !is.na(e) & e > 0
    df <- data.frame(expected = e[keep], measured = o[keep],
                     group = SummarizedExperiment::rowData(se)$spikeGroup[keep])
    df$detected <- !is.na(df$measured) & df$measured > detect
    if (!any(df$detected))
        stop("no spike is above detect = ", detect, " in run ", run,
             ", so there is nothing to plot on the log axes")
    if (noAuto)
        warning("no expected amount is detected in every baseline run; ",
                "no fit line is drawn. Set minAmount or name good runs in ",
                "'baseline'")
    floor <- min(df$measured[df$detected], 1, na.rm = TRUE)
    df$measured[!df$detected] <- floor
    df$group[is.na(df$group)] <- "spike"
    label <- if (is.character(run)) run else colnames(se)[run]
    if (is.null(label) || is.na(label)) label <- as.character(run)
    fitted <- if (is.na(minAmount)) df[0, ] else
        df[df$detected & df$expected >= minAmount, ]
    # The rule runMetrics() fits by: three or more spikes over two or more
    # amounts, or no line (code review round 8).
    if (nrow(fitted) < 3L || length(unique(fitted$expected)) < 2L)
        fitted <- df[0, ]
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$expected, .data$measured)) +
        ggplot2::geom_point(ggplot2::aes(shape = .data$detected,
                                         colour = .data$group)) +
        ggplot2::geom_smooth(data = fitted, method = "lm",
                             formula = y ~ x, se = FALSE, colour = "grey30",
                             linewidth = 0.4) +
        ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
        ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
        ggplot2::labs(x = "Expected amount", y = "Measured",
                      title = label,
                      shape = "Detected", colour = "Group") +
        ggplot2::theme_bw()
    if (!is.na(minAmount) && minAmount > 0)
        p <- p + ggplot2::geom_vline(xintercept = minAmount, linetype = 2,
                                     colour = "grey50")
    p
}

#' Plot a control chart
#'
#' @param chart A chart from [driftChart()].
#' @param title Plot title.
#'
#' @return A `ggplot` object: the chart statistic by run, its limits, the
#'   baseline runs shaded and signalling runs in red. EWMA and CUSUM charts
#'   start after the baseline, so their baseline runs show no statistic.
#'
#' @examples
#' ch <- driftChart(c(rep(c(1, 1.01, 0.99, 1.02), 5), 0.9, 0.85),
#'                  baseline = 20, method = "ewma")
#' plotDrift(ch)
#' @importFrom rlang .data
#' @export
plotDrift <- function(chart, title = NULL) {
    df <- chart
    df$i <- seq_len(nrow(df))
    b <- which(df$baseline)
    ggplot2::ggplot(df, ggplot2::aes(.data$i, .data$stat)) +
        ggplot2::annotate("rect", xmin = b - 0.5, xmax = b + 0.5,
                          ymin = -Inf, ymax = Inf, alpha = 0.1) +
        ggplot2::geom_line(colour = "grey40", na.rm = TRUE) +
        ggplot2::geom_point(ggplot2::aes(colour = .data$signal),
                            na.rm = TRUE) +
        ggplot2::geom_line(ggplot2::aes(y = .data$upper), linetype = 2,
                           na.rm = TRUE) +
        ggplot2::geom_line(ggplot2::aes(y = .data$lower), linetype = 2,
                           na.rm = TRUE) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = .data$center),
                            colour = "grey60") +
        ggplot2::scale_colour_manual(values = c(`FALSE` = "black",
                                                `TRUE` = "red")) +
        ggplot2::labs(x = "Run", y = attr(chart, "method"),
                      colour = "Signal", title = title) +
        ggplot2::theme_bw()
}
