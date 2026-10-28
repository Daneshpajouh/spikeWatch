#' Chart limits for a target false-alarm probability
#'
#' Sets the CUSUM decision interval `h` or the EWMA limit width `L` so that
#' a chart from [driftChart()] signals by mistake within `horizon`
#' in-control runs with probability about `alpha`, counting the error in
#' the centre and sigma estimated from a short baseline. The fixed
#' defaults (`h = 5`, `L = 3`) ignore that error, and with 20 baseline
#' runs they can give many more false alarms than their textbook rates.
#'
#' Bootstrap: each of `B` replicates resamples the baseline values with
#' replacement, re-estimates the centre and sigma from the resample as
#' [driftChart()] does, draws `horizon` future in-control values from the
#' baseline values, and records the largest chart statistic over the
#' horizon under the re-estimated centre and sigma. The limit is the
#' `1 - alpha` quantile of those maxima.
#'
#' This calibrates the false-alarm probability averaged over baselines.
#' It is not the conditional guarantee of Gandy and Kvaløy (2013), which
#' bounds the false-alarm rate for the baseline actually drawn with a
#' stated probability. It treats baseline runs as independent and draws
#' future values only from the baseline values seen, so a very short
#' baseline understates the tails.
#'
#' @inheritParams driftChart
#' @param method `"cusum"` or `"ewma"`.
#' @param alpha Target probability of any false alarm within `horizon`
#'   runs.
#' @param horizon Number of monitored runs the target applies to.
#' @param B Bootstrap replicates. Call [set.seed()] first for a
#'   repeatable limit.
#'
#' @return A named number, `h` for CUSUM or `L` for EWMA, to pass to
#'   [driftChart()]. Attributes hold `alpha`, `horizon`, `B` and `method`.
#'
#' @references Gandy A, Kvaløy JT (2013). Guaranteed conditional
#'   performance of control charts via bootstrap methods. Scandinavian
#'   Journal of Statistics 40(4), 647-668. \doi{10.1002/sjos.12006}
#'
#' @examples
#' set.seed(1)
#' x <- c(rnorm(20, 1, 0.02), rnorm(10, 0.97, 0.02))
#' h <- calibrateChart(x, baseline = 20, method = "cusum", horizon = 10,
#'                     B = 500)
#' ch <- driftChart(x, baseline = 20, method = "cusum", h = h)
#' firstSignal(ch)
#' @export
calibrateChart <- function(x, metric = NULL, baseline = 20L,
                           method = c("cusum", "ewma"), alpha = 0.05,
                           horizon = 20L, B = 1000L, lambda = 0.2, k = 0.5,
                           sigma = c("mr", "sd")) {
    method <- match.arg(method)
    sigma <- match.arg(sigma)
    if (alpha <= 0 || alpha >= 1)
        stop("'alpha' must be in (0, 1)")
    if (horizon < 1L || B < 1L / alpha)
        stop("'horizon' must be at least 1 and 'B' at least 1 / alpha")
    if (is.data.frame(x)) {
        if (is.null(metric) || !metric %in% names(x))
            stop("name the column to chart with 'metric'")
        x <- x[[metric]]
    }
    x <- as.numeric(x)
    base <- if (length(baseline) == 1L) seq_len(baseline) else baseline
    if (any(base < 1L | base > length(x)))
        stop("baseline indices fall outside the series")
    xb <- x[sort(unique(as.integer(base)))]
    xb <- xb[!is.na(xb)]
    n <- length(xb)
    if (n < 3L)
        stop("the baseline needs at least three non-missing runs")
    xs <- matrix(xb[sample.int(n, n * B, replace = TRUE)], B)
    mu <- rowMeans(xs)
    s <- if (sigma == "mr")
        rowMeans(abs(xs[, -1L, drop = FALSE] - xs[, -n, drop = FALSE])) /
            1.128 else apply(xs, 1L, stats::sd)
    fut <- matrix(xb[sample.int(n, B * horizon, replace = TRUE)], B)
    u <- (fut - mu) / s
    u[s <= 0, ] <- Inf
    m <- if (method == "cusum") .cusumMax(u, k) else .ewmaMax(u, lambda)
    thr <- sort(m)[ceiling((1 - alpha) * B)]
    out <- stats::setNames(thr, if (method == "cusum") "h" else "L")
    attr(out, "alpha") <- alpha
    attr(out, "horizon") <- horizon
    attr(out, "B") <- B
    attr(out, "method") <- method
    out
}

# Largest two-sided CUSUM sum over the columns of u, one value per row.
.cusumMax <- function(u, k) {
    cu <- cl <- m <- numeric(nrow(u))
    for (t in seq_len(ncol(u))) {
        cu <- pmax(0, cu + u[, t] - k)
        cl <- pmax(0, cl - u[, t] - k)
        m <- pmax(m, cu, cl)
    }
    m
}

# Largest EWMA distance from the centre in units of its exact standard
# deviation at each step, one value per row.
.ewmaMax <- function(u, lambda) {
    z <- m <- numeric(nrow(u))
    for (t in seq_len(ncol(u))) {
        z <- lambda * u[, t] + (1 - lambda) * z
        w <- sqrt(lambda / (2 - lambda) * (1 - (1 - lambda)^(2 * t)))
        m <- pmax(m, abs(z) / w)
    }
    m
}
