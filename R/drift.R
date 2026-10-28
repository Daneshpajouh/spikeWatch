#' Control chart of one accuracy figure across runs
#'
#' Tracks a per-run figure from [runMetrics()], in run order, against its
#' behaviour over a baseline of runs known to be good. The centre is the
#' baseline mean. The spread `sigma` is estimated from the baseline by the
#' average moving range divided by 1.128 (`"mr"`, the usual estimate for
#' individual values) or by the standard deviation (`"sd"`). Moving ranges
#' are taken only between baseline runs that are next to each other in run
#' order and both observed, so a gap or a missing run never makes two
#' distant runs neighbours.
#'
#' * `"shewhart"`: an individuals chart. Signals when a run lies more than
#'   `L` sigma from the centre.
#' * `"ewma"`: exponentially weighted moving average with weight `lambda`
#'   and exact time-varying limits at `L` sigma.
#' * `"cusum"`: two-sided tabular CUSUM with reference value `k` sigma and
#'   decision interval `h` sigma. `stat` is the larger of the upper and
#'   lower sums.
#'
#' The EWMA and CUSUM statistics start fresh at the first run after the
#' last baseline run (the EWMA at the centre, the CUSUM sums at zero), so
#' monitoring never inherits what the baseline runs built up. Their rows up
#' to the last baseline run have `stat`, `lower` and `upper` set to `NA`
#' and never signal. The Shewhart chart has no memory and is shown for
#' every run.
#'
#' Runs with a missing value are carried: the chart statistic does not
#' move and the run cannot signal.
#'
#' @param x Numeric vector of the figure in run order, or a `data.frame`
#'   from [runMetrics()] together with `metric`.
#' @param metric Column of `x` to chart when `x` is a `data.frame`.
#' @param baseline Either the number of leading runs that form the
#'   baseline, or their indices. Fewer than 20 non-missing baseline runs
#'   give a warning: limits from a short baseline are themselves uncertain.
#' @param method One of `"ewma"`, `"cusum"` or `"shewhart"`.
#' @param lambda EWMA weight, between 0 and 1.
#' @param L Limit width in sigma for the Shewhart and EWMA charts.
#' @param k,h CUSUM reference value and decision interval, in sigma.
#' @param sigma Spread estimate: `"mr"` or `"sd"`.
#'
#' @return A `data.frame` with one row per run: `run`, `value`, `stat`,
#'   `center`, `lower`, `upper`, `baseline` and `signal`; a CUSUM chart
#'   also has `cusumUp` and `cusumLow`, the upper and lower sums in sigma
#'   units. Attributes hold the method, its parameters and the baseline
#'   `sigma`.
#'
#' @examples
#' set.seed(1)
#' x <- c(rnorm(20, 1, 0.02), rnorm(10, 0.95, 0.02))
#' ch <- driftChart(x, baseline = 20, method = "ewma")
#' firstSignal(ch)
#' @export
driftChart <- function(x, metric = NULL, baseline = 20L,
                       method = c("ewma", "cusum", "shewhart"),
                       lambda = 0.2, L = 3, k = 0.5, h = 5,
                       sigma = c("mr", "sd")) {
    method <- match.arg(method)
    sigma <- match.arg(sigma)
    # Limits may arrive named and with attributes.
    lambda <- as.vector(unname(lambda))
    L <- as.vector(unname(L))
    k <- as.vector(unname(k))
    h <- as.vector(unname(h))
    runs <- NULL
    if (is.data.frame(x)) {
        if (is.null(metric) || !metric %in% names(x))
            stop("name the column to chart with 'metric'")
        runs <- x$run
        x <- x[[metric]]
    }
    x <- as.numeric(x)
    n <- length(x)
    runs <- runs %||% as.character(seq_len(n))
    base <- if (length(baseline) == 1L) seq_len(baseline) else baseline
    if (any(base < 1L | base > n))
        stop("baseline indices fall outside the series")
    base <- sort(unique(as.integer(base)))
    xb <- x[base]
    ok <- !is.na(xb)
    pair <- which(ok[-length(ok)] & ok[-1L] & diff(base) == 1L)
    xb <- xb[ok]
    if (length(xb) < 3L)
        stop("the baseline needs at least three non-missing runs")
    if (length(xb) < 20L)
        warning("the baseline has fewer than 20 non-missing runs; ",
                "its limits are uncertain")
    mu <- mean(xb)
    if (sigma == "mr" && !length(pair))
        stop("the baseline has no two adjacent observed runs for a ",
             "moving range")
    s <- if (sigma == "mr")
        mean(abs(x[base][pair + 1L] - x[base][pair])) / 1.128 else
        stats::sd(xb)
    if (!is.finite(s) || s <= 0)
        stop("the baseline has no spread; cannot set limits")
    out <- if (method == "shewhart") .shewhart(x, mu, s, L) else
        .afterBaseline(x, max(base), function(y) switch(method,
            ewma = .ewma(y, mu, s, lambda, L),
            cusum = .cusum(y, mu, s, k, h)))
    res <- data.frame(run = runs, value = x, out,
                      baseline = seq_len(n) %in% base,
                      stringsAsFactors = FALSE)
    attr(res, "method") <- method
    attr(res, "params") <- switch(method,
        shewhart = c(L = L), ewma = c(lambda = lambda, L = L),
        cusum = c(k = k, h = h))
    attr(res, "sigma") <- s
    res
}

# Runs the chart on the runs after the last baseline run only; the rows up
# to it keep the column layout with NA statistics and no signal.
.afterBaseline <- function(x, last, chart) {
    n <- length(x)
    shape <- chart(x[1L])
    pad <- shape[rep(1L, last), , drop = FALSE]
    pad[] <- lapply(pad, function(col)
        if (is.logical(col)) rep(FALSE, last) else rep(NA_real_, last))
    pad$center <- rep(shape$center[1L], last)
    if (last >= n)
        return(pad)
    out <- rbind(pad, chart(x[(last + 1L):n]))
    rownames(out) <- NULL
    out
}

.shewhart <- function(x, mu, s, L) {
    data.frame(stat = x, center = mu, lower = mu - L * s,
               upper = mu + L * s,
               signal = !is.na(x) & abs(x - mu) > L * s)
}

.ewma <- function(x, mu, s, lambda, L) {
    if (lambda <= 0 || lambda > 1)
        stop("'lambda' must be in (0, 1]")
    n <- length(x)
    z <- numeric(n)
    w <- numeric(n)
    prev <- mu
    t <- 0
    for (i in seq_len(n)) {
        if (!is.na(x[i])) {
            t <- t + 1
            prev <- lambda * x[i] + (1 - lambda) * prev
        }
        z[i] <- prev
        w[i] <- L * s * sqrt(lambda / (2 - lambda) *
                             (1 - (1 - lambda)^(2 * max(t, 1))))
    }
    data.frame(stat = z, center = mu, lower = mu - w, upper = mu + w,
               signal = !is.na(x) & abs(z - mu) > w)
}

.cusum <- function(x, mu, s, k, h) {
    n <- length(x)
    up <- numeric(n)
    lo <- numeric(n)
    cu <- 0
    cl <- 0
    for (i in seq_len(n)) {
        if (!is.na(x[i])) {
            cu <- max(0, cu + (x[i] - mu) / s - k)
            cl <- max(0, cl - (x[i] - mu) / s - k)
        }
        up[i] <- cu
        lo[i] <- cl
    }
    data.frame(stat = pmax(up, lo), center = 0, lower = 0, upper = h,
               cusumUp = up, cusumLow = lo,
               signal = !is.na(x) & pmax(up, lo) > h)
}

#' First run that signals on a control chart
#'
#' @param chart A chart from [driftChart()].
#' @param after Ignore signals at or before this run index. `NULL` (the
#'   default) uses the last baseline run.
#'
#' @return The index of the first signalling run, or `NA` if none.
#'
#' @examples
#' ch <- driftChart(c(rep(c(1, 1.01, 0.99, 1.02), 5), 0.5),
#'                  baseline = 20, method = "shewhart")
#' firstSignal(ch)
#' @export
firstSignal <- function(chart, after = NULL) {
    if (!any(chart$baseline))
        stop("the chart has no baseline runs")
    if (is.null(after))
        after <- max(which(chart$baseline))
    hit <- which(chart$signal & seq_len(nrow(chart)) > after)
    if (length(hit)) hit[1] else NA_integer_
}
