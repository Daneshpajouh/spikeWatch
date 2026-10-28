#' Build a spike-in design
#'
#' A design lists every spike-in feature and the amount of it in each mix.
#' Amounts are in one unit for the whole design. A design with one mix
#' suits equimolar spikes such as UPS1, where the amount per sample comes
#' from a dose column in the sample data (see [addTruth()]).
#'
#' @param feature Character. Feature identifiers, matched against the row
#'   names of the experiment.
#' @param amount Numeric. Known amount of each feature in its mix. Must be
#'   positive.
#' @param mix Character. Mix label for each row. Defaults to one mix
#'   called `"mix"`.
#' @param group Optional character. A label for each feature, such as the
#'   ERCC ratio subgroup. Carried into the results.
#' @param unit Character. The unit of `amount`.
#'
#' @return A `data.frame` with columns `feature`, `mix`, `amount` and
#'   `group`, and an attribute `unit`.
#'
#' @examples
#' d <- spikeDesign(c("s1", "s2", "s3"), amount = c(1, 10, 100),
#'                  unit = "fmol")
#' d
#' @export
spikeDesign <- function(feature, amount, mix = "mix", group = NA_character_,
                        unit = "amount") {
    feature <- as.character(feature)
    n <- length(feature)
    if (n == 0L)
        stop("a design needs at least one feature")
    if (!is.numeric(amount) || length(amount) != n)
        stop("'amount' must be numeric with one value per feature")
    if (anyNA(amount) || any(amount <= 0))
        stop("every amount must be positive")
    mix <- rep_len(as.character(mix), n)
    group <- rep_len(as.character(group), n)
    d <- data.frame(feature = feature, mix = mix, amount = amount,
                    group = group, stringsAsFactors = FALSE)
    if (anyDuplicated(d[, c("feature", "mix")]))
        stop("each feature may appear only once per mix")
    attr(d, "unit") <- unit
    d
}

#' The ERCC RNA spike-in mixes
#'
#' Returns the design of the ERCC ExFold RNA Spike-In Mixes 1 and 2: 92
#' transcripts in four subgroups whose Mix 1 to Mix 2 ratios are 4:1, 1:1,
#' 1:1.5 and 1:2; `group` holds that label. Concentrations are in
#' attomoles per microlitre. The table
#' was taken from the Bioconductor package erccdashboard; see
#' `inst/script/make_ercc_mixes.R` for its source and the checks run on it.
#'
#' @return A design as returned by [spikeDesign()], with mixes `"Mix1"` and
#'   `"Mix2"`.
#'
#' @examples
#' ercc <- erccDesign()
#' table(ercc$group, ercc$mix)
#' @export
erccDesign <- function() {
    f <- system.file("extdata", "ercc_mixes.tsv", package = "spikeWatch")
    t <- utils::read.delim(f, stringsAsFactors = FALSE)
    spikeDesign(rep(t$feature, 2),
                amount = c(t$Mix1, t$Mix2),
                mix = rep(c("Mix1", "Mix2"), each = nrow(t)),
                group = rep(t$subgroup, 2),
                unit = "attomol/uL")
}
