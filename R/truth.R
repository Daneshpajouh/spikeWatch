#' Attach the known spike-in amounts to an experiment
#'
#' Adds an assay of expected amounts, the same shape as the measured
#' assays, so that every later step compares a measured value with its
#' known truth cell by cell. The expected amount of feature `f` in sample
#' `s` is the design amount of `f` in the mix of `s`, times the dose of
#' `s`. Rows that are not spike-ins hold `NA`.
#'
#' @param se A `SummarizedExperiment`. Row names must hold the feature
#'   identifiers used in the design.
#' @param design A design from [spikeDesign()] or [erccDesign()].
#' @param mix Name of the `colData` column giving each sample's mix. May
#'   be `NULL` when the design has one mix.
#' @param dose Optional name of a numeric `colData` column that scales the
#'   design amounts per sample, for example the fmol of UPS1 spiked into
#'   each sample, or a dilution factor.
#' @param name Name of the new assay.
#'
#' @return `se` with an added assay `name`, and `rowData` columns `spike`
#'   (logical) and `spikeGroup`.
#'
#' @examples
#' library(SummarizedExperiment)
#' d <- spikeDesign(c("s1", "s2"), amount = c(1, 10))
#' m <- matrix(c(2, 20, 5, 4, 41, 9), nrow = 3,
#'             dimnames = list(c("s1", "s2", "gene"), c("r1", "r2")))
#' se <- SummarizedExperiment(list(counts = m))
#' se <- addTruth(se, d)
#' assay(se, "expected")
#' @export
#' @importFrom SummarizedExperiment assay assay<- colData rowData rowData<-
addTruth <- function(se, design, mix = NULL, dose = NULL,
                     name = "expected") {
    if (!methods::is(se, "SummarizedExperiment"))
        stop("'se' must be a SummarizedExperiment")
    if (is.null(rownames(se)))
        stop("'se' needs row names that match the design features")
    mixes <- unique(design$mix)
    cd <- colData(se)
    if (is.null(mix)) {
        if (length(mixes) != 1L)
            stop("the design has ", length(mixes),
                 " mixes; name the colData column with 'mix'")
        sampleMix <- rep(mixes, ncol(se))
    } else {
        if (!mix %in% colnames(cd))
            stop("colData has no column '", mix, "'")
        sampleMix <- as.character(cd[[mix]])
        bad <- setdiff(unique(sampleMix[!is.na(sampleMix)]), mixes)
        if (length(bad))
            stop("mix labels not in the design: ",
                 paste(bad, collapse = ", "))
    }
    sampleDose <- rep(1, ncol(se))
    if (!is.null(dose)) {
        if (!dose %in% colnames(cd))
            stop("colData has no column '", dose, "'")
        sampleDose <- as.numeric(cd[[dose]])
        if (any(!is.na(sampleDose) & sampleDose < 0))
            stop("doses must not be negative")
    }
    rows <- match(design$feature, rownames(se))
    if (all(is.na(rows)))
        stop("no design feature matches a row name of 'se'")
    miss <- unique(design$feature[is.na(rows)])
    if (length(miss))
        warning(length(miss), " design feature(s) not in the row names of ",
                "'se', so they get no expected amount: ",
                paste(miss[seq_len(min(5L, length(miss)))], collapse = ", "),
                if (length(miss) > 5L) ", ..." else "")
    expected <- matrix(NA_real_, nrow(se), ncol(se),
                       dimnames = dimnames(se))
    for (m in mixes) {
        d <- design[design$mix == m & !is.na(rows), , drop = FALSE]
        cols <- which(sampleMix == m)
        if (!length(cols) || !nrow(d))
            next
        r <- match(d$feature, rownames(se))
        expected[r, cols] <- outer(d$amount, sampleDose[cols])
    }
    assay(se, name, withDimnames = FALSE) <- expected
    group <- design$group[match(rownames(se), design$feature)]
    rowData(se)$spike <- rownames(se) %in% design$feature
    rowData(se)$spikeGroup <- group
    se
}
