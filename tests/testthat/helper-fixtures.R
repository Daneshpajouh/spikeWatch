# A tiny experiment: three spikes and one endogenous gene, two samples.
tinySE <- function(m) {
    SummarizedExperiment::SummarizedExperiment(
        list(counts = m),
        colData = S4Vectors::DataFrame(mix = c("A", "B"), dose = c(1, 2),
                                       row.names = colnames(m)))
}
