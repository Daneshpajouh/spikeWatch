# V3: real targeted proteomics, msqc1. PLAN.md.
suppressMessages({
    library(spikeWatch)
    library(SummarizedExperiment)
    library(msqc1)
})
pep <- peptides
design <- spikeDesign(as.character(pep$Peptide.Sequence),
                      amount = pep$SIL.peptide.per.vial,
                      unit = "per vial [unit not stated in msqc1]")

# Heavy areas summed over transitions, one column per run.
toSE <- function(x) {
    h <- x[x$Isotope.Label.Type == "heavy" &
           x$Peptide.Sequence %in% design$feature, ]
    h$run <- paste(h$instrument, h$Replicate.Name, sep = "|")
    a <- tapply(h$Area, list(as.character(h$Peptide.Sequence), h$run),
                sum, na.rm = TRUE)
    runs <- colnames(a)
    inst <- sub("\\|.*", "", runs)
    dose <- tapply(h$relative.amount, h$run, function(v) unique(v)[1])
    cd <- S4Vectors::DataFrame(instrument = inst,
                               dose = as.numeric(dose[runs]),
                               row.names = runs)
    SummarizedExperiment(list(area = a), colData = cd)
}

dir.create("validation/out", showWarnings = FALSE)

# 8 replicates per instrument: per-run dose response.
se8 <- addTruth(toSE(msqc1_8rep), design)
m8 <- runMetrics(se8, assay = "area")
m8$instrument <- se8$instrument
write.csv(m8, "validation/out/v3_8rep_run_metrics.csv", row.names = FALSE)
agg <- do.call(rbind, lapply(split(m8, m8$instrument), function(d)
    data.frame(instrument = d$instrument[1], runs = nrow(d),
               slopeMedian = stats::median(d$slope),
               slopeMin = min(d$slope), slopeMax = max(d$slope),
               r2Median = stats::median(d$r2),
               residSDMedian = stats::median(d$residSD),
               lodMedian = stats::median(d$lod),
               lodCensoredRuns = sum(d$lodCensored))))
write.csv(agg, "validation/out/v3_8rep_by_instrument.csv", row.names = FALSE)
cat("msqc1_8rep, per-run dose response by instrument:\n")
print(agg, digits = 3, row.names = FALSE)

# Dilution series: the sample-level dose enters the truth.
sed <- toSE(msqc1_dil)
sed <- sed[, !is.na(sed$dose) & sed$dose > 0]
sed <- addTruth(sed, design, dose = "dose")
md <- runMetrics(sed, assay = "area")
md$instrument <- sed$instrument
md$dose <- sed$dose
# offset is log2 measured minus log2 expected, and expected already holds
# the dose; so add log2 dose back to test linearity in dose.
md$log2Response <- md$offset + log2(md$dose)
write.csv(md, "validation/out/v3_dil_run_metrics.csv", row.names = FALSE)
lin <- do.call(rbind, lapply(split(md, md$instrument), function(d) {
    f <- stats::lm(log2Response ~ log2(dose), data = d)
    ci <- stats::confint(f)[2, ]
    data.frame(instrument = d$instrument[1], runs = nrow(d),
               doses = length(unique(d$dose)),
               doseSlope = unname(stats::coef(f)[2]),
               ciLow = unname(ci[1]), ciHigh = unname(ci[2]))
}))
write.csv(lin, "validation/out/v3_dil_linearity.csv", row.names = FALSE)
cat("\nmsqc1_dil, per-run response against known dose (slope 1 = linear):\n")
print(lin, digits = 3, row.names = FALSE)
