# spikeWatch

spikeWatch scores every run against the known amounts of its spike-in controls, and charts those scores across runs so that a run that drifts is flagged. It is an R package built on SummarizedExperiment.

The known amounts sit in an assay beside the measurements, cell by cell. That one idea covers feature-level designs, where each spike has its own amount (the ERCC RNA mixes), and sample-level designs, where the amount changes from sample to sample (a peptide or UPS1 dilution series). Per run, spikeWatch reports the dose-response slope and fit, the residual spread, the overall response and the detection limit. Between two groups it reports fold-change accuracy against the designed ratios. Across runs it puts any of those figures on a Shewhart, EWMA or CUSUM chart with limits from a baseline of good runs.

## Install

```r
# not yet on Bioconductor; from GitHub:
BiocManager::install("Daneshpajouh/spikeWatch")
```

R 4.6 or later. Imports SummarizedExperiment, S4Vectors and ggplot2.

## Use

```r
library(spikeWatch)
se <- addTruth(se, erccDesign(), mix = "mix")   # colData column "mix" says Mix1 or Mix2
m <- runMetrics(se, assay = "counts", minAmount = 1)
ch <- driftChart(m, "slope", baseline = 20, method = "cusum")
firstSignal(ch)
plotDrift(ch)
```

For a dilution series, give the design one amount per spike and each sample a dose column: `addTruth(se, design, dose = "fmol")`. `simulateSpikeRuns()` makes series with a planted drift, to test a chart and its limits before trusting them. The vignette walks through both designs.

`minAmount` matters. Spikes near the detection limit are seen only when their count happens to be high, so fitting every detected spike pulls the slope down and ties it to depth. Choosing spikes by their known amount, not by what was measured, avoids that. In 200 simulated runs with no drift (`validation/v0_min_amount.R`; true slope 1):

<!-- V0 -->
| minAmount | Mean slope | SD of slope | Correlation with depth |
|---|---|---|---|
| 0 | 0.985 | 0.0150 | 0.226 |
| 1 | 1.003 | 0.0102 | -0.182 |
| 10 | 1.001 | 0.0123 | -0.038 |
<!-- /V0 -->

A floor of 1 attomol/uL fixes the mean and narrows the spread but leaves a weaker link to depth, of the opposite sign; a floor of 10 removes that too, at some cost in spread. The validation below uses 1, as frozen in its plan before this table was made.

## Neighbours

Bioconductor already has spike-in and QC packages, read from the release index on 2026-10-27:

- erccdashboard: "Technical performance metrics for differential gene expression experiments using External RNA Controls Consortium (ERCC) spike-in ratio mixtures."
- Anaquin: analysis of sequins, "synthetic sequencing spike-in controls".
- MSstatsQC: "longitudinal system suitability monitoring and quality control tools for proteomic experiments."
- qcmetrics: "a framework for generic quality control of data."

spikeWatch takes a different shape: one data model for RNA and proteomics spike-ins, figures computed against the known amounts for every run, and charts of those figures over time. Its ERCC design table comes from erccdashboard.

## Validation

The plan, the detectors and the rule for claiming a gain were committed before any result was computed (`validation/PLAN.md`). Every figure below is filled in from `validation/out/` by `validation/render_readme.py`.

### V1. Planted drift in simulated runs

200 series per scenario, 40 runs each, ERCC Mix 1, sequencing depth changing from run to run. Runs 1 to 20 are the baseline; the drift starts at run 31. Every chart is a two-sided CUSUM (k = 0.5, h = 5) with limits from the baseline. The truth-based charts are spikeWatch's `slope`, `residSD` and `offset`. The two charts without truth use only what a lab can compute from spike counts when it knows which rows are spikes but not their amounts: B1 is the log2 total spike count and B2 is the spread of log2 spike counts, which also shrinks under compression. "Detected" counts series with no false alarm before run 31. All scenarios share seeds 1 to 200, so runs before 31 are identical across scenarios and each chart's false alarm rate repeats by design: the comparison is paired.

<!-- V1 -->
| Scenario | Chart | False alarm, runs 21 to 30 | Detected by run 40 | Median delay (runs) |
|---|---|---|---|---|
| null | slope | 8.5% |  |  |
| null | residSD | 12.0% |  |  |
| null | offset | 10.5% |  |  |
| null | B1_total | 9.5% |  |  |
| null | B2_spread | 7.0% |  |  |
| compression small | slope | 8.5% | 100.0% | 1.0 |
| compression small | residSD | 12.0% | 9.7% | 4.0 |
| compression small | offset | 10.5% | 21.8% | 5.0 |
| compression small | B1_total | 9.5% | 65.2% | 5.0 |
| compression small | B2_spread | 7.0% | 95.7% | 3.0 |
| compression | slope | 8.5% | 100.0% | 0.0 |
| compression | residSD | 12.0% | 8.5% | 5.0 |
| compression | offset | 10.5% | 36.9% | 5.0 |
| compression | B1_total | 9.5% | 100.0% | 2.0 |
| compression | B2_spread | 7.0% | 100.0% | 1.0 |
| sensitivity loss | slope | 8.5% | 18.0% | 4.0 |
| sensitivity loss | residSD | 12.0% | 77.3% | 3.5 |
| sensitivity loss | offset | 10.5% | 100.0% | 2.0 |
| sensitivity loss | B1_total | 9.5% | 100.0% | 2.0 |
| sensitivity loss | B2_spread | 7.0% | 94.6% | 3.0 |
| precision loss | slope | 8.5% | 36.6% | 4.0 |
| precision loss | residSD | 12.0% | 100.0% | 0.0 |
| precision loss | offset | 10.5% | 11.7% | 2.0 |
| precision loss | B1_total | 9.5% | 11.6% | 5.0 |
| precision loss | B2_spread | 7.0% | 15.6% | 6.0 |
<!-- /V1 -->

What this shows, under the rule frozen in the plan (a gain is claimed only when the truth-based chart detects more often and its false alarm rate is no higher):

- Small compression: the slope chart detects it more often than total counts (B1), at a false alarm rate no higher. That gain over B1 is claimed. Against the spread of counts (B2) the slope chart also detects more often and sooner, but its false alarm rate was higher, so no gain over B2 is claimed.
- Larger compression: every chart except offset and residSD detects it in every clean series, so no chart detects more often and no gain is claimed. The slope chart's median delay is the shortest, but delay is not part of the frozen rule.
- Precision loss: the residual spread chart detects it in every clean series, while the two charts without truth stay near their false alarm rates. Its own false alarm rate was higher than B1's and B2's, so under the frozen rule no gain is claimed, though the gap in detection is large.
- Sensitivity loss: the offset chart and total counts do equally well. Known amounts add nothing here.
- False alarms: every chart raised one within 10 runs after the baseline in 7 to 12 per cent of series. That is the cost of limits estimated from 20 runs. A longer baseline or a larger `h` lowers it; this was not tested here.

### V2. Real RNA-seq: SEQC ERCC counts

Lab 5 of the SEQC interlaboratory study, as shipped in erccdashboard: 4 runs of Universal Human Reference RNA with ERCC Mix 1 and 4 runs of Human Brain Reference RNA with Mix 2 (erccdashboard vignette). Per run, on counts per kilobase of ERCC length, with `minAmount = 1`; then fold-change accuracy on raw counts, centred on the spikes.

<!-- V2 -->
| Run | Detected | Slope | R squared | Detection limit (attomol/uL) |
|---|---|---|---|---|
| UHRR_1 | 83 of 92 | 1.010 | 0.952 | 1.079 |
| UHRR_2 | 86 of 92 | 1.016 | 0.953 | 0.593 |
| UHRR_3 | 87 of 92 | 1.014 | 0.959 | 0.666 |
| UHRR_4 | 90 of 92 | 1.017 | 0.950 | 0.119 |
| HBRR_1 | 89 of 92 | 1.004 | 0.959 | 0.068 |
| HBRR_2 | 89 of 92 | 1.002 | 0.947 | 0.238 |
| HBRR_3 | 87 of 92 | 1.007 | 0.958 | 0.264 |
| HBRR_4 | 88 of 92 | 1.010 | 0.953 | 0.123 |

| Expected log2 ratio | Spikes | Median error (log2) | RMSE (log2) |
|---|---|---|---|
| 2.000 | 18 | -0.030 | 0.257 |
| 0.000 | 18 | 0.016 | 0.161 |
| -0.585 | 21 | -0.016 | 0.279 |
| -1.000 | 21 | 0.011 | 0.300 |
<!-- /V2 -->

The slopes sit close to 1 and the measured ratios match the designed 4:1, 1:1, 1:1.5 and 1:2 with small median errors. The detection limit follows depth: the run with the fewest reads has the highest limit.

### V3. Real targeted proteomics: msqc1

The msqc1 data package: 14 heavy labelled peptides at known amounts per vial, measured in 8 replicate runs on each of 5 instruments, and a dilution series of the same mix. Peptide areas are summed over transitions. A peptide missing from an instrument's method is left out for that instrument; the first pass counted it as undetected, and its log is kept (`validation/v3_first_pass.log`).

Within a run, across the 14 peptides:

<!-- V3 -->
| Instrument | Runs | Slope (median) | Slope range | R squared (median) |
|---|---|---|---|---|
| QExactive | 8 | 0.784 | 0.779 to 0.801 | 0.852 |
| QExactiveHF | 8 | 0.747 | 0.715 to 0.959 | 0.660 |
| QTRAP | 8 | 0.788 | 0.780 to 0.794 | 0.580 |
| TRIPLETOF | 8 | 0.718 | 0.704 to 0.778 | 0.876 |
| TSQVantage | 8 | 0.996 | 0.983 to 1.002 | 0.835 |

| Instrument | Runs | Doses | Slope against dose | 95% interval |
|---|---|---|---|---|
| QExactive | 18 | 6 | 0.903 | 0.836 to 0.970 |
| QExactiveHF | 18 | 6 | 1.009 | 0.993 to 1.025 |
| QTRAP | 18 | 6 | 0.943 | 0.908 to 0.978 |
| TRIPLETOF | 18 | 6 | 0.862 | 0.834 to 0.891 |
| TSQVantage | 18 | 6 | 1.066 | 1.031 to 1.101 |
<!-- /V3 -->

The first table is not a clean linearity test. Each peptide ionises differently, so a slope across different peptides mixes response factors with amounts; four instruments gave slopes well below 1. The dilution series (second table) changes the dose of the same peptides, and there the response follows dose with slopes near 1 on every instrument, some intervals including 1 and some not. For designs where each spike is a different molecule, chart the per-run figures against their own baseline rather than reading the slope as absolute linearity.

## Checks

`R CMD check`, BiocCheck and the testthat suite run on every push (`.github/workflows/check.yaml`). Tests use small fixtures whose answers can be worked by hand: an exact power law, a hand-computed ratio table, and drift series whose chart values follow from the formulas.

## Licence

GPL-3. The ERCC design table is taken from erccdashboard (GPL (>= 2)); `inst/script/make_ercc_mixes.R` records how and the checks run on it.
