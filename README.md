# spikeWatch

spikeWatch scores every run against the known amounts of its spike-in controls, and charts those scores across runs to show drift. The charts are exploratory, not validated run-release alarms: at the default limits with a 20-run baseline, 15.7% to 35.2% of simulated series with no drift raised a false alarm within 20 monitored runs (Default limits, below), so set and check the limits on your own runs before acting on an alarm. It is an R package built on SummarizedExperiment.

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
m <- runMetrics(se, assay = "counts", baseline = 20)   # fit range chosen from runs 1 to 20
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

A floor of 1 attomol/uL fixes the mean and narrows the spread but leaves a weaker link to depth, of the opposite sign; a floor of 10 removes that too, at some cost in spread. The default, `minAmount = "auto"`, takes the smallest amount at and above which every spike was seen in every run named by `baseline`, so the range is chosen from good runs and is the same for every run. Without `baseline` it uses every run, which suits one experiment scored once; for runs you will chart, name the good ones, or one failed run will raise the range for all. Version 1 of the validation used a fixed 1, as frozen in its plan; version 2 uses the default.

## Neighbours

Read at source on 2026-10-27:

- erccdashboard (Bioconductor, installed copy 1.46.0): "Technical performance metrics for differential gene expression experiments using External RNA Controls Consortium (ERCC) spike-in ratio mixtures." It scores one experiment, between two samples (AUC, LODR, ratio bias).
- MSstatsQC (Bioconductor release 3.23, vignette): "longitudinal system suitability monitoring and quality control tools for proteomic experiments." It charts per-peptide run metrics such as retention time and peak area with X, moving range and CUSUM charts and change point estimation, with limits from a guide set of good runs or from a known mean and SD; its vignette recommends a CUSUM reference value of 0.5 and decision interval of 5 for standardized data.
- Panorama AutoQC (Bereman and others, Journal of Proteome Research 15(12), 4763-4769, 2016, doi:10.1021/acs.jproteome.6b00744): an automated pipeline that extracts "identification-free metrics" from targeted peptides and shows them with "Levey-Jennings and Pareto plots" across instruments.
- Anaquin (Bioconductor): analysis of sequins, "synthetic sequencing spike-in controls".

Truth-based QC and control charts are not new. What spikeWatch adds is narrower: one R object model in which every run is scored against the known spike amounts, for RNA and proteomics designs alike. Its validation compares it with the best charts we could build without the known amounts. In simulation, under the frozen rule, the amounts gave one gain that was not confirmed on fresh seeds (EWMA, compression -0.05) and no better alarm elsewhere; what they add is per-run figures in units a lab can read. Its ERCC design table comes from erccdashboard.

## Validation

Two rounds, each with its plan, detectors and rule for claiming a gain committed before any of its results were computed. Every figure is filled in from the output files by a script; none is typed.

## Validation, version 2 (current code)

These results were computed at commit d7dfeb0. That commit also held `calibrateChart()`, which is not in this release (see "Setting limits from a short baseline" below); to rerun `validation/pv2_2_score.R` as it ran, use the `calibrate-chart` branch.

Plan: `validation/PLAN_v2.md`, committed before any v2 number, with one dated addendum that corrected a sentence and changed no rule. Scripts `validation/pv2_*.R`; outputs `validation/out-v2/`; tables by `validation/render_readme_v2.py`.

**Setup.** Simulated series of 40 runs of ERCC Mix 1, with run depth changing from run to run. Runs 1 to 20 are the baseline; runs 21 to 40 are monitored, and any planted drift starts at run 21. Six simulation families: the base model (F0), and five departures from the model the fit assumes: compression about the top amount (F1), overdispersed counts (F2), saturation of the highest spikes (F3), fixed per-spike biases (F4) and correlated runs (F5). Five drifts: compression (two sizes, and a ramp), sensitivity loss and precision loss.

**What is compared.** The truth-based figures (`slope`, `residSD`, `offset`, with the fit range chosen from the baseline) against the best charts we could build without the known amounts: B1, the log2 total spike count; B2, the spread of log2 spike counts; and B3, the profile slope from `profileMetrics()`, which compares each run with the average profile of the same spikes over the baseline. U-T and U-F are the unions (slope or residSD; B2 or B3). Every detector is calibrated on separate null series to the same 5% chance of any false alarm within 10 or 20 monitored runs, so a detector cannot look sensitive by alarming more. This is an oracle calibration, from 1,000 null series per family, which no lab can do: it makes the comparison fair between detectors, not the limits usable in practice. Each comparison is a paired bootstrap interval on the difference in detection, over 400 series per drift and family.

**The gain rule.** A gain or loss is claimed for a pair, drift and chart only when the interval excludes zero at both horizons in at least five of the six families.

<!-- PV2-HEADLINE -->
| Pair (with truth vs without) | Drift | CUSUM | EWMA |
|---|---|---|---|
| U-T vs U-F | compression -0.05 | none (gain in 1, loss in 1) | none (gain in 1, loss in 1) |
| U-T vs U-F | compression -0.10 | none (gain in 0, loss in 0) | none (gain in 0, loss in 0) |
| U-T vs U-F | compression ramp | none (gain in 0, loss in 0) | none (gain in 0, loss in 0) |
| U-T vs U-F | sensitivity loss | **loss** (gain in 0, loss in 5) | **loss** (gain in 0, loss in 5) |
| U-T vs U-F | precision loss | **gain** (gain in 6, loss in 0) | **gain** (gain in 6, loss in 0) |
| slope vs B3 | compression -0.05 | none (gain in 0, loss in 1) | none (gain in 1, loss in 1) |
| slope vs B3 | compression -0.10 | none (gain in 0, loss in 0) | none (gain in 0, loss in 0) |
| slope vs B3 | compression ramp | none (gain in 0, loss in 0) | none (gain in 0, loss in 0) |
| slope vs B3 | sensitivity loss | none (gain in 3, loss in 0) | none (gain in 2, loss in 0) |
| slope vs B3 | precision loss | none (gain in 3, loss in 0) | none (gain in 4, loss in 0) |
| residSD vs B2 | compression -0.05 | **loss** (gain in 0, loss in 6) | **loss** (gain in 0, loss in 6) |
| residSD vs B2 | compression -0.10 | **loss** (gain in 0, loss in 5) | **loss** (gain in 0, loss in 6) |
| residSD vs B2 | compression ramp | **loss** (gain in 0, loss in 6) | **loss** (gain in 0, loss in 6) |
| residSD vs B2 | sensitivity loss | **loss** (gain in 0, loss in 5) | **loss** (gain in 0, loss in 5) |
| residSD vs B2 | precision loss | **gain** (gain in 6, loss in 0) | **gain** (gain in 6, loss in 0) |
| offset vs B1 | compression -0.05 | **loss** (gain in 1, loss in 5) | **loss** (gain in 1, loss in 5) |
| offset vs B1 | compression -0.10 | **loss** (gain in 1, loss in 5) | **loss** (gain in 1, loss in 5) |
| offset vs B1 | compression ramp | **loss** (gain in 1, loss in 5) | **loss** (gain in 1, loss in 5) |
| offset vs B1 | sensitivity loss | none (gain in 1, loss in 0) | none (gain in 2, loss in 0) |
| offset vs B1 | precision loss | none (gain in 0, loss in 2) | none (gain in 0, loss in 2) |

Each gain the rule claims in this table, with the range of each side's false-alarm rate on fresh null seeds over the six families and two horizons, and the number of those cells where either side's false-alarm interval excludes the 5% target:

| Pair | Drift | Chart | False alarm, with truth | False alarm, without | Cells off the 5% target |
|---|---|---|---|---|---|
| U-T vs U-F | precision loss | CUSUM | 3.6% to 6.9% | 3.2% to 6.1% | 3 of 12 |
| residSD vs B2 | precision loss | CUSUM | 4.6% to 6.9% | 4.2% to 6.5% | 4 of 12 |
| U-T vs U-F | precision loss | EWMA | 3.7% to 8.2% | 3.9% to 6.4% | 3 of 12 |
| residSD vs B2 | precision loss | EWMA | 3.8% to 7.2% | 4.2% to 6.9% | 5 of 12 |
<!-- /PV2-HEADLINE -->

**Superseded in part.** The precision-loss gains in this table (U-T vs U-F, residSD vs B2) do not hold against the truth-free residual spread `profileResidSD`; see addendum 2 below.

**F1 counted twice.** A review pointed out that F1 differs from F0 only in the pivot, which cancels when the slope is 1, so for sensitivity and precision loss F1 repeats F0 and the five-of-six rule counts that series twice. The rule stays as planned; here is each claim for those two drifts with F1 left out:

<!-- PV2-DEDUP -->
F1 had the same detection rates as F0 in 16 of 16 sensitivity-loss cells and 16 of 16 precision-loss cells.

| Pair | Drift | Chart | Claim | Families with F1 left out |
|---|---|---|---|---|
| U-T vs U-F | sensitivity loss | CUSUM | loss | 4 of 5 |
| U-T vs U-F | sensitivity loss | EWMA | loss | 4 of 5 |
| U-T vs U-F | precision loss | CUSUM | gain | 5 of 5 |
| U-T vs U-F | precision loss | EWMA | gain | 5 of 5 |
| residSD vs B2 | sensitivity loss | CUSUM | loss | 4 of 5 |
| residSD vs B2 | sensitivity loss | EWMA | loss | 4 of 5 |
| residSD vs B2 | precision loss | CUSUM | gain | 5 of 5 |
| residSD vs B2 | precision loss | EWMA | gain | 5 of 5 |
<!-- /PV2-DEDUP -->

**Addendum 2: a truth-free residual comparator.** Both reviewers of the v2 results pointed out that the precision-loss gain was scored against B2, while `profileMetrics()` already gives `profileResidSD`, the residual spread of each run's fit on the baseline profile, which needs no amounts. `validation/PLAN_v2_ADDENDUM_2.md` was committed before any of the following was computed. It adds B4 (`profileResidSD`), U-F2 (B3 or B4), C5 (`residSD` vs B4) and C6 (U-T vs U-F2), under the same seeds, calibration and gain rule, and one confirmation on untouched seeds. The series were simulated again and their six earlier figures checked equal to the saved ones, and C1 to C4 reproduced exactly, before C5 and C6 were scored.

<!-- PV2-A2 -->
| Pair (with truth vs without) | Drift | CUSUM | EWMA |
|---|---|---|---|
| residSD vs B2 | compression -0.05 | **loss** (gain in 0, loss in 6) | **loss** (gain in 0, loss in 6) |
| residSD vs B2 | compression -0.10 | **loss** (gain in 0, loss in 5) | **loss** (gain in 0, loss in 6) |
| residSD vs B2 | compression ramp | **loss** (gain in 0, loss in 6) | **loss** (gain in 0, loss in 6) |
| residSD vs B2 | sensitivity loss | **loss** (gain in 0, loss in 5) | **loss** (gain in 0, loss in 5) |
| residSD vs B2 | precision loss | **gain** (gain in 6, loss in 0) | **gain** (gain in 6, loss in 0) |
| residSD vs B4 | compression -0.05 | none (gain in 3, loss in 0) | none (gain in 4, loss in 0) |
| residSD vs B4 | compression -0.10 | none (gain in 4, loss in 0) | none (gain in 4, loss in 0) |
| residSD vs B4 | compression ramp | none (gain in 4, loss in 0) | none (gain in 2, loss in 0) |
| residSD vs B4 | sensitivity loss | none (gain in 0, loss in 3) | none (gain in 2, loss in 3) |
| residSD vs B4 | precision loss | none (gain in 0, loss in 0) | none (gain in 0, loss in 0) |
| U-T vs U-F | compression -0.05 | none (gain in 1, loss in 1) | none (gain in 1, loss in 1) |
| U-T vs U-F | compression -0.10 | none (gain in 0, loss in 0) | none (gain in 0, loss in 0) |
| U-T vs U-F | compression ramp | none (gain in 0, loss in 0) | none (gain in 0, loss in 0) |
| U-T vs U-F | sensitivity loss | **loss** (gain in 0, loss in 5) | **loss** (gain in 0, loss in 5) |
| U-T vs U-F | precision loss | **gain** (gain in 6, loss in 0) | **gain** (gain in 6, loss in 0) |
| U-T vs U-F2 | compression -0.05 | none (gain in 3, loss in 0) | **gain** (gain in 6, loss in 0) |
| U-T vs U-F2 | compression -0.10 | none (gain in 0, loss in 0) | none (gain in 1, loss in 0) |
| U-T vs U-F2 | compression ramp | none (gain in 0, loss in 0) | none (gain in 1, loss in 0) |
| U-T vs U-F2 | sensitivity loss | none (gain in 0, loss in 3) | none (gain in 1, loss in 3) |
| U-T vs U-F2 | precision loss | none (gain in 0, loss in 0) | none (gain in 0, loss in 0) |

Each gain the rule claims here, with the range of each side's false-alarm rate on fresh null seeds over the six families and two horizons, and the number of those cells where either side's false-alarm interval excludes the 5% target:

| Pair | Drift | Chart | False alarm, with truth | False alarm, without | Cells off the 5% target |
|---|---|---|---|---|---|
| U-T vs U-F | precision loss | CUSUM | 3.6% to 6.9% | 3.2% to 6.1% | 3 of 12 |
| residSD vs B2 | precision loss | CUSUM | 4.6% to 6.9% | 4.2% to 6.5% | 4 of 12 |
| U-T vs U-F2 | compression -0.05 | EWMA | 3.7% to 8.2% | 2.7% to 6.2% | 6 of 12 |
| U-T vs U-F | precision loss | EWMA | 3.7% to 8.2% | 3.9% to 6.4% | 3 of 12 |
| residSD vs B2 | precision loss | EWMA | 3.8% to 7.2% | 4.2% to 6.9% | 5 of 12 |

Cells of the claimed gain on the addendum's own pairs (C5, C6), each family and horizon:

| Pair | Drift | Chart | Family | Runs | Detected, with truth | Detected, without | False alarm, with truth | False alarm, without | Off the 5% target |
|---|---|---|---|---|---|---|---|---|---|
| U-T vs U-F2 | compression -0.05 | EWMA | F0 | 10 | 100% | 99% | 5.2% | 4.3% | no |
| U-T vs U-F2 | compression -0.05 | EWMA | F0 | 20 | 100% | 98% | 5.3% | 3.1% | yes |
| U-T vs U-F2 | compression -0.05 | EWMA | F1 | 10 | 100% | 99% | 5.2% | 4.3% | no |
| U-T vs U-F2 | compression -0.05 | EWMA | F1 | 20 | 100% | 98% | 5.3% | 3.1% | yes |
| U-T vs U-F2 | compression -0.05 | EWMA | F2 | 10 | 100% | 88% | 5.8% | 5.3% | no |
| U-T vs U-F2 | compression -0.05 | EWMA | F2 | 20 | 100% | 85% | 5.5% | 4.1% | no |
| U-T vs U-F2 | compression -0.05 | EWMA | F3 | 10 | 84% | 24% | 5.5% | 4.4% | no |
| U-T vs U-F2 | compression -0.05 | EWMA | F3 | 20 | 90% | 27% | 3.7% | 2.7% | yes |
| U-T vs U-F2 | compression -0.05 | EWMA | F4 | 10 | 100% | 99% | 6.9% | 4.6% | yes |
| U-T vs U-F2 | compression -0.05 | EWMA | F4 | 20 | 100% | 99% | 8.2% | 4.9% | yes |
| U-T vs U-F2 | compression -0.05 | EWMA | F5 | 10 | 84% | 57% | 7.2% | 6.2% | yes |
| U-T vs U-F2 | compression -0.05 | EWMA | F5 | 20 | 88% | 61% | 4.8% | 3.9% | no |

Detection in the base family (F0), CUSUM, 20 monitored runs, calibrated limits:

| Detector | No drift | Precision loss |
|---|---|---|
| residSD | 6% | 100% |
| B2 | 4% | 7% |
| B4 (profileResidSD) | 3% | 100% |
| U-T | 4% | 100% |
| U-F | 5% | 13% |
| U-F2 (B3 or B4) | 3% | 100% |

Confirmation on untouched seeds 910001 to 910400, precision loss:

| Family | Pair | Chart | Runs | With truth | Without | Difference, points (95% interval) | Verdict |
|---|---|---|---|---|---|---|---|
| F0 | U-T vs U-F | CUSUM | 10 | 100% | 13% | +87.0 (+83.5 to +90.2) | gain |
| F0 | U-T vs U-F | CUSUM | 20 | 100% | 16% | +84.2 (+80.8 to +87.5) | gain |
| F0 | U-T vs U-F | EWMA | 10 | 100% | 14% | +86.5 (+83.0 to +89.5) | gain |
| F0 | U-T vs U-F | EWMA | 20 | 100% | 18% | +81.5 (+77.5 to +85.3) | gain |
| F0 | residSD vs B2 | CUSUM | 10 | 100% | 6% | +94.0 (+91.5 to +96.2) | gain |
| F0 | residSD vs B2 | CUSUM | 20 | 100% | 6% | +93.8 (+91.2 to +96.0) | gain |
| F0 | residSD vs B2 | EWMA | 10 | 100% | 7% | +92.8 (+90.2 to +95.2) | gain |
| F0 | residSD vs B2 | EWMA | 20 | 100% | 8% | +91.5 (+88.8 to +94.2) | gain |
| F0 | residSD vs B4 | CUSUM | 10 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F0 | residSD vs B4 | CUSUM | 20 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F0 | residSD vs B4 | EWMA | 10 | 100% | 100% | +0.5 (+0.0 to +1.2) | no clear difference |
| F0 | residSD vs B4 | EWMA | 20 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F0 | U-T vs U-F2 | CUSUM | 10 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F0 | U-T vs U-F2 | CUSUM | 20 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F0 | U-T vs U-F2 | EWMA | 10 | 100% | 100% | +0.5 (+0.0 to +1.2) | no clear difference |
| F0 | U-T vs U-F2 | EWMA | 20 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F5 | U-T vs U-F | CUSUM | 10 | 100% | 5% | +94.8 (+92.5 to +96.8) | gain |
| F5 | U-T vs U-F | CUSUM | 20 | 100% | 4% | +95.8 (+93.8 to +97.5) | gain |
| F5 | U-T vs U-F | EWMA | 10 | 98% | 6% | +92.8 (+90.2 to +95.2) | gain |
| F5 | U-T vs U-F | EWMA | 20 | 99% | 6% | +93.5 (+91.0 to +95.8) | gain |
| F5 | residSD vs B2 | CUSUM | 10 | 100% | 8% | +92.0 (+89.5 to +94.5) | gain |
| F5 | residSD vs B2 | CUSUM | 20 | 100% | 6% | +93.8 (+91.2 to +96.0) | gain |
| F5 | residSD vs B2 | EWMA | 10 | 100% | 6% | +94.2 (+91.8 to +96.5) | gain |
| F5 | residSD vs B2 | EWMA | 20 | 100% | 7% | +93.0 (+90.5 to +95.2) | gain |
| F5 | residSD vs B4 | CUSUM | 10 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F5 | residSD vs B4 | CUSUM | 20 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F5 | residSD vs B4 | EWMA | 10 | 100% | 100% | +0.3 (+0.0 to +0.8) | no clear difference |
| F5 | residSD vs B4 | EWMA | 20 | 100% | 100% | +0.3 (+0.0 to +0.8) | no clear difference |
| F5 | U-T vs U-F2 | CUSUM | 10 | 100% | 100% | -0.3 (-0.8 to +0.0) | no clear difference |
| F5 | U-T vs U-F2 | CUSUM | 20 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F5 | U-T vs U-F2 | EWMA | 10 | 98% | 99% | -0.7 (-2.0 to +0.5) | no clear difference |
| F5 | U-T vs U-F2 | EWMA | 20 | 99% | 100% | -0.5 (-1.2 to +0.0) | no clear difference |
<!-- /PV2-A2 -->

**What this says.** Under the frozen rule, against the strongest truth-free charts we know how to build, the known amounts gave one gain: U-T against U-F2 under compression -0.05, on the EWMA chart, in all six families. The CUSUM chart does not agree, and that cell was not part of the confirmation on fresh seeds, so it is a gain under the rule that is not yet confirmed. Its truth side also had the higher false-alarm rate in every family and horizon. Where both sides caught nearly every series, the detection difference is about the size of that gap, so there the gain may be extra alarms; where the profile chart missed many series, the difference is far larger than the gap (addendum 2, cells of the claimed gain). No other pair meets the rule for any planted drift. The precision-loss gain in the first table came from comparing with B2; against `profileResidSD` (B4), which needs no amounts, there is no difference under the rule, and the confirmation on untouched seeds agrees (addendum 2, above). Rows marked "at ceiling" there are ones where both charts caught every series: they show no advantage at the tested drift size, not that the two charts are equal in power. Under compression, the truth-based slope and the profile slope both catch every planted drift within 20 runs in the base family, and no pair meets the rule at -0.10 or for the ramp. Under sensitivity loss the truth-free union wins: the spread of counts reacts to depth directly, while slope and residSD react less, though residSD still catches much of it; the truth-based `offset` catches that drift as well as the total count does.

**What the package is for, given that.** The known amounts give per-run figures in units a lab understands (a slope of 1, a detection limit in attomoles), and the real-data checks below test those figures. Apart from the one unconfirmed EWMA gain, they did not give a better alarm at the drift sizes tested. These per-run figures are model-dependent QC summaries: the residual spread also carries each molecule's own response, not only the run's precision. A lab that only wants alarms can chart the profile figures, which need no amounts, after checking their limits on its own runs. Where a truth-based figure loses a pair, it is paired with a truth-free chart that reacts to a different drift.

Detection in the base family (F0) at calibrated limits, share of series with any signal within 20 monitored runs (the null column is the false-alarm rate on fresh seeds; other families and the 10-run horizon are in `validation/out-v2/scores.csv`):

<!-- PV2-DETECT -->
CUSUM, 20 monitored runs:

| Detector | no drift | compression -0.05 | compression -0.10 | compression ramp | sensitivity loss | precision loss |
|---|---|---|---|---|---|---|
| slope | 6% | 100% | 100% | 100% | 13% | 25% |
| residSD | 6% | 4% | 7% | 4% | 75% | 100% |
| offset | 5% | 10% | 26% | 20% | 100% | 7% |
| U-T (slope or residSD) | 4% | 100% | 100% | 100% | 68% | 100% |
| B1 | 5% | 69% | 100% | 99% | 100% | 9% |
| B2 | 4% | 98% | 100% | 100% | 98% | 7% |
| B3 | 6% | 100% | 100% | 100% | 8% | 15% |
| U-F (B2 or B3) | 5% | 100% | 100% | 100% | 94% | 13% |

EWMA, 20 monitored runs:

| Detector | no drift | compression -0.05 | compression -0.10 | compression ramp | sensitivity loss | precision loss |
|---|---|---|---|---|---|---|
| slope | 5% | 100% | 100% | 100% | 12% | 25% |
| residSD | 6% | 5% | 5% | 4% | 62% | 100% |
| offset | 5% | 9% | 20% | 15% | 99% | 8% |
| U-T (slope or residSD) | 5% | 100% | 100% | 100% | 57% | 100% |
| B1 | 5% | 57% | 100% | 98% | 99% | 10% |
| B2 | 6% | 94% | 100% | 100% | 94% | 8% |
| B3 | 5% | 100% | 100% | 100% | 8% | 15% |
| U-F (B2 or B3) | 5% | 100% | 100% | 100% | 82% | 15% |
<!-- /PV2-DETECT -->

**Default limits.** Version 1 reported its charts at the textbook limits. With a 20-run baseline they alarm far more often than those limits suggest. False alarms with no drift at all, base family:

<!-- PV2-DEFAULTS -->
| Detector | CUSUM h = 5, 10 runs | CUSUM h = 5, 20 runs | EWMA L = 3, 10 runs | EWMA L = 3, 20 runs |
|---|---|---|---|---|
| slope | 6.4% | 16.2% | 9.0% | 16.8% |
| residSD | 8.2% | 18.7% | 11.1% | 19.2% |
| offset | 6.1% | 15.7% | 8.3% | 16.6% |
| U-T (slope or residSD) | 13.4% | 30.4% | 17.9% | 31.5% |
| B1 | 5.7% | 16.4% | 7.7% | 17.4% |
| B2 | 6.8% | 15.7% | 9.5% | 16.5% |
| B3 | 9.0% | 23.1% | 11.8% | 24.0% |
| U-F (B2 or B3) | 15.3% | 33.7% | 20.0% | 35.2% |
<!-- /PV2-DEFAULTS -->

**Setting limits from a short baseline.** The plan also tested `calibrateChart()`, a bootstrap that sets a chart's limit for a stated false-alarm rate from the baseline runs. It failed its own pass rule in 5 of the 6 families, so it is not in this release. The function, its tests and its full results are on the `calibrate-chart` branch. Until a version passes, expect the default limits to give more false alarms than they promise, as the table above shows.

**Correlated runs.** When runs drift slowly from day to day (family F5), every chart in this package alarms far too often. With no drift at all, 36.7% to 83.4% of series raised a false alarm within 20 monitored runs at the default limits (slope, residSD, offset and B3, both charts). The charts assume independent runs. If yours are not, model the correlation first.

**Confirmation on untouched seeds**, run once after everything above was written down (seeds 900001 to 900400, thresholds from the calibration seeds):

<!-- PV2-CONFIRM -->
| Family | Pair | Drift | Chart | Runs | With truth | Without | Difference, points (95% interval) | Verdict |
|---|---|---|---|---|---|---|---|---|
| F0 | U-T vs U-F | compression -0.10 | CUSUM | 10 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F0 | U-T vs U-F | compression -0.10 | CUSUM | 20 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F0 | U-T vs U-F | compression -0.10 | EWMA | 10 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F0 | U-T vs U-F | compression -0.10 | EWMA | 20 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F0 | slope vs B3 | compression -0.10 | CUSUM | 10 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F0 | slope vs B3 | compression -0.10 | CUSUM | 20 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F0 | slope vs B3 | compression -0.10 | EWMA | 10 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F0 | slope vs B3 | compression -0.10 | EWMA | 20 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F0 | U-T vs U-F | sensitivity loss | CUSUM | 10 | 53% | 78% | -25.5 (-31.5 to -19.2) | loss |
| F0 | U-T vs U-F | sensitivity loss | CUSUM | 20 | 70% | 94% | -24.5 (-29.5 to -19.5) | loss |
| F0 | U-T vs U-F | sensitivity loss | EWMA | 10 | 50% | 71% | -21.0 (-27.5 to -14.8) | loss |
| F0 | U-T vs U-F | sensitivity loss | EWMA | 20 | 54% | 83% | -28.7 (-34.8 to -23.0) | loss |
| F0 | slope vs B3 | sensitivity loss | CUSUM | 10 | 12% | 8% | +3.8 (+0.2 to +7.2) | gain |
| F0 | slope vs B3 | sensitivity loss | CUSUM | 20 | 14% | 11% | +2.8 (-1.2 to +6.8) | no clear difference |
| F0 | slope vs B3 | sensitivity loss | EWMA | 10 | 12% | 8% | +4.0 (+0.5 to +7.8) | gain |
| F0 | slope vs B3 | sensitivity loss | EWMA | 20 | 13% | 8% | +4.5 (+0.8 to +8.5) | gain |
| F5 | U-T vs U-F | compression -0.10 | CUSUM | 10 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F5 | U-T vs U-F | compression -0.10 | CUSUM | 20 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F5 | U-T vs U-F | compression -0.10 | EWMA | 10 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F5 | U-T vs U-F | compression -0.10 | EWMA | 20 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F5 | slope vs B3 | compression -0.10 | CUSUM | 10 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F5 | slope vs B3 | compression -0.10 | CUSUM | 20 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F5 | slope vs B3 | compression -0.10 | EWMA | 10 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F5 | slope vs B3 | compression -0.10 | EWMA | 20 | 100% | 100% | +0.0 (+0.0 to +0.0) | no clear difference (at ceiling) |
| F5 | U-T vs U-F | sensitivity loss | CUSUM | 10 | 27% | 43% | -16.5 (-22.5 to -11.0) | loss |
| F5 | U-T vs U-F | sensitivity loss | CUSUM | 20 | 36% | 50% | -14.8 (-20.8 to -8.8) | loss |
| F5 | U-T vs U-F | sensitivity loss | EWMA | 10 | 27% | 40% | -13.2 (-19.0 to -7.5) | loss |
| F5 | U-T vs U-F | sensitivity loss | EWMA | 20 | 30% | 39% | -8.5 (-14.0 to -2.8) | loss |
| F5 | slope vs B3 | sensitivity loss | CUSUM | 10 | 7% | 6% | +0.8 (-1.2 to +3.0) | no clear difference |
| F5 | slope vs B3 | sensitivity loss | CUSUM | 20 | 6% | 5% | +0.5 (-1.2 to +2.2) | no clear difference |
| F5 | slope vs B3 | sensitivity loss | EWMA | 10 | 6% | 8% | -1.5 (-3.8 to +0.8) | no clear difference |
| F5 | slope vs B3 | sensitivity loss | EWMA | 20 | 6% | 5% | +1.2 (-0.5 to +3.0) | no clear difference |
<!-- /PV2-CONFIRM -->

**Real data at the v2 code.** The SEQC and msqc1 checks of version 1 rerun with `minAmount = "auto"` over every run (`validation/out-v2/real_*.csv`). The msqc1 cross-peptide slope is reported as relative linearity only: each peptide has its own amount and its own response, so it is not a test of absolute linearity. No public, licensed series of many runs of one spike-in mix with run order was found for this version, so **the control charts are tested on simulation only**; the real data test the per-run figures, not the charts.

<!-- PV2-REAL -->
msqc1 dilution series at the v2 code (`minAmount = "auto"` over every run):

| Instrument | Runs | Doses | Slope against dose | 95% interval |
|---|---|---|---|---|
| QExactive | 18 | 6 | 0.988 | 0.929 to 1.048 |
| QExactiveHF | 18 | 6 | 0.975 | 0.950 to 1.000 |
| QTRAP | 18 | 6 | 0.766 | 0.711 to 0.822 |
| TRIPLETOF | 18 | 6 | 0.982 | 0.936 to 1.028 |
| TSQVantage | 18 | 6 | 0.882 | 0.830 to 0.933 |
<!-- /PV2-REAL -->

## Validation, version 1 (commit 061caed)

The plan, the detectors and the rule for claiming a gain were committed before any result was computed (`validation/PLAN.md`). Every figure below is filled in from `validation/out/` by `validation/render_readme.py`. Version 2 replaces its scoring: v1 compared detectors at the textbook limits, scored each detector only on series where it had not alarmed early, and had weaker truth-free baselines. The v1 gain over B1 claimed below is superseded by the version 2 results above.

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

## Data and licences

The validation reads its data at run time from Bioconductor packages; no data from them is copied into this repository.

| Package | Version read | Licence (Bioconductor index) | Used for |
|---|---|---|---|
| erccdashboard | 1.46.0 | GPL (>= 2) | SEQC lab 5 ERCC counts (V2); the ERCC Mix 1 and Mix 2 concentrations |
| msqc1 | 1.40.0 | GPL | heavy peptide areas, 8 replicates and dilution series (V3) |

One table does ship with the package: the 92 ERCC concentrations in `inst/extdata/ercc_mixes.tsv`, taken from erccdashboard. erccdashboard is licensed GPL (>= 2) and this package GPL-3, which that licence allows; `inst/script/make_ercc_mixes.R` records the source and the checks.

## Checks

`R CMD check`, BiocCheck and the testthat suite run on every push (`.github/workflows/check.yaml`). Tests use small fixtures whose answers can be worked by hand: an exact power law, a hand-computed ratio table, and drift series whose chart values follow from the formulas.

## Licence

GPL-3. See "Data and licences" for the ERCC table.
