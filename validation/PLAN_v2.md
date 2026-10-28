# spikeWatch validation plan v2

Written and committed before any v2 result was computed. It replaces the scoring rules of PLAN.md (v1) for every number reported against the v2 code. v1 results stay in the README, labelled as v1 and as computed at commit 061caed. Results that break these rules are reported as they come, not re-scored. Any change to this plan after the first v2 number is computed goes in a dated addendum below, says what changed and why, and leaves the numbers already computed under the old rule in place.

## What changed from v1, and why

- v1 compared detectors at fixed limits, so a detector with more false alarms could look more sensitive. v2 calibrates every detector to the same in-control false-alarm probability before comparing detection.
- v1 dropped series with a false alarm from the detection score, separately for each detector, so the detectors were scored on different series. v2 scores every series for every detector.
- v1 had two weak baselines that do not use the known amounts. v2 adds a third, the profile slope (B3), which is the strongest truth-free chart we know how to build.
- v1 simulated from the same model the estimator assumes. v2 adds five departures from that model and reports each apart.
- v1 published the union detector's false alarm rate only in a log. v2 publishes every detector at both horizons, at default and at calibrated limits.

## Series

Design `erccDesign()`, Mix1. Each series has 40 runs. Runs 1 to 20 are the baseline; runs 21 to 40 are monitored. Drift, when planted, starts at run 21 (the first monitored run), so every chart starts from its fresh state (zero-state). Drift that starts after a stretch of in-control monitoring (steady-state) is not studied in v2, and the README says so.

## Simulation families

Each family is reported apart. Parameters not named take the `simulateSpikeRuns()` defaults (`scale = 10`, `slope = 1`, `noise = 0.25`, `depthSD = 0.5`, Poisson counts, pivot at the median amount).

| Family | Departure from the base model |
|---|---|
| F0 base | none (the v1 model) |
| F1 top pivot | the slope turns about the highest design amount, so compression lowers every spike below the top |
| F2 overdispersed | negative binomial counts with size 10 in place of Poisson |
| F3 saturation | the log2 mean is compressed by half above the log2 mean of the 90th percentile design amount at depth 0 |
| F4 spike biases | each spike carries a fixed log2 bias, normal with SD 0.5, drawn once per series and kept in every run |
| F5 correlated runs | run depth and a run-level slope wobble (SD 0.02) both follow AR(1) with correlation 0.6, at the stationary SDs given |

The simulator gains arguments for these departures. Their names and defaults are fixed in code before any v2 number is computed; the defaults reproduce the v1 model exactly (same draws for the same seed).

## Drift scenarios

| Scenario | Type | Size | Shape |
|---|---|---|---|
| null | none | 0 | |
| compression small | slope | -0.05 | step |
| compression | slope | -0.10 | step |
| compression ramp | slope | -0.10 | ramp over 10 runs |
| sensitivity loss | offset | -1 | step |
| precision loss | noise | +0.25 | step |

## Figures charted

Truth-based (spikeWatch): `slope`, `residSD` and `offset` from `runMetrics(minAmount = "auto", baseline = 20)`, so the fit range is chosen from the baseline runs only.

Truth-free (what a lab can chart knowing which rows are spikes but not their amounts):

- B1: log2 total spike count per run.
- B2: SD of log2 counts across spikes detected in the run.
- B3: profile slope from `profileMetrics()`: the least-squares slope of log2 (count + 0.5) on the baseline-mean log2 (count + 0.5) of the same spikes, over spikes detected in every baseline run. The spike set and the reference profile come from runs 1 to 20 only.

Unions: truth-based U-T (slope or residSD); truth-free U-F (B2 or B3).

Charts: two-sided tabular CUSUM with `k = 0.5`, and EWMA with `lambda = 0.2`, each with centre and sigma from the baseline as in `driftChart()` (moving range). A detector is one figure (or union) on one chart.

## Horizons and target

Horizon H = 10 (runs 21 to 30) and H = 20 (runs 21 to 40). Target in-control false-alarm probability alpha = 0.05 within the horizon. Each horizon has its own threshold.

## Calibration (equal false alarms before any comparison)

For each family, detector and horizon, the threshold is set on 1,000 null calibration series (seeds 100001 to 101000) as the 950th smallest of the per-series maxima of the chart statistic over the horizon: for CUSUM the larger of the two sums in sigma units; for EWMA the distance from the centre divided by the exact time-varying EWMA standard deviation. A union uses one threshold for both of its members, set so that the union's maxima meet the same rule.

These seeds are used for calibration only and never for scoring.

## Scoring

Evaluation seeds: null, 1,000 series per family (seeds 200001 to 201000); each drift scenario, 400 series per family (seeds 300001 to 300400). Every detector sees the same series.

- Achieved false-alarm probability: share of evaluation null series with any signal within the horizon, with a Clopper-Pearson 95% interval.
- Detection probability: share of drift series with any signal within the horizon. No series is dropped for any reason.
- Delay: runs from drift start to first signal (first signal run minus 20), censored at H + 1 when there is none; reported as the mean of the censored values and as the median among detected series.

## Comparisons and the gain rule

Pre-specified pairs, for every family, drift scenario, chart and horizon:

- C1: U-T against U-F.
- C2: slope against B3.
- C3: residSD against B2.
- C4: offset against B1.

Each difference in detection probability gets a paired bootstrap 95% percentile interval: 2,000 resamples of series, the same resample for both detectors.

- A gain is reported for a cell only when the interval's lower end is above 0; a loss only when its upper end is below 0; otherwise "no clear difference".
- With hundreds of cells, some intervals will exclude 0 by chance. A headline claim about a pair, drift and chart is made only when it holds in at least five of the six families, at both horizons. Every cell is published.
- Both detectors in a pair are calibrated to the same alpha, and their achieved evaluation false-alarm rates are printed beside the difference. If either achieved rate's interval excludes 0.05, the cell is marked.

## Default limits, published in full

For every detector, chart, family and horizon: the achieved false-alarm probability and detection probability at the package defaults (CUSUM `h = 5`; EWMA `L = 3`), on the evaluation seeds. This includes the v1 union detector, so its false-alarm rate is in the README beside its detection rate.

## calibrateChart() check

`calibrateChart()` sets a chart threshold from one series' own baseline, by bootstrap: each replicate resamples the 20 baseline values, re-estimates centre and sigma from the resample, resamples H future in-control values from the baseline, and records the chart's maximum statistic over the horizon under the re-estimated centre and sigma. The threshold is the (1 - alpha) quantile of those maxima, so estimation error in the baseline is inside the limit. This calibrates the false-alarm probability averaged over baselines (unconditional). It does not give the conditional guarantee of Gandy and Kvaløy (a bound that holds with stated probability for the baseline actually drawn); the README says so and cites that method.

Check, for each family, chart and horizon, on the evaluation null series: `calibrateChart(alpha = 0.05, B = 500)` on each series' own runs 1 to 20, for slope, residSD, offset and B3; then the achieved false-alarm probability across series with its Clopper-Pearson interval, against 0.05 and against the default limits. Detection probability with these limits on the drift series is reported too. Bootstrap seeds are fixed per series (the series seed plus 5,000,000).

Pass rule, stated now: the method "holds its target" in a family when the upper end of the achieved interval is at or below 0.08 for every figure, chart and horizon in that family. Families where it fails are named in the README.

## Confirmation on untouched seeds

After every number above is computed and written down, one confirmation run repeats the C1 and C2 comparisons for families F0 and F5, compression and sensitivity loss, both charts and both horizons, on seeds 900001 to 900400 with calibration reused from the calibration seeds. It runs once. Its results are reported as they come, beside the main ones, whether they agree or not.

## Real data

- V2 (SEQC lab 5 ERCC counts, erccdashboard) and V3 (msqc1) are rerun at the v2 code with the v1 rules, except that `minAmount` is `"auto"` with every run as the baseline, and the chosen amount is printed.
- V3 `msqc1_8rep`: each peptide has one amount per vial, so the cross-peptide slope says how response varies between peptides at different amounts. It is reported as relative linearity, never as absolute linearity, and never for a design where all spikes share one amount. Absolute accuracy is the `offset`, labelled with its unit [CONFIRM] as in v1.
- V3 `msqc1_dil`: the offset is a depth-like figure; the dose-adjusted check (offset regressed on log2 dose, expected slope 1) is reported as in v1.

## A real longitudinal series

Before the README is finished, a search for a public, licensed series of many runs of one spike-in mix with run order or dates (ERCC, sequins, SIRV, UPS or the MSQC1 kit) is logged with what was checked. If one is found, a dated addendum to this plan names it and its analysis before any of its numbers are computed. If none is found, the README says plainly that the control charts are tested on simulation only, and that the real data sets test the per-run figures, not the charts.

## Addenda

None yet.
