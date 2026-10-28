# spikeWatch validation plan

Written and committed before any validation result was computed. Results that break these rules are reported as they come, not re-scored.

## V1. Planted drift in simulated series

Design: `erccDesign()`, Mix1. Each series has 40 runs; runs 1 to 20 are the baseline; drift starts at run 31 as a step. Package defaults otherwise (`scale = 10`, `slope = 1`, `noise = 0.25`, `depthSD = 0.5`).

Scenarios, 200 series each, seeds 1 to 200:

| Scenario | driftType | driftSize |
|---|---|---|
| null | none | 0 |
| compression small | slope | -0.05 |
| compression | slope | -0.10 |
| sensitivity loss | offset | -1 (half the signal) |
| precision loss | noise | +0.25 |

Detectors, all two-sided tabular CUSUM with `k = 0.5`, `h = 5`, sigma from the baseline moving range:

- Truth-based (spikeWatch, `minAmount = 1`): `slope`, `residSD`, `offset`, and "any of slope or residSD" (union).
- Without truth (what a lab can chart from spike counts alone, knowing which rows are spikes but not their amounts): log2 total spike count per run (B1), and the standard deviation of log2 counts across detected spikes (B2), which also shrinks under compression.

Scores per detector and scenario:

- False alarm: any signal in runs 21 to 30 (before the drift). In the null scenario, also any signal in runs 21 to 40.
- Detection: a first signal at run 31 or later, among series with no false alarm. Reported as the share detected by run 40 and the median delay (first signal minus 31).

No claim of a gain is made unless the truth-based detector both detects more often and has a false alarm rate no higher than the baseline it is compared with, in the same scenario.

## V2. Real RNA-seq: SEQC lab 5 ERCC counts

Data: `UHRR.HBRR.countDat` from erccdashboard (GPL (>= 2)); UHRR runs carry ERCC Mix 1 and HBRR runs Mix 2 (erccdashboard vignette, lines 271 to 272 of erccdashboard.Rmd, version 1.46.0). 4 runs each.

- Per run: `runMetrics` on counts per kilobase (counts divided by ERCC length from `ERCCDef`, times 1000), `minAmount = 1`.
- Between groups: `ratioAccuracy(UHRR, HBRR, center = "spikes")` on raw counts (length cancels in a ratio).
- Known answer: expected log2 ratios 2, 0, -0.585 and -1. Report median error and RMSE per subgroup as they come.

## V3. Real targeted proteomics: msqc1

Data: msqc1 (GPL). Heavy peptide areas summed over transitions per run.

- `msqc1_8rep`: 14 heavy peptides at known amounts per vial (`peptides$SIL.peptide.per.vial`, unit not stated in the package [CONFIRM]), 8 runs per instrument, 5 instruments. Per run: slope, R squared, residual SD, detection limit. Known answer: a linear response gives slope 1. Report per instrument as they come.
- `msqc1_dil`: the same peptides at a known per-sample dose (`relative.amount`). Truth by `addTruth(..., dose = "relative.amount")`. Known answer: per-run `offset` regressed on log2 dose has slope 1 within each instrument. Report the fitted slope per instrument.
