# Plan v2, addendum 2: a truth-free residual comparator

Written at 2026-10-28 00:29 PT, before any series or score below was computed. The
rules of PLAN_v2.md and its first addendum stay as they are, and so do the
results already reported under them.

## Why

A hostile review of the v2 results (round 2, two reviewers, both HIGH)
found that the one headline gain, precision loss, was scored against B2,
the spread of log2 counts across spikes. `profileMetrics()` already
returns `profileResidSD`, the residual spread of each run's fit on the
baseline profile: the truth-free twin of `residSD`. It was left out, so
the gain was never tested against its natural competitor.

## What is added

- Figure B4: `profileResidSD` from `profileMetrics(se, baseline = 20)`,
  the same call that gives B3.
- Detectors: B4 alone, and U-F2, the union of B3 and B4, calibrated as
  U-F is (one threshold for the union's maxima).
- Pairs: C5, `residSD` against B4; C6, U-T against U-F2.
- Everything else as in PLAN_v2.md: the same seeds, families, drifts,
  charts, horizons, calibration on seeds 100001 to 101000, paired
  bootstrap, and the gain rule (five of six families at both horizons).

## Checks before scoring

The series are simulated again with the same seeds and code, with B4 as a
seventh column. The six existing figures must equal the saved v2 series
exactly, in every family; if any differs, scoring stops and the README
says so. C1 to C4 must reproduce the reported comparisons exactly.

## Confirmation

After the C5 and C6 results are written down, one run on untouched seeds
910001 to 910400, families F0 and F5, precision loss only, pairs C1, C3,
C5 and C6, both charts and both horizons, with thresholds from the
calibration seeds. It runs once and is reported as it comes.

## Reporting

C5 and C6 are reported beside C1 to C4. If either shows no gain under the
rule, the README says the precision-loss gain does not hold against the
strongest truth-free chart we know how to build. F1 repeats F0 for this
drift; the count with F1 left out is reported too.

`calibrateChart()` is not in this release, so its block in
`pv2_2_score.R` is skipped in this run.
