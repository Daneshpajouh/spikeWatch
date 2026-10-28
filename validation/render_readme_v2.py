"""Fill the v2 validation tables in README.md from validation/out-v2/*.csv.

Run from the repository root after the v2 scripts. Every figure in the
README's v2 section comes from here; none is typed.
"""
import csv
import re

OUT = "validation/out-v2/"
FAMILIES = ["F0", "F1", "F2", "F3", "F4", "F5"]
DRIFTS = ["compression_small", "compression", "compression_ramp",
          "sensitivity_loss", "precision_loss"]
DRIFT_NAMES = {"null": "no drift", "compression_small": "compression -0.05",
               "compression": "compression -0.10",
               "compression_ramp": "compression ramp",
               "sensitivity_loss": "sensitivity loss",
               "precision_loss": "precision loss"}
PAIRS = {"C1": "U-T vs U-F", "C2": "slope vs B3", "C3": "residSD vs B2",
         "C4": "offset vs B1"}
DETS = ["slope", "residSD", "offset", "UT", "B1", "B2", "B3", "UF"]
DET_NAMES = {"UT": "U-T (slope or residSD)", "UF": "U-F (B2 or B3)"}


def rows(name):
    with open(OUT + name, newline="") as f:
        return list(csv.DictReader(f))


def pct(x):
    return "%.0f%%" % (100 * float(x))


def pct1(x):
    return "%.1f%%" % (100 * float(x))


def table(head, body):
    lines = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    lines += ["| " + " | ".join(r) + " |" for r in body]
    return "\n".join(lines)


def headline():
    hl = rows("headline.csv")
    get = {(r["pair"], r["scenario"], r["chart"]): r for r in hl}
    body = []
    for p in PAIRS:
        for d in DRIFTS:
            cells = []
            for ch in ("cusum", "ewma"):
                r = get[(p, d, ch)]
                g, lo = int(r["families_gain_both_H"]), int(r["families_loss_both_H"])
                word = {"gain": "**gain**", "loss": "**loss**", "none": "none"}[r["headline"]]
                cells.append("%s (gain in %d, loss in %d)" % (word, g, lo))
            body.append([PAIRS[p], DRIFT_NAMES[d]] + cells)
    return table(["Pair (with truth vs without)", "Drift", "CUSUM", "EWMA"], body)


def dedup():
    """Families with a gain or loss at both horizons, F1 left out where it
    repeats F0 (every drift but compression, which F1's pivot changes)."""
    c = rows("comparisons.csv")
    same = {}
    for sc in ("sensitivity_loss", "precision_loss"):
        f0 = {(r["pair"], r["chart"], r["H"]): r for r in c if r["family"] == "F0" and r["scenario"] == sc}
        f1 = {(r["pair"], r["chart"], r["H"]): r for r in c if r["family"] == "F1" and r["scenario"] == sc}
        same[sc] = sum(f0[k]["detect_truth"] == f1[k]["detect_truth"] and
                       f0[k]["detect_without"] == f1[k]["detect_without"] for k in f0)
    hl = {(r["pair"], r["scenario"], r["chart"]): r["headline"] for r in rows("headline.csv")}
    body = []
    for p in PAIRS:
        for sc in ("sensitivity_loss", "precision_loss"):
            for ch in ("cusum", "ewma"):
                if hl[(p, sc, ch)] == "none":
                    continue
                verdict = hl[(p, sc, ch)]
                fams = sorted({r["family"] for r in c if r["family"] != "F1"})
                n = 0
                for f in fams:
                    v = [r["verdict"] for r in c if r["family"] == f and r["pair"] == p
                         and r["scenario"] == sc and r["chart"] == ch]
                    n += len(v) == 2 and all(x == verdict for x in v)
                body.append([PAIRS[p], DRIFT_NAMES[sc], ch.upper(), verdict,
                             "%d of %d" % (n, len(fams))])
    note = ("F1 had the same detection rates as F0 in %d of 16 sensitivity-loss cells and %d of 16 "
            "precision-loss cells." % (same["sensitivity_loss"], same["precision_loss"]))
    return note + "\n\n" + table(["Pair", "Drift", "Chart", "Claim", "Families with F1 left out"], body)


def detection(chart, H):
    s = rows("scores.csv")
    get = {(r["detector"], r["scenario"]): r for r in s
           if r["family"] == "F0" and r["chart"] == chart and r["H"] == str(H)
           and r["limits"] == "calibrated"}
    body = []
    for d in DETS:
        body.append([DET_NAMES.get(d, d)] +
                    [pct(get[(d, sc)]["rate"]) for sc in ["null"] + DRIFTS])
    return table(["Detector"] + [DRIFT_NAMES[x] for x in ["null"] + DRIFTS], body)


def defaults():
    s = rows("scores.csv")
    get = {(r["detector"], r["chart"], r["H"]): r for r in s
           if r["family"] == "F0" and r["limits"] == "default"
           and r["scenario"] == "null"}
    body = []
    for d in DETS:
        body.append([DET_NAMES.get(d, d)] +
                    [pct1(get[(d, ch, H)]["rate"])
                     for ch in ("cusum", "ewma") for H in ("10", "20")])
    return table(["Detector", "CUSUM h = 5, 10 runs", "CUSUM h = 5, 20 runs",
                  "EWMA L = 3, 10 runs", "EWMA L = 3, 20 runs"], body)


def confirm():
    try:
        cf = rows("confirm.csv")
    except FileNotFoundError:
        return "Not yet run."
    body = []
    for r in cf:
        body.append([r["family"], PAIRS[r["pair"]], DRIFT_NAMES[r["scenario"]],
                     r["chart"].upper(), r["H"], pct(r["detect_truth"]),
                     pct(r["detect_without"]),
                     "%+.1f (%+.1f to %+.1f)" % tuple(100 * float(r[k]) for k in
                                                       ("diff", "ci_low", "ci_high")),
                     r["verdict"]])
    return table(["Family", "Pair", "Drift", "Chart", "Runs", "With truth", "Without",
                  "Difference, points (95% interval)", "Verdict"], body)


def fill(text, tag, body):
    pat = re.compile(r"(<!-- %s -->\n).*?(<!-- /%s -->)" % (tag, tag), re.S)
    if not pat.search(text):
        raise SystemExit("README has no %s block" % tag)
    return pat.sub(lambda m: m.group(1) + body + "\n" + m.group(2), text)


readme = open("README.md").read()
readme = fill(readme, "PV2-HEADLINE", headline())
readme = fill(readme, "PV2-DEDUP", dedup())
readme = fill(readme, "PV2-DETECT", "CUSUM, 20 monitored runs:\n\n" + detection("cusum", 20)
              + "\n\nEWMA, 20 monitored runs:\n\n" + detection("ewma", 20))
readme = fill(readme, "PV2-DEFAULTS", defaults())
readme = fill(readme, "PV2-CONFIRM", confirm())
open("README.md", "w").write(readme)
print("README v2 tables filled")
