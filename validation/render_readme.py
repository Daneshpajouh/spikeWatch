"""Fill the validation tables in README.md from validation/out/*.csv.

Run from the repository root after the three validation scripts. Every
figure in the README's Validation section comes from here; none is typed.
"""
import csv
import re
import sys

OUT = "validation/out/"


def rows(name):
    with open(OUT + name, newline="") as f:
        return list(csv.DictReader(f))


def pct(x):
    return "%.1f%%" % (100 * float(x))


def num(x, d=2):
    return ("%." + str(d) + "f") % float(x)


def table(head, body):
    lines = ["| " + " | ".join(head) + " |",
             "|" + "---|" * len(head)]
    lines += ["| " + " | ".join(r) + " |" for r in body]
    return "\n".join(lines)


def v0():
    body = [[num(r["minAmount"], 0), num(r["meanSlope"], 3),
             num(r["sdSlope"], 4), num(r["corWithDepth"], 3)]
            for r in rows("v0_min_amount.csv")]
    return table(["minAmount", "Mean slope", "SD of slope",
                  "Correlation with depth"], body)


def v1():
    order = ["null", "compression_small", "compression", "sensitivity_loss",
             "precision_loss"]
    dets = ["slope", "residSD", "offset", "B1_total", "B2_spread"]
    s = {(r["scenario"], r["detector"]): r for r in rows("v1_summary.csv")}
    body = []
    for sc in order:
        for d in dets:
            r = s[(sc, d)]
            delay = r["medianDelay"]
            body.append([sc.replace("_", " "), d, pct(r["falseAlarm21to30"]),
                         pct(r["detectedBy40"]) if sc != "null" else "",
                         num(delay, 1) if sc != "null" and delay != "NA"
                         else ""])
    return table(["Scenario", "Chart", "False alarm, runs 21 to 30",
                  "Detected by run 40", "Median delay (runs)"], body)


def v2():
    m = rows("v2_run_metrics.csv")
    body = [[r["run"], r["nDetected"] + " of " + r["nSpikes"],
             num(r["slope"], 3), num(r["r2"], 3), num(r["lod"], 3)]
            for r in m]
    t1 = table(["Run", "Detected", "Slope", "R squared",
                "Detection limit (attomol/uL)"], body)
    s = rows("v2_ratio_summary.csv")
    body = [[num(r["expectedLog2Ratio"], 3), r["nSpikes"],
             num(r["medianError"], 3), num(r["rmse"], 3)] for r in s]
    t2 = table(["Expected log2 ratio", "Spikes", "Median error (log2)",
                "RMSE (log2)"], body)
    return t1 + "\n\n" + t2


def v3():
    a = rows("v3_8rep_by_instrument.csv")
    body = [[r["instrument"], r["runs"], num(r["slopeMedian"], 3),
             num(r["slopeMin"], 3) + " to " + num(r["slopeMax"], 3),
             num(r["r2Median"], 3)] for r in a]
    t1 = table(["Instrument", "Runs", "Slope (median)", "Slope range",
                "R squared (median)"], body)
    b = rows("v3_dil_linearity.csv")
    body = [[r["instrument"], r["runs"], r["doses"], num(r["doseSlope"], 3),
             num(r["ciLow"], 3) + " to " + num(r["ciHigh"], 3)] for r in b]
    t2 = table(["Instrument", "Runs", "Doses", "Slope against dose",
                "95% interval"], body)
    return t1 + "\n\n" + t2


def main():
    path = "README.md"
    text = open(path).read()
    for key, fn in (("V0", v0), ("V1", v1), ("V2", v2), ("V3", v3)):
        pat = re.compile(r"(<!-- %s -->\n).*?(<!-- /%s -->)" % (key, key),
                         re.S)
        if not pat.search(text):
            sys.exit("README has no %s block" % key)
        text = pat.sub(lambda m: m.group(1) + fn() + "\n" + m.group(2),
                       text)
    open(path, "w").write(text)


if __name__ == "__main__":
    main()
