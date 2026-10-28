"""Plant each bug in mutants.tsv into a fresh copy of the committed head
and run the test suite on it. A mutant is killed when at least one test
fails or errors. Run from the repository root; writes
validation/mutants/results.tsv.

Each line of mutants.tsv: file, exact text to replace, replacement. The
text must occur exactly once, or the mutant is reported as not applied.
"""
import os
import subprocess
import sys
import tempfile

root = os.getcwd()
head = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"],
                               text=True).strip()
lines = open("validation/mutants/mutants.tsv").read().strip().split("\n")
rows = []
with tempfile.TemporaryDirectory() as tmp:
    for i, line in enumerate(lines, 1):
        f, old, new = line.split("\t")
        d = os.path.join(tmp, "m%02d" % i)
        os.mkdir(d)
        subprocess.run("git archive HEAD | tar -x -C '%s'" % d, shell=True,
                       check=True)
        p = os.path.join(d, f)
        s = open(p).read()
        if s.count(old) != 1:
            rows.append(("m%02d" % i, f, "NOT_APPLIED", ""))
            continue
        open(p, "w").write(s.replace(old, new))
        code = ("r <- as.data.frame(testthat::test_local('%s', "
                "reporter = 'silent', stop_on_failure = FALSE)); "
                "cat('FAILED', sum(r$failed) + sum(r$error), '\\n')" % d)
        r = subprocess.run(["Rscript", "-e", code], capture_output=True,
                           text=True, timeout=600)
        k = [x for x in (r.stdout + r.stderr).split("\n")
             if x.startswith("FAILED")]
        n = int(k[0].split()[1]) if k else -1
        verdict = "KILLED" if n > 0 or r.returncode != 0 else "SURVIVED"
        rows.append(("m%02d" % i, f, verdict, "failing tests: %d" % n))
with open("validation/mutants/results.tsv", "w") as out:
    out.write("# head %s\n" % head)
    for r in rows:
        out.write("\t".join(r) + "\n")
print(open("validation/mutants/results.tsv").read())
sys.exit(0 if all(r[2] == "KILLED" for r in rows) else 1)
