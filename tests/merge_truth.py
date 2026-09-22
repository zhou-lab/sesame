#!/usr/bin/env python3
"""Merge the C and R sides into one frozen truth file, row by row.

Every row records where its value came from: `R` where the two agree, or
`C:<code>` where they do not -- and a `C:` row is only allowed when <code> is
a reason already written down in DIVERGENCES.md. That is the whole discipline
here: the file may record a disagreement, but never an unexplained one.

    merge_truth.py <test> <c.tsv> <r.tsv> <ordering|""> <out.tsv>
"""
import gzip, os, re, sys

test, cf, rf, ordf, out = sys.argv[1:6]
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def load(fn):
    d = {}
    if not (fn and os.path.exists(fn) and os.path.getsize(fn)):
        return d
    with open(fn) as f:
        head = next(f).rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            d[p[0]] = p[1:]
    return d, head

def reasons():
    """The codes DIVERGENCES.md documents, so an undocumented one cannot ship."""
    p = os.path.join(ROOT, "DIVERGENCES.md")
    if not os.path.exists(p):
        return set()
    return set(re.findall(r"^###\s+`([A-Za-z0-9:_-]+)`", open(p).read(), re.M))

KNOWN = reasons()
C, chead = load(cf) if load(cf) else ({}, [])
R, rhead = load(rf) if load(rf) else ({}, [])

def emit(rows, cols):
    with open(out, "w") as f:
        f.write("\t".join(["Probe_ID"] + cols + ["src"]) + "\n")
        for pid, vals, src in rows:
            f.write("\t".join([pid] + list(vals) + [src]) + "\n")
    n_c = sum(1 for _, _, s in rows if s.startswith("C:"))
    print(f"  {os.path.basename(out)}: {len(rows)} rows, {n_c} C-sourced")

def check(code):
    if code not in KNOWN:
        sys.exit(f"merge_truth: reason `{code}` is not documented in DIVERGENCES.md.\n"
                 f"  A frozen value nobody can explain is a bug with a test around it.\n"
                 f"  Add a `### \\`{code}\\`` section saying what differs and why, then re-run.")
    return code

if test == "vcf":
    ## the channel each side used decides whether a disagreement is annotation
    ocol = {}
    with gzip.open(ordf, "rt") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t"); ocol[p[0]] = p[3]
    ri = {k: i for i, k in enumerate(rhead[1:])}
    rows, unexplained = [], []
    for pid in sorted(set(C) | set(R)):
        c, r = C.get(pid), R.get(pid)
        if c and r:
            rcol = r[ri["Rcol"]] if "Rcol" in ri and len(r) > ri["Rcol"] else ""
            agree_gt = c[0] == r[0]
            agree_pvf = abs(float(c[2]) - float(r[2])) < 1e-6
            if agree_gt and agree_pvf:
                rows.append((pid, [r[0], r[1], r[2]], "R"))
            elif rcol and ocol.get(pid) and rcol != ocol[pid]:
                rows.append((pid, [c[0], c[1], c[2]], check("C:vcf-channel")))
            else:
                unexplained.append(pid)
        elif c:
            rows.append((pid, c[:3], check("C:vcf-extra-probe")))
    if unexplained:
        sys.exit(f"merge_truth: {len(unexplained)} probes differ for no recorded reason "
                 f"(first: {unexplained[:3]}). Explain them before freezing.")
    emit(rows, ["GT", "GS", "PVF"])

elif test == "pneg":
    rows = []
    if not R:                      ## R could not compute it at all
        code = check("C:pneg-no-negctl")
        rows = [(pid, v[:1], code) for pid, v in sorted(C.items())]
    else:
        unexplained = []
        def num(x):
            ## R writes NA for a probe it reports but cannot score; C writes NA
            ## too. NA==NA is agreement, NA against a number is not.
            return None if x in ("NA", "NaN", "") else float(x)
        for pid in sorted(set(C) | set(R)):
            c, r = C.get(pid), R.get(pid)
            if c and r:
                cv, rv = num(c[0]), num(r[0])
                if cv is None and rv is None:
                    rows.append((pid, ["NA"], "R"))
                elif cv is not None and rv is not None and abs(cv - rv) <= 1e-9:
                    rows.append((pid, [r[0]], "R"))
                else:
                    unexplained.append(pid)
            elif c:                ## positional ordering carries control probes R omits
                rows.append((pid, c[:1], check("C:pneg-control-probe")))
        if unexplained:
            sys.exit(f"merge_truth: {len(unexplained)} pneg values differ for no recorded "
                     f"reason (first: {unexplained[:3]}).")
    emit(rows, ["pval"])
else:
    sys.exit(f"merge_truth: unknown test {test}")
