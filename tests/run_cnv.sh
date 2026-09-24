#!/bin/sh
# cnv: validate the copy-number numeric core (target~normals OLS + log2 ratio)
# against R's lm on IDENTICAL inputs, and sanity-check the genome binning. Feeds R
# the same per-probe totals C uses (extracted from C's own .cg via describe-probe),
# so this is a pure numeric check -- no data-lineage or channel ambiguity. Needs a
# store with EPICv2 ordering+coord+cnvnormals, the hg38 genome, an EPICv2 IDAT, and
# Rscript; skips cleanly otherwise.
set -eu

## The R oracle. Overridable because the binary is not called the same
## thing everywhere: on the lab HPC plain `Rscript` is a different R
## without a usable sesame, and the one to use is Rscript-4.6.0.
RSCRIPT=${RSCRIPT:-Rscript}

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
## The shared store yame fetch fills. It is keyed on the browser
## hierarchy -- species/platform -- so a platform's assets sit directly
## at $store/$plat/ and a genome build's at $store/<build>/.
yhome=${YAME_DATA_HOME:-${XDG_DATA_HOME:-$HOME/.local/share}/yame}
store=$yhome
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "SKIP cnv: no Rscript"; exit 0; }

plat=EPICv2
ord="$store/$plat/$plat.ordering.tsv.gz"
crd="$store/$plat/$plat.hg38.coord.tsv.gz"
nrm="$store/$plat/$plat.cnvnormals.cg"
seq="$yhome/hg38/seqinfo.tsv.gz"
gap="$yhome/hg38/gaps.tsv.gz"
pfx="$idats/$plat/206909630042_R08C01"
for f in "$ord" "$crd" "$nrm" "$seq" "$gap"; do
    [ -f "$f" ] || { echo "SKIP cnv: missing $f (fetch $plat + genome hg38, make cnv-normals)"; exit 0; }
done
if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
    echo "SKIP cnv: no IDAT $pfx"; exit 0; fi

export YAME_DATA_HOME="$yhome"

# target total intensity (raw), then C's CNV probe signals
"$bin" preprocess --prep "" --raw-signal --output total_intensity --out "$work" "$pfx" 2>/dev/null
"$bin" cnv --probes --platform $plat --genome hg38 "$work/total_intensity.cg" "$work/c_seg.tsv" "$work/c_probes.tsv" 2>/dev/null
"$bin" cnv          --platform $plat --genome hg38 "$work/total_intensity.cg" "$work/c_seg.tsv" "$work/c_bins.tsv"   2>/dev/null

# extract C's own totals/normals/coords (Probe_ID-keyed) for the R oracle
"$bin" describe-probe --index "$ord" "$work/total_intensity.cg" 2>/dev/null > "$work/tgt.tsv"
"$bin" describe-probe --all --index "$ord" "$nrm" 2>/dev/null > "$work/norm.tsv"
"$bin" describe-probe --index "$ord" "$crd" 2>/dev/null > "$work/coord.tsv"
"$RSCRIPT" "$here/compare_cnv.R" "$work/tgt.tsv" "$work/norm.tsv" "$work/coord.tsv" "$work/r_probes.tsv" >/dev/null 2>&1

python3 - "$work/c_probes.tsv" "$work/r_probes.tsv" "$work/c_bins.tsv" <<'PY'
import sys, math
c = {}
for l in open(sys.argv[1]):
    p = l.rstrip("\n").split("\t")
    if p[0] == "sample": continue
    c[p[1]] = float(p[4])
o = {}
for l in open(sys.argv[2]):
    p = l.rstrip("\n").split("\t")
    if p[0] == "Probe_ID": continue
    try: o[p[0]] = float(p[1])
    except ValueError: pass
comm = set(c) & set(o)
if not comm: print("FAIL: no shared probes"); sys.exit(1)
d = sorted(abs(c[k] - o[k]) for k in comm)
mx, med = max(d), d[len(d)//2]

# binning sanity: bins exist, log2 finite, span multiple chromosomes
chroms, nb, badv = set(), 0, 0
for l in open(sys.argv[3]):
    p = l.rstrip("\n").split("\t")
    if p[0] == "sample": continue
    nb += 1; chroms.add(p[1])
    v = p[5]
    if v not in ("nan", "") and not math.isfinite(float(v)): badv += 1

ok = (mx < 1e-3 and med < 1e-5 and nb > 1000 and len(chroms) > 20 and badv == 0)
print(f"cnv OLS+log2 vs R lm: n={len(comm)} max={mx:.2e} median={med:.2e}")
print(f"bins: {nb} over {len(chroms)} chromosomes, {badv} non-finite")
if not ok:
    print(f"FAIL: max={mx:.2e} (<1e-3) med={med:.2e} (<1e-5) bins={nb} chroms={len(chroms)} bad={badv}")
    sys.exit(1)
print(f"ok   cnv: probe log2 matches R lm to {mx:.1e}; {nb} bins over {len(chroms)} chroms")
PY
rc=$?

echo
echo "passed $([ $rc -eq 0 ] && echo 1 || echo 0), failed $([ $rc -eq 0 ] && echo 0 || echo 1)"
[ $rc -eq 0 ]
