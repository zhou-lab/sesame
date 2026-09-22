#!/bin/sh
# vcf: validate SNP genotyping (sesame formatVCF) against R on real EPICv2 data.
# GT (the call) and PVF (variant fraction) must match R exactly; GS (the quality
# score) matches except in the deep tail, where R's dbinom is more precise than
# the log-pmf here. Needs an EPICv2 ordering, the SNP annotation, an EPICv2 IDAT,
# and Rscript with sesame; skips cleanly otherwise.
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
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "SKIP vcf: no Rscript"; exit 0; }

plat=EPICv2
ord="$store/$plat/$plat.ordering.tsv.gz"
[ -f "$ord" ] || ord="$root/testdata/$plat.ordering.tsv.gz"
pfx="$idats/$plat/206909630042_R08C01"

# locate the SNP annotation (not yet in the fetch set): store, $SESAME_SNP_DIR, or a sibling InfiniumAnnotation checkout
snp=""
for c in "$store/$plat/$plat.hg38.snp.tsv.gz" \
         "${SESAME_SNP_DIR:-}/$plat.hg38.snp.tsv.gz" \
         "$HOME/repo/InfiniumAnnotation/$plat/$plat.hg38.snp.tsv.gz"; do
    [ -f "$c" ] && { snp="$c"; break; }
done
[ -n "$snp" ] || { echo "SKIP vcf: no $plat.hg38.snp.tsv.gz (set SESAME_SNP_DIR)"; exit 0; }
[ -f "$ord" ] || { echo "SKIP vcf: no ordering $ord"; exit 0; }
if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
    echo "SKIP vcf: no IDAT $pfx"; exit 0; fi

# C: full VCF, parse INFO -> Probe_ID GT GS PVF
"$bin" vcf "$pfx" --index "$ord" --snp "$snp" 2>/dev/null | grep -v '^#' | \
  python3 -c "
import sys
print('Probe_ID\tGT\tGS\tPVF')
for l in sys.stdin:
    p=l.rstrip().split('\t'); d=dict(kv.split('=',1) for kv in p[7].split(';') if '=' in kv)
    print(d['Probe_ID'],d['GT'],p[5],d['PVF'],sep='\t')" > "$work/c.tsv"

# R oracle
"$RSCRIPT" "$here/compare_vcf.R" "$pfx" $plat "$snp" "$work/r.tsv" >/dev/null 2>&1

python3 - "$work/c.tsv" "$work/r.tsv" "$ord" <<'PY'
import sys
import gzip
def load(fn):
    d={}
    for l in open(fn):
        p=l.rstrip('\n').split('\t')
        if p[0]=='Probe_ID': continue
        d[p[0]]=p[1:]                       # GT, GS, PVF [, Rcol]
    return d
C=load(sys.argv[1]); R=load(sys.argv[2])
common=set(C)&set(R)
if not common: print("FAIL: no shared probes"); sys.exit(1)

## The two sides read DIFFERENT manifests: C the store's ordering, R the one
## inside sesameData. Where they disagree about a Type-I probe's colour
## channel, getAFTypeIbySumAlleles is computed from the other channel and the
## fraction comes out as the exact complement -- an annotation-version
## difference, not a divergence in the port (the formula is identical; see
## NUMERICS.md). Exclude those probes from the gate, but only on the evidence
## of the two manifests, never on whether the answer matched, and report how
## many so drift shows up as a CHANGED count rather than a wider blind spot.
##
## WHICH MANIFEST IS RIGHT IS DELIBERATELY NOT DECIDED HERE. This test asks
## one question -- does the C port reproduce formatVCF -- and answering it does
## not require knowing whose channel assignment is correct. Adjudicating would
## bind the test to re-adjudicate on every annotation release, and would put a
## claim about Illumina's manifest inside a fidelity gate, which is the wrong
## place for it. The excluded probes are listed by the harness; anyone who
## wants the answer can take it up with the manifests themselves.
ocol={}
with gzip.open(sys.argv[3],'rt') as f:
    next(f)
    for l in f:
        p=l.rstrip('\n').split('\t'); ocol[p[0]]=p[3]
excl={k for k in common
      if len(R[k])>3 and R[k][3] and R[k][3]!='NA' and ocol.get(k) and R[k][3]!=ocol[k]}
EXPECTED_EXCL=22          # EPICv2 v8.1 ordering vs sesameData 1.29.10, 2026-09-20
graded=common-excl

gt_mis=0; pvf_max=0.0; gs_off=0
for k in graded:
    cg,cs,cp=C[k][0],C[k][1],C[k][2]; rg,rs,rp=R[k][0],R[k][1],R[k][2]
    if cg!=rg: gt_mis+=1
    try:
        if abs(int(cs)-int(rs))>1: gs_off+=1
    except ValueError: pass
    try: pvf_max=max(pvf_max, abs(float(cp)-float(rp)))
    except ValueError: pass
n=len(graded)
print(f"vcf vs R formatVCF: n={n} GT-mismatch={gt_mis} PVF-maxdiff={pvf_max:.1e} "
      f"GS-off(>1)={gs_off}  [excluded {len(excl)} channel-disagreeing probes]")
ok = (gt_mis==0 and pvf_max < 1e-6 and gs_off <= n*0.001)
if not ok:
    print(f"FAIL: GT={gt_mis} PVF={pvf_max:.1e} GS-off={gs_off} (>{int(n*0.001)})"); sys.exit(1)
if len(excl)!=EXPECTED_EXCL:
    print(f"FAIL: {len(excl)} probes have a channel disagreement between the store "
          f"ordering and sesameData, expected {EXPECTED_EXCL}. The manifests moved; "
          f"re-check which is right and update EXPECTED_EXCL with a note in NUMERICS.md.")
    sys.exit(1)
print(f"ok   vcf: {n} probes, GT + PVF exact vs R, GS off on {gs_off} (deep-tail dbinom); "
      f"{len(excl)} excluded on a manifest channel disagreement")
PY
rc=$?

echo
echo "passed $([ $rc -eq 0 ] && echo 1 || echo 0), failed $([ $rc -eq 0 ] && echo 0 || echo 1)"
[ $rc -eq 0 ]
