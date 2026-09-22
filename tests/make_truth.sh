#!/bin/sh
# tests/make_truth.sh -- regenerate the frozen ground truth in tests/truth/.
#
# For each (test, platform) it runs BOTH sides and writes one row per probe:
# the value, and where it came from -- `R` where the two agree, `C:<code>`
# where they do not or where R cannot compute it at all. A `C:` row without an
# entry in DIVERGENCES.md is refused, because a frozen value nobody can explain
# is a bug with a test around it.
#
#   tests/make_truth.sh [test ...]      default: every test
#
# Needs the store, the test IDATs and an R with the latest sesame/sesameData
# (RSCRIPT=, e.g. Rscript-4.6.0). Never run in CI: the point of the committed
# files is that CI does not need any of this.
set -eu

RSCRIPT=${RSCRIPT:-Rscript}
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
truth="$here/truth"
yhome=${YAME_DATA_HOME:-${XDG_DATA_HOME:-$HOME/.local/share}/yame}
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
bin="$root/sesame"
dump="$root/pipeline_dump"

command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "make_truth: no $RSCRIPT" >&2; exit 1; }
"$RSCRIPT" -e 'quit(status = !requireNamespace("sesame", quietly = TRUE))' >/dev/null 2>&1 \
    || { echo "make_truth: $RSCRIPT has no sesame -- the oracle must be the latest R/Bioc" >&2; exit 1; }
[ -x "$bin" ]  || { echo "make_truth: $bin not built" >&2; exit 1; }
[ -x "$dump" ] || { echo "make_truth: $dump not built (make pipeline_dump)" >&2; exit 1; }

mkdir -p "$truth"

## Provenance: the four versions without which "matches R" means nothing.
prov="$truth/PROVENANCE.tsv"
{
  printf 'generated\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'store\t%s\n' "$yhome"
  "$RSCRIPT" -e 'cat(sprintf("R\t%s\n", R.version.string));
     for (p in c("sesame","sesameData","BiocVersion"))
       cat(sprintf("%s\t%s\n", p, tryCatch(as.character(packageVersion(p)),
                                           error=function(e) "NOT INSTALLED")))' 2>/dev/null
  printf 'sesame-cli\t%s\n' "$("$bin" version 2>/dev/null | head -1 | cut -d' ' -f2)"
  printf 'yame\t%s\n' "$("$bin" version 2>/dev/null | sed -n 's/.*YAME \(v[0-9.]*\).*/\1/p')"
} > "$prov"
echo "make_truth: provenance ->"
sed 's/^/    /' "$prov"

## ---------------------------------------------------------------- vcf
# One row per SNP probe: GT, GS, PVF. GT and PVF are the exact-gated
# quantities; GS carries a documented deep-tail tolerance and is stored as R's.
make_vcf() {
    plat=EPICv2; rel=EPICv2/206909630042_R08C01
    ord="$yhome/$plat/$plat.ordering.tsv.gz"; snp="$yhome/$plat/$plat.hg38.snp.tsv.gz"
    pfx="$idats/$rel"
    [ -f "$ord" ] && [ -f "$snp" ] || { echo "make_truth: vcf needs $ord and $snp" >&2; return 1; }
    w=$(mktemp -d); trap 'rm -rf "$w"' RETURN 2>/dev/null || true

    "$bin" vcf "$pfx" --platform "$ord" --snp "$snp" 2>/dev/null | grep -v '^#' | python3 -c "
import sys
print('Probe_ID\tGT\tGS\tPVF')
for l in sys.stdin:
    p=l.rstrip().split('\t'); d=dict(kv.split('=',1) for kv in p[7].split(';') if '=' in kv)
    print(d['Probe_ID'],d['GT'],p[5],d['PVF'],sep='\t')" > "$w/c.tsv"
    "$RSCRIPT" "$here/compare_vcf.R" "$pfx" $plat "$snp" "$w/r.tsv" >/dev/null 2>&1

    python3 "$here/merge_truth.py" vcf "$w/c.tsv" "$w/r.tsv" "$ord" \
            "$truth/vcf.$plat.tsv" || { rm -rf "$w"; return 1; }
    gzip -f "$truth/vcf.$plat.tsv"; rm -rf "$w"
}

## ---------------------------------------------------------------- pneg
# The negative-control ECDF detection p. R cannot compute it on EPIC or HM450
# with sesameData 1.29.10 (no usable negative controls in its SigDF), so those
# platforms are wholly C-sourced, with the reason recorded.
make_pneg() {
    for pair in "EPICv2 EPICv2/206909630040_R03C01" "MSA MSA/207760740030_R01C03" \
                "EPIC EPIC/GSM2995280_201868590258_R01C01" "HM450 HM450/3999492009_R01C01"; do
        set -- $pair; plat=$1; rel=$2
        ord="$yhome/$plat/$plat.ordering.tsv.gz"; pfx="$idats/$rel"
        [ -f "$ord" ] || { echo "  skip $plat: no ordering"; continue; }
        w=$(mktemp -d)
        YAME_DATA_HOME="$yhome" "$dump" --prep "" --what pval --detection pneg "$pfx" \
            2>/dev/null | /usr/bin/awk 'BEGIN{OFS="\t"; print "Probe_ID","pval"} {print $1,$2}' > "$w/c.tsv"
        ## R's failure must be THE documented one. `|| true` here would let any
        ## failure -- a loaded node, a missing package, an interrupted run --
        ## fall through to "R cannot compute it" and freeze our output under a
        ## code that says something untrue. It did exactly that once: EPICv2
        ## came out C-sourced when it agrees with R to 0.00e+00.
        rerr="$w/r.err"
        "$RSCRIPT" --vanilla - "$plat" "$pfx" "$w/r.tsv" <<'PY' >/dev/null 2>"$rerr" || true
suppressMessages(library(sesame))
a <- commandArgs(TRUE)
sdf <- readIDATpair(a[2], platform = a[1])
rp  <- detectionPnegEcdf(sdf, return.pval = TRUE)
write.table(data.frame(Probe_ID = names(rp), pval = as.numeric(rp)),
            a[3], sep = "\t", quote = FALSE, row.names = FALSE)
PY
        if [ ! -s "$w/r.tsv" ] && ! grep -q "must have 1 or more non-missing values" "$rerr"; then
            echo "make_truth: $plat: R produced nothing, and not for the documented reason:" >&2
            sed 's/^/    /' "$rerr" | head -4 >&2
            echo "    refusing to freeze sesame's own output as truth here." >&2
            rm -rf "$w"; return 1
        fi
        python3 "$here/merge_truth.py" pneg "$w/c.tsv" "$w/r.tsv" "" \
                "$truth/pneg.$plat.tsv" || { rm -rf "$w"; return 1; }
        gzip -f "$truth/pneg.$plat.tsv"; rm -rf "$w"
    done
}

for t in ${*:-vcf pneg}; do
    echo "make_truth: $t"
    make_$t || { echo "make_truth: $t FAILED" >&2; exit 1; }
done
echo "make_truth: done; review the diff before committing"
