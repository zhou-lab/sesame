#!/bin/sh
# imputeBetasMatrixByMean: `impute --method mean` fills NA with a probe's mean
# across samples (--axis probe, R axis=1) or a sample's mean across probes
# (--axis sample, R axis=2). Compared to R's imputeBetasMatrixByMean on C's own
# betas matrix (read back via describe-probe --all) -- exact up to float32 .cg
# rounding. Needs a store with the ordering, yame's describe path, R, and >=2 IDATs.
#
# (impute --method neighbors is NOT gated here: it needs R's empirically-mapped
# manifest coordinates, which the store's design-coord table does not provide.)
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
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "SKIP impute: no Rscript"; exit 0; }

plat=EPICv2
ord="$store/$plat/$plat.ordering.tsv.gz"
p1="$idats/$plat/206909630040_R03C01"
p2="$idats/$plat/206909630042_R08C01"
[ -f "$ord" ] || { echo "SKIP impute: no $ord"; exit 0; }
for p in "$p1" "$p2"; do
    if [ ! -f "$p"_Grn.idat ] && [ ! -f "$p"_Grn.idat.gz ]; then
        echo "SKIP impute: no IDAT $p"; exit 0; fi
done

YAME_DATA_HOME="$yhome" "$bin" preprocess --platform $plat --output beta \
    --out "$work/pp" "$p1" "$p2" 2>/dev/null
"$bin" describe-probe --all --index "$ord" "$work/pp/beta.cg" 2>/dev/null > "$work/in.tsv"

PASS=0; FAIL=0
one_axis() {
    ax=$1; rax=$2
    "$bin" impute --method mean --axis "$ax" "$work/pp/beta.cg" "$work/out.cg" 2>/dev/null
    "$bin" describe-probe --all --index "$ord" "$work/out.cg" 2>/dev/null > "$work/out.tsv"
    if "$RSCRIPT" --vanilla - "$work/in.tsv" "$work/out.tsv" "$rax" "$ax" <<'PY' 2>"$work/r.err"
suppressMessages(library(sesame)); a<-commandArgs(TRUE)
inp<-as.matrix(read.table(a[1],header=TRUE,sep="\t",row.names=1,check.names=FALSE))
out<-as.matrix(read.table(a[2],header=TRUE,sep="\t",row.names=1,check.names=FALSE))
imp<-imputeBetasMatrixByMean(inp, axis=as.integer(a[3])); co<-out[rownames(imp),colnames(imp)]
d<-abs(imp-co); d<-d[!is.na(d)]
cat(sprintf("ok   mean axis=%-6s: max|diff|=%.2e NAmis=%d\n", a[4],
            if(length(d))max(d)else 0, sum(xor(is.na(imp),is.na(co)))))
if ((length(d)&&max(d)>1e-5) || sum(xor(is.na(imp),is.na(co)))!=0) { cat("FAIL\n"); quit(status=1) }
PY
    then PASS=$((PASS+1)); else sed 's/^/    /' "$work/r.err" | head -3; FAIL=$((FAIL+1)); fi
}

one_axis probe  1
one_axis sample 2

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
