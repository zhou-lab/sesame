#!/bin/sh
# mLiftOver: lift a beta.cg across platforms by the probe-ID prefix join. Tested
# by lifting C's OWN raw betas and comparing to R's mLiftOver of the same source
# vector -- so the mapping (target set, order, first-match choice, NA pattern) is
# checked exactly, independent of any raw/prep divergence. Values must match
# within the float32 .cg product precision. Needs source+target orderings (store
# or testdata/), the yame + pipeline_dump binaries, the R oracle, and IDATs.
set -eu

## The R oracle. Overridable because the binary is not called the same
## thing everywhere: on the lab HPC plain `Rscript` is a different R
## without a usable sesame, and the one to use is Rscript-4.6.0.
RSCRIPT=${RSCRIPT:-Rscript}

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
dump="$root/pipeline_dump"
yame="$root/YAME/yame"
## The shared store yame fetch fills. It is keyed on the browser
## hierarchy -- species/platform -- so a platform's assets sit directly
## at $store/$plat/ and a genome build's at $store/<build>/.
yhome=${YAME_DATA_HOME:-${XDG_DATA_HOME:-$HOME/.local/share}/yame}
store=$yhome
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ]  || { echo "FAIL: $bin not built"; exit 1; }
[ -x "$dump" ] || { echo "FAIL: $dump not built (make pipeline_dump)"; exit 1; }
[ -x "$yame" ] || { echo "SKIP mLiftOver: no $yame"; exit 0; }
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "SKIP mLiftOver: no Rscript"; exit 0; }

find_ord() {   # echo the first existing ordering for platform $1
    for c in "$store/$1/$1.ordering.tsv.gz" "$root/testdata/$1.ordering.tsv.gz"; do
        [ -f "$c" ] && { echo "$c"; return 0; }
    done
    ## Explicitly 0: without it the last failed `[ -f ]` is this function's
    ## status, and `so=$(find_ord ...)` then kills the script under `set -e` --
    ## which is how `make test` aborted at test-liftover with no message and
    ## skipped the twelve targets after it, instead of printing SKIP like every
    ## other driver.
    return 0
}

PASS=0; FAIL=0
run_pair() {
    sp=$1; rel=$2; tp=$3
    so=$(find_ord "$sp"); to=$(find_ord "$tp")
    pfx="$idats/$rel"
    [ -n "$so" ] && [ -n "$to" ] || { echo "SKIP $sp->$tp: missing ordering"; return; }
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "SKIP $sp->$tp: no IDAT $pfx"; return; fi

    "$bin" preprocess --platform "$sp" --index "$so" --output beta --prep "" \
        --out "$work/pp" "$pfx" 2>/dev/null
    "$bin" mliftover --to "$tp" --platform "$sp" --index "$so" --index-to "$to" \
        "$work/pp/beta.cg" "$work/out.cg" 2>/dev/null
    "$yame" unpack "$work/out.cg" 2>/dev/null > "$work/vals.txt"
    zcat < "$to" | tail -n +2 | cut -f1 > "$work/tids.txt"
    paste "$work/tids.txt" "$work/vals.txt" > "$work/c.txt"
    "$dump" --index "$so" --prep "" --what beta "$pfx" 2>/dev/null > "$work/src.txt"

    if "$RSCRIPT" --vanilla - "$sp" "$tp" "$work/src.txt" "$work/c.txt" <<'PY' 2>"$work/r.err"
suppressMessages(library(sesame))
a <- commandArgs(TRUE); sp<-a[1]; tp<-a[2]
sv <- read.table(a[3], colClasses=c("character","numeric"))
rl <- mLiftOver(setNames(sv$V2, sv$V1), tp, sp)
cc <- read.table(a[4], colClasses=c("character","numeric")); cL<-setNames(cc$V2,cc$V1)

## Control probes are named differently by the two orderings, and only by
## name: the store's EPIC ordering calls them ctl_<address>_<TYPE>
## (ctl_10609447_NEGATIVE), sesameData's EPIC.address$ordering calls them
## ctl_<address>. Measured 2026-09-20: 635 such probes on EPIC, and with the
## type suffix stripped the two sets are IDENTICAL -- same addresses, none of
## them a real probe. Comparing the raw strings reported 1270 set differences
## for a difference in spelling, so normalise before comparing. A real probe
## missing from either side still fails, which is what this gate is for.
ctl_norm <- function(v) sub("^(ctl_[0-9]+)_.*$", "\\1", v)
names(rl) <- ctl_norm(names(rl)); names(cL) <- ctl_norm(names(cL))
rl <- rl[!duplicated(names(rl))]; cL <- cL[!duplicated(names(cL))]

ro <- length(setdiff(names(rl),names(cL))); co <- length(setdiff(names(cL),names(rl)))
ids <- intersect(names(rl), names(cL))
namis <- sum(xor(is.na(rl[ids]), is.na(cL[ids])))
d <- abs(rl[ids]-cL[ids]); d <- d[!is.na(d)]; mx <- if(length(d)) max(d) else 0
cat(sprintf("ok   %-6s -> %-6s: %d probes, set-diff=%d NAmis=%d max|diff|=%.2e\n",
            sp, tp, length(ids), ro+co, namis, mx))
if (ro+co != 0 || namis != 0 || mx > 2e-3) { cat("FAIL: mLiftOver diverges\n"); quit(status=1) }
PY
    then PASS=$((PASS+1)); else sed 's/^/    /' "$work/r.err" | head -4; FAIL=$((FAIL+1)); fi
}

run_pair EPICv2 EPICv2/206909630040_R03C01           EPIC
run_pair EPIC   EPIC/GSM2995280_201868590258_R01C01  EPICv2

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
