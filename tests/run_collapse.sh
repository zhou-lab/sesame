#!/bin/sh
# betasCollapseToPfx: preprocess --collapse averages replicate probes to their
# cg-prefix. To test the collapse arithmetic in isolation from every pre-existing
# raw/ordering divergence, R's betasCollapseToPfx is run on C's OWN uncollapsed
# raw betas (pipeline_dump --prep "" --what beta) and compared to C's collapsed
# beta.cg -- identical input, so the prefix set and NA pattern must match exactly
# and values must match within the float32 .cg product precision. Needs a store
# with the ordering, the yame + pipeline_dump binaries, the R oracle, and IDATs.
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
[ -x "$yame" ] || { echo "SKIP collapse: no $yame"; exit 0; }
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "SKIP collapse: no Rscript"; exit 0; }

PASS=0; FAIL=0
run_one() {
    plat=$1; rel=$2
    ## SESAME_ONLY: an ERE over "<platform> <prefix> ..."; cases that do not
    ## match are skipped so a parallel gate can run one case per job.
    if [ -n "${SESAME_ONLY:-}" ] && ! echo "$*" | grep -Eq "$SESAME_ONLY"; then return; fi
    ord="$store/$plat/$plat.ordering.tsv.gz"
    pfx="$idats/$rel"
    [ -f "$ord" ] || { echo "SKIP $plat collapse: no $ord"; return; }
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "SKIP $plat collapse: no IDAT $pfx"; return; fi

    YAME_DATA_HOME="$yhome" "$dump" --prep "" --what beta "$pfx" 2>/dev/null > "$work/raw.txt"
    YAME_DATA_HOME="$yhome" "$bin" preprocess --platform "$plat" --output beta \
        --prep "" --collapse --out "$work/pp" "$pfx" 2>/dev/null
    "$yame" unpack "$work/pp/beta.cg" 2>/dev/null > "$work/vals.txt"
    paste "$work/pp/beta.cg.probes" "$work/vals.txt" > "$work/c.txt"

    if "$RSCRIPT" --vanilla - "$work/raw.txt" "$work/c.txt" <<'PY' 2>"$work/r.err"
suppressMessages(library(sesame))
a <- commandArgs(TRUE)
rv <- read.table(a[1], colClasses=c("character","numeric"))
bc <- betasCollapseToPfx(setNames(rv$V2, rv$V1))   # R collapse of C's own betas
cc <- read.table(a[2], colClasses=c("character","numeric")); cC<-setNames(cc$V2,cc$V1)
extra <- length(setdiff(names(bc),names(cC))) + length(setdiff(names(cC),names(bc)))
ids <- intersect(names(bc), names(cC))
namis <- sum(xor(is.na(bc[ids]), is.na(cC[ids])))
d <- abs(bc[ids]-cC[ids]); d <- d[!is.na(d)]
mx <- if (length(d)) max(d) else 0
cat(sprintf("ok   collapse: %d prefixes, set-diff=%d NAmis=%d max|diff|=%.2e\n",
            length(ids), extra, namis, mx))
if (extra != 0 || namis != 0 || mx > 2e-3) { cat("FAIL: collapse diverges\n"); quit(status=1) }
PY
    then PASS=$((PASS+1)); else echo "  ($plat)"; sed 's/^/    /' "$work/r.err" | head -4; FAIL=$((FAIL+1)); fi
}

run_one EPICv2 EPICv2/206909630040_R03C01
run_one MSA    MSA/207760740030_R01C03

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
