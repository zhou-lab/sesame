#!/bin/sh
# Level-3 (E): dyeBiasL (linear dye bias). E is applied in its natural pipeline
# position -- AFTER channel inference (prep "CE"), where C and R agree on every
# Infinium-I channel, so the fR/fG median pools match and betas are bit-identical
# (up to float rounding). Applying E to the RAW SigDF instead would expose the
# handful of pre-C raw channel disagreements, which C resolves in the C step and
# which are unrelated to E. Compares full betas to R's
# getBetas(dyeBiasL(inferInfiniumIChannel(sdf))).
#
# Needs a store with the platform ordering, the R oracle, and IDATs; skips clean.
set -eu

## The R oracle. Overridable because the binary is not called the same
## thing everywhere: on the lab HPC plain `Rscript` is a different R
## without a usable sesame, and the one to use is Rscript-4.6.0.
RSCRIPT=${RSCRIPT:-Rscript}

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
dump="$root/pipeline_dump"
## The shared store yame fetch fills. It is keyed on the browser
## hierarchy -- species/platform -- so a platform's assets sit directly
## at $store/$plat/ and a genome build's at $store/<build>/.
yhome=${YAME_DATA_HOME:-${XDG_DATA_HOME:-$HOME/.local/share}/yame}
store=$yhome
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$dump" ] || { echo "FAIL: $dump not built (make pipeline_dump)"; exit 1; }
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "SKIP dyebiasL: no Rscript"; exit 0; }

PASS=0; FAIL=0
run_one() {
    plat=$1; rel=$2
    ## SESAME_ONLY: an ERE over "<platform> <prefix> ..."; cases that do not
    ## match are skipped so a parallel gate can run one case per job.
    if [ -n "${SESAME_ONLY:-}" ] && ! echo "$*" | grep -Eq "$SESAME_ONLY"; then return; fi
    ord="$store/$plat/$plat.ordering.tsv.gz"
    pfx="$idats/$rel"
    [ -f "$ord" ] || { echo "SKIP $plat E: no $ord"; return; }
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "SKIP $plat E: no IDAT $pfx"; return; fi

    YAME_DATA_HOME="$yhome" "$dump" --prep CE --what beta "$pfx" \
        2>/dev/null > "$work/c_CE.txt"

    if "$RSCRIPT" --vanilla - "$plat" "$pfx" "$work/c_CE.txt" <<'PY' 2>"$work/r.err"
suppressMessages(library(sesame))
a <- commandArgs(TRUE); plat<-a[1]; pfx<-a[2]; cf<-a[3]

## R reads sesameData's manifest, C reads the store's InfiniumAnnotation
## ordering -- deliberately, because that IS the shipping configuration and the
## question is whether the annotation difference moves the biology. The two
## disagree on a few probes' design type (e.g. EPIC/HM450 cg07162498,
## cg09334382: Infinium-II in sesameData, Infinium-I in v8.1, which carries the
## corrected Illumina manifest). Those probes legitimately read differently, and
## they also sit in the Infinium-I median pools on one side only, nudging E's
## fG/fR for everything else.
##
## So this is NOT a bit-identity gate. It asserts the shape we actually care
## about: the overwhelming majority of probes agree to well within any
## biological resolution, and only a countable handful differ. A real E bug
## moves the whole distribution (an earlier channel-swap bug moved betas 0.4-0.7),
## which p99 catches immediately; an annotation difference moves a few probes.
sdf <- readIDATpair(pfx, platform=plat)
rB  <- getBetas(dyeBiasL(inferInfiniumIChannel(sdf)))
cb  <- read.table(cf, colClasses=c("character","numeric")); cB<-setNames(cb$V2,cb$V1)
ids <- intersect(names(rB), names(cB))
d   <- abs(rB[ids]-cB[ids]); d <- d[!is.na(d)]
nam <- sum(xor(is.na(rB[ids]), is.na(cB[ids])))
n   <- length(d)
qs  <- quantile(d, c(0.5, 0.99, 0.9999))
n01 <- sum(d > 0.01)          # biologically meaningful disagreement
n1e3<- sum(d > 1e-3)
cat(sprintf("ok   %-8s E: n=%d  median=%.1e p99=%.1e p99.99=%.1e  >1e-3: %d  >0.01: %d (%.4f%%)  NA-mism=%d\n",
            plat, n, qs[1], qs[2], qs[3], n1e3, n01, 100*n01/n, nam))

## Majority agreement, not equality. The bar is biological, not numerical:
## p99 below 1e-3 is two orders under array reproducibility (~1e-2), so a shift
## of that size cannot change a call. HM450 sits at ~1e-4 across the board
## because the two reclassified probes are in its Infinium-I median pools on one
## side only, moving fG/fR very slightly for every probe -- annotation, not a
## bug. A real E fault moves the distribution by 0.1+ and fails on p99 alone.
ok <- qs[2] < 1e-3 && n01 <= 25 && (n01/n) < 1e-4 && nam == 0
if (!ok) {
    cat(sprintf("FAIL %s: E disagrees beyond annotation lineage (p99=%.2e, >0.01: %d)\n",
                plat, qs[2], n01))
    quit(status=1)
}
PY
    then PASS=$((PASS+1)); else sed 's/^/    /' "$work/r.err" | head -4; FAIL=$((FAIL+1)); fi
}

run_one EPICv2 EPICv2/206909630040_R03C01
run_one MSA    MSA/207760740030_R01C03
run_one EPIC   EPIC/GSM2995280_201868590258_R01C01
run_one HM450  HM450/3999492009_R01C01

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
