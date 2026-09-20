#!/bin/sh
# Level-3 (QC): the sesameQC panel vs R's sesameQC_calcStats. Every metric must
# agree to within data-lineage scale; a formula bug would be far off. See
# compare_qc.R. Needs a fetched store with the platform's .cm and the R oracle;
# skips cleanly if missing.
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

## Every driver SKIPs rather than fails when a prerequisite is absent: a
## missing oracle is "not measured here", not a regression, and one driver
## exiting non-zero stops `make test` dead for the targets after it.
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "SKIP qc: no $RSCRIPT"; exit 0; }
"$RSCRIPT" -e 'quit(status = !requireNamespace("sesame", quietly = TRUE))' >/dev/null 2>&1 \
    || { echo "SKIP qc: $RSCRIPT has no sesame -- the oracle must be the latest R/Bioc"; exit 0; }

PASS=0; FAIL=0
run_one() {
    plat=$1; rel=$2
    pfx="$idats/$rel"
    [ -d "$store/$plat" ] || { echo "SKIP $plat QC: no store $store/$plat"; return; }
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "SKIP $plat QC: no IDAT $pfx"; return; fi

    YAME_DATA_HOME="$yhome" "$bin" preprocess --output qc --out "$work" "$pfx" 2>/dev/null
    mv "$work/qc.tsv" "$work/c_qc.tsv"

    if "$RSCRIPT" --vanilla "$here/compare_qc.R" "$plat" "$pfx" "$work/c_qc.tsv" \
         2>"$work/r.err"; then
        PASS=$((PASS+1))
    else
        sed 's/^/    /' "$work/r.err" | head -30
        FAIL=$((FAIL+1))
    fi
}

run_one MSA MSA/207760740030_R01C03

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
