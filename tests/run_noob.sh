#!/bin/sh
# Level-3 (B): noob. Unit-tests normExpSignal + MASS::huber against the C
# primitives (must agree to a few ULP), then compares full betas after --prep B
# to getBetas(noob(sdf)) (lineage-bounded, like Q/P/D). See compare_noob.R.
#
# Needs a fetched store with the platform's .cm, the R oracle, and the
# normexp_test unit harness (built by `make test-noob`). Skips cleanly if the
# store or IDATs are missing.
set -eu

## The R oracle. Overridable because the binary is not called the same
## thing everywhere: on the lab HPC plain `Rscript` is a different R
## without a usable sesame, and the one to use is Rscript-4.6.0.
RSCRIPT=${RSCRIPT:-Rscript}

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
dump="$root/pipeline_dump"
unit="$root/normexp_test"
## The shared store yame fetch fills. It is keyed on the browser
## hierarchy -- species/platform -- so a platform's assets sit directly
## at $store/$plat/ and a genome build's at $store/<build>/.
yhome=${YAME_DATA_HOME:-${XDG_DATA_HOME:-$HOME/.local/share}/yame}
store=$yhome
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ]  || { echo "FAIL: $bin not built"; exit 1; }
[ -x "$unit" ] || { echo "FAIL: $unit not built (make test-noob)"; exit 1; }

## Every driver SKIPs rather than fails when a prerequisite is absent: a
## missing oracle is "not measured here", not a regression, and one driver
## exiting non-zero stops `make test` dead for the targets after it.
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "SKIP noob: no $RSCRIPT"; exit 0; }
"$RSCRIPT" -e 'quit(status = !requireNamespace("sesame", quietly = TRUE))' >/dev/null 2>&1 \
    || { echo "SKIP noob: $RSCRIPT has no sesame -- the oracle must be the latest R/Bioc"; exit 0; }

PASS=0; FAIL=0
run_one() {
    plat=$1; rel=$2
    pfx="$idats/$rel"
    [ -d "$store/$plat" ] || { echo "SKIP $plat B: no store $store/$plat"; return; }
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "SKIP $plat B: no IDAT $pfx"; return; fi

    # raw and noob betas, so the oracle can isolate B on raw-identical probes
    YAME_DATA_HOME="$yhome" "$dump" --prep "" --what beta "$pfx" 2>/dev/null > "$work/c_raw.txt"
    YAME_DATA_HOME="$yhome" "$dump" --prep B --what beta "$pfx" 2>/dev/null > "$work/c_noob.txt"

    if "$RSCRIPT" --vanilla "$here/compare_noob.R" "$plat" "$pfx" \
         "$work/c_noob.txt" "$unit" "$work/c_raw.txt" 2>"$work/r.err"; then
        PASS=$((PASS+1))
    else
        sed 's/^/    /' "$work/r.err" | head -8
        FAIL=$((FAIL+1))
    fi
}

run_one MSA MSA/207760740030_R01C03

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
