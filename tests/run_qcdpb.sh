#!/bin/sh
# Level-3 (Q + the default QCDPB): sesame vs R on every platform, R given the
# CLI's quality mask (see compare_qcdpb.R for why, and for the pins).
set -eu

## The R oracle. Overridable because the binary is not called the same
## thing everywhere: on the lab HPC plain `Rscript` is a different R
## without a usable sesame, and the one to use is Rscript-4.6.0.
RSCRIPT=${RSCRIPT:-Rscript}

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
dump="$root/pipeline_dump"
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

## Every driver SKIPs rather than fails when a prerequisite is absent: a
## missing oracle is "not measured here", not a regression, and one driver
## exiting non-zero stops `make test` dead for the targets after it.
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "SKIP qcdpb: no $RSCRIPT"; exit 0; }
"$RSCRIPT" -e 'quit(status = !requireNamespace("sesame", quietly = TRUE))' >/dev/null 2>&1 \
    || { echo "SKIP qcdpb: $RSCRIPT has no sesame -- the oracle must be the latest R/Bioc"; exit 0; }

## the probes NA in a dump: $1 = prep
na_ids() {
    YAME_DATA_HOME="$yhome" "$dump" --prep "$1" --what beta "$pfx" 2>/dev/null \
      | python3 -c "import sys; print('\n'.join(sorted(l.split('\t')[0] for l in sys.stdin if l.rstrip('\n').split('\t')[1]=='NA')))"
}

PASS=0; FAIL=0
run_one() {
    plat=$1; rel=$2
    ## SESAME_ONLY: an ERE over "<platform> <prefix>"; cases that do not
    ## match are skipped so a parallel gate can run one case per job.
    if [ -n "${SESAME_ONLY:-}" ] && ! echo "$*" | grep -Eq "$SESAME_ONLY"; then return; fi
    pfx="$idats/$rel"
    ## a platform with no store or IDAT is a FAIL: none may drop out quietly
    [ -d "$store/$plat" ] || {
        echo "FAIL $plat QCDPB: no store $store/$plat"; FAIL=$((FAIL+1)); return; }
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "FAIL $plat QCDPB: no IDAT $pfx"; FAIL=$((FAIL+1)); return; fi

    ## Q's own set: NA under Q, less what is NA with no prep at all
    na_ids Q > "$work/na_Q.txt"
    na_ids "" > "$work/na_.txt"
    LC_ALL=C comm -23 "$work/na_Q.txt" "$work/na_.txt" > "$work/c_q.txt"
    YAME_DATA_HOME="$yhome" "$dump" --prep QCDPB --what beta "$pfx" \
        2>/dev/null > "$work/c_b.txt"

    if "$RSCRIPT" --vanilla "$here/compare_qcdpb.R" "$plat" "$pfx" \
         "$work/c_q.txt" "$work/c_b.txt" 2>"$work/r.err"; then
        PASS=$((PASS+1))
    else
        sed 's/^/    /' "$work/r.err" | head -4
        FAIL=$((FAIL+1))
    fi
}

run_one HM450  HM450/3999492009_R01C01
run_one EPIC   EPIC/GSM2995280_201868590258_R01C01
run_one EPICv2 EPICv2/206909630040_R03C01
run_one MSA    MSA/207760740030_R01C03

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
