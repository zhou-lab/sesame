#!/bin/sh
# Golden level-1: the C IDAT reader must agree with R's readIDAT() bit-for-bit
# on IlluminaID / Mean / SD / NBeads. No tolerance -- these are integers.
#
# Corpora:
#   1. sesameData's extdata (small HM450 subset, always available)
#   2. $SESAME_TEST_IDATS (default ~/repo/InfiniumTestIDATs) -- real
#      full-size arrays across every platform, plain and gzipped.
#
# Requires Rscript with sesame + sesameData installed (the R oracle).
set -eu

## The R oracle. Overridable because the binary is not called the same
## thing everywhere: on the lab HPC plain `Rscript` is a different R
## without a usable sesame, and the one to use is Rscript-4.6.0.
RSCRIPT=${RSCRIPT:-Rscript}

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }

## Every driver SKIPs rather than fails when a prerequisite is absent: a
## missing oracle is "not measured here", not a regression, and one driver
## exiting non-zero stops `make test` dead for the targets after it.
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "SKIP golden: no $RSCRIPT"; exit 0; }
"$RSCRIPT" -e 'quit(status = !requireNamespace("sesame", quietly = TRUE))' >/dev/null 2>&1 \
    || { echo "SKIP golden: $RSCRIPT has no sesame -- the oracle must be the latest R/Bioc"; exit 0; }

fail=0
pass=0

check_one() {
    f=$1
    label=$2

    # R oracle: dump addr/mean/sd/nbeads in file order, no header.
    if ! "$RSCRIPT" --vanilla -e '
        suppressMessages(library(sesame))
        f <- commandArgs(trailingOnly=TRUE)[1]
        r <- suppressWarnings(sesame:::readIDAT(f))
        q <- r$Quants
        cat(sprintf("%s\t%d\t%d\t%d",
            rownames(q), q[,"Mean"], q[,"SD"], q[,"NBeads"]), sep="\n")
    ' "$f" > "$work/r.tsv" 2>"$work/r.err"; then
        echo "FAIL $label: R oracle errored"
        sed 's/^/    /' "$work/r.err" | head -3
        fail=$((fail+1)); return
    fi

    if ! "$bin" idat-dump --tsv "$f" > "$work/c.tsv" 2>"$work/c.err"; then
        echo "FAIL $label: sesame errored"
        sed 's/^/    /' "$work/c.err" | head -3
        fail=$((fail+1)); return
    fi

    if cmp -s "$work/r.tsv" "$work/c.tsv"; then
        n=$(wc -l < "$work/c.tsv" | tr -d ' ')
        printf 'ok   %-52s %8s records\n' "$label" "$n"
        pass=$((pass+1))
    else
        echo "FAIL $label: C and R differ"
        diff "$work/r.tsv" "$work/c.tsv" | head -6 | sed 's/^/    /'
        fail=$((fail+1))
    fi
}

echo "== corpus 1: sesameData extdata =="
extdata=$("$RSCRIPT" -e 'cat(system.file("extdata","",package="sesameData"))' 2>/dev/null || true)
if [ -n "${extdata:-}" ] && [ -d "$extdata" ]; then
    for f in "$extdata"/*.idat; do
        [ -e "$f" ] || continue
        check_one "$f" "$(basename "$f")"
    done
else
    echo "SKIP: sesameData not installed"
fi

echo
echo "== corpus 2: real arrays, all platforms =="
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
if [ -d "$idats" ]; then
    # List first, then loop with a redirect (not a pipe) so the counters stay
    # in this shell rather than a forked subshell.
    find "$idats" \( -iname '*.idat' -o -iname '*.idat.gz' \) 2>/dev/null \
      | grep -v '/\.git/' | sort > "$work/list"
    while IFS= read -r f; do
        check_one "$f" "${f#"$idats"/}"
    done < "$work/list"
else
    echo "SKIP: $idats not found (set SESAME_TEST_IDATS)"
fi

echo
echo "== gz round-trip =="
src=$(ls "$extdata"/*_Grn.idat 2>/dev/null | head -1 || true)
if [ -n "${src:-}" ]; then
    gzip -c "$src" > "$work/t.idat.gz"
    "$bin" idat-dump --tsv "$src"            > "$work/plain.tsv"
    "$bin" idat-dump --tsv "$work/t.idat.gz" > "$work/gz.tsv"
    if cmp -s "$work/plain.tsv" "$work/gz.tsv"; then
        echo "ok   plain == gz"
        pass=$((pass+1))
    else
        echo "FAIL gz round-trip"
        fail=$((fail+1))
    fi
fi

echo
echo "passed $pass, failed $fail"
[ "$fail" -eq 0 ]
