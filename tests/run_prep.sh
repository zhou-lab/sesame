#!/bin/sh
# Level-3: each prep step applied in isolation must match R exactly.
set -eu

## The R oracle. Overridable because the binary is not called the same
## thing everywhere: on the lab HPC plain `Rscript` is a different R
## without a usable sesame, and the one to use is Rscript-4.6.0.
RSCRIPT=${RSCRIPT:-Rscript}

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
dump="$root/pipeline_dump"
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }

## Every driver SKIPs rather than fails when a prerequisite is absent: a
## missing oracle is "not measured here", not a regression, and one driver
## exiting non-zero stops `make test` dead for the targets after it.
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "SKIP prep: no $RSCRIPT"; exit 0; }
"$RSCRIPT" -e 'quit(status = !requireNamespace("sesame", quietly = TRUE))' >/dev/null 2>&1 \
    || { echo "SKIP prep: $RSCRIPT has no sesame -- the oracle must be the latest R/Bioc"; exit 0; }
[ -d "$idats" ] || { echo "SKIP: $idats not found"; exit 0; }

fail=0
pass=0

# platform <space> prefix <space> prep-code
while IFS=' ' read -r plat rel code; do
    [ -n "$plat" ] || continue
    idx="$root/testdata/$plat.ordering.tsv.gz"
    pfx="$idats/$rel"
    [ -f "$idx" ] || { echo "SKIP $plat: no $idx (run: make index)"; continue; }

    if ! "$dump" --index "$idx" --prep "$code" --what col "$pfx" \
           > "$work/col.tsv" 2>"$work/c.err"; then
        echo "FAIL $plat $code: sesame errored"; sed 's/^/    /' "$work/c.err" | head -3
        fail=$((fail+1)); continue
    fi
    if "$RSCRIPT" --vanilla "$here/compare_prep.R" "$plat" "$pfx" "$code" \
         "$work/col.tsv" 2>"$work/r.err"; then
        pass=$((pass+1))
    else
        sed 's/^/    /' "$work/r.err" | head -3
        fail=$((fail+1))
    fi
done <<'EOF'
HM450 HM450/3999492009_R01C01 C
EPIC EPIC/GSM2995280_201868590258_R01C01 C
EPICv2 EPICv2/206909630040_R03C01 C
EPICv2 EPICv2/MaleBlood/206909630014_R05C01 C
MSA MSA/207760740030_R01C03 C
EOF

echo
echo "passed $pass, failed $fail"
[ "$fail" -eq 0 ]
