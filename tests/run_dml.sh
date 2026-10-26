#!/bin/sh
# Level (DML): per-probe differential methylation vs R's DML. No data-lineage
# caveat -- the betas matrix is the input -- so the OLS estimates, t/F p-values,
# effect sizes, and BH adjustment must match R to ~1e-8. The R oracle generates a
# small dataset, runs `sesame dml`, and compares (see compare_dml.R).
# Before that, and without R: the .cg form takes its Probe_IDs from --platform,
# given as an ordering path or as a name resolved in the store.
set -eu

## The R oracle. Overridable because the binary is not called the same
## thing everywhere: on the lab HPC plain `Rscript` is a different R
## without a usable sesame, and the one to use is Rscript-4.6.0.
RSCRIPT=${RSCRIPT:-Rscript}

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
yame="$root/YAME/yame"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }

## --- the .cg form: --platform <path>, --platform <name>, --index alias ------
## Six samples of five probes, group A ~0.2 and group B ~0.7, packed one
## format-4 record per sample, plus a 5-probe ordering in a throwaway store.
fail=0
if [ -x "$yame" ]; then
    printf 'Probe_ID\tM\tU\tcol\n' > "$work/ord.tsv"
    for p in 1 2 3 4 5; do
        printf 'cg%07d_BC21\tNA\t%d\t2\n' "$p" "$((1000+p))" >> "$work/ord.tsv"
    done
    gzip "$work/ord.tsv"
    mkdir -p "$work/store/EPICv2"
    cp "$work/ord.tsv.gz" "$work/store/EPICv2/EPICv2.ordering.tsv.gz"
    printf 'sample\tgroup\n' > "$work/meta.tsv"
    : > "$work/beta.cg"
    for s in 1 2 3 4 5 6; do
        /usr/bin/awk -v s="$s" 'BEGIN { for (p = 1; p <= 5; p++)
            printf "%.3f\n", (s <= 3 ? 0.2 : 0.7) + 0.01*p + 0.003*s }' \
            > "$work/b.txt"
        "$yame" pack -f n "$work/b.txt" "$work/b.cg" 2>/dev/null
        cat "$work/b.cg" >> "$work/beta.cg"
        echo "S$s" >> "$work/names.txt"
        printf 'S%d\t%s\n' "$s" "$([ "$s" -le 3 ] && echo A || echo B)" \
            >> "$work/meta.tsv"
    done
    "$yame" index -s "$work/names.txt" "$work/beta.cg" 2>/dev/null

    run() { "$bin" dml --betas "$work/beta.cg" --meta "$work/meta.tsv" \
                --formula '~ group' "$@"; }
    if run --platform "$work/ord.tsv.gz" > "$work/p.tsv" 2>"$work/p.err" &&
       [ "$(sed -n '2p' "$work/p.tsv" | cut -f1)" = cg0000001_BC21 ] &&
       [ "$(wc -l < "$work/p.tsv")" -eq 6 ]; then
        echo "ok   dml: .cg with --platform <ordering path>"
    else echo "FAIL dml: .cg with --platform <ordering path>"; fail=1; fi

    if YAME_DATA_HOME="$work/store" run --platform EPICv2 \
           > "$work/n.tsv" 2>/dev/null &&
       cmp -s "$work/p.tsv" "$work/n.tsv"; then
        echo "ok   dml: .cg with --platform <name> from the store"
    else echo "FAIL dml: .cg with --platform <name> from the store"; fail=1; fi

    ## --index is the same ordering under its old name: identical output, and
    ## no notice sending the reader to a flag
    if run --index "$work/ord.tsv.gz" > "$work/i.tsv" 2>"$work/i.err" &&
       cmp -s "$work/p.tsv" "$work/i.tsv" && [ ! -s "$work/i.err" ]; then
        echo "ok   dml: --index is a silent alias of --platform <path>"
    else echo "FAIL dml: --index alias differs or warns"; fail=1; fi

    if run > /dev/null 2>"$work/none.err" ||
       ! grep -q -- '--platform' "$work/none.err"; then
        echo "FAIL dml: a .cg with no ordering must fail naming --platform"
        fail=1
    else echo "ok   dml: a .cg with no ordering fails naming --platform"; fi
else
    echo "SKIP dml .cg form: $yame not built"
fi
[ "$fail" -eq 0 ] || { echo "passed 0, failed 1"; exit 1; }

## Every driver SKIPs rather than fails when a prerequisite is absent: a
## missing oracle is "not measured here", not a regression, and one driver
## exiting non-zero stops `make test` dead for the targets after it.
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "SKIP dml: no $RSCRIPT"; exit 0; }
"$RSCRIPT" -e 'quit(status = !requireNamespace("sesame", quietly = TRUE))' >/dev/null 2>&1 \
    || { echo "SKIP dml: $RSCRIPT has no sesame -- the oracle must be the latest R/Bioc"; exit 0; }

if "$RSCRIPT" --vanilla "$here/compare_dml.R" "$bin"; then
    echo "passed 1, failed 0"
else
    echo "passed 0, failed 1"; exit 1
fi
