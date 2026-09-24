#!/bin/sh
# describe-probe, both modes: prepend the ordering's Probe_ID to a positional
# file's rows, and resolve probe IDs to <chrm>_<beg1> coordinates.
# Fully self-contained -- builds a tiny ordering, a format-3 .cg (via mu2cg), a
# text table and a coordinate table, then checks the labeled output, the
# lineage-mismatch guard, replicate expansion and the dropped-probe counts. No
# IDATs, no store, no R oracle.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
mu2cg="$root/mu2cg"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ]   || { echo "FAIL: $bin not built"; exit 1; }
[ -x "$mu2cg" ] || { echo "SKIP describe: $mu2cg not built"; exit 0; }

# --- fixtures: 5 probes, a 5-col ordering (mask inline) and a 4-col one -------
printf 'Probe_ID\tM\tU\tcol\tmask\n' > "$work/ord5.tsv"
printf 'Probe_ID\tM\tU\tcol\n'       > "$work/ord4.tsv"
i=1
while [ $i -le 5 ]; do
    printf 'cg%07d_BC21\tNA\t%d\t2\t0\n' "$i" "$((1000+i))" >> "$work/ord5.tsv"
    printf 'cg%07d_BC21\tNA\t%d\t2\n'    "$i" "$((1000+i))" >> "$work/ord4.tsv"
    i=$((i+1))
done
gzip -f "$work/ord5.tsv" "$work/ord4.tsv"

# a format-3 .cg from a tiny M/U table (one sample "S1")
printf 'Probe_ID\tS1_M\tS1_U\n' > "$work/mu.tsv"
i=1
while [ $i -le 5 ]; do
    printf 'cg%07d_BC21\t%d\t%d\n' "$i" "$((100*i))" "$((10*i))" >> "$work/mu.tsv"
    i=$((i+1))
done
"$mu2cg" "$work/mu.tsv" "$work/s1.cg" >/dev/null 2>&1

# a text table (header + 5 data rows), coord-style, no Probe_ID column
printf 'chrm\tpos\n' > "$work/coord.tsv"
i=1
while [ $i -le 5 ]; do printf 'chr1\t%d\n' "$((5000+i))" >> "$work/coord.tsv"; i=$((i+1)); done

fail=0

# 1) text: Probe_ID prepended, header kept, order preserved
"$bin" describe-probe --index "$work/ord5.tsv.gz" "$work/coord.tsv" > "$work/coord.out" 2>/dev/null
head1=$(sed -n '1p' "$work/coord.out")
[ "$head1" = "$(printf 'Probe_ID\tchrm\tpos')" ] || { echo "FAIL: text header '$head1'"; fail=1; }
[ "$(sed -n '2p' "$work/coord.out")" = "$(printf 'cg0000001_BC21\tchr1\t5001')" ] || { echo "FAIL: text row1"; fail=1; }
[ "$(wc -l < "$work/coord.out")" -eq 6 ] || { echo "FAIL: text nrow"; fail=1; }

# 2) text against a 4-column ordering (mask moved to .cm) parses the same
"$bin" describe-probe --index "$work/ord4.tsv.gz" "$work/coord.tsv" > "$work/coord4.out" 2>/dev/null
cmp -s "$work/coord.out" "$work/coord4.out" || { echo "FAIL: 4-col vs 5-col ordering differ"; fail=1; }

# 3) fmt3 .cg default -> M<TAB>U columns, header "<name>_M <name>_U"
"$bin" describe-probe --index "$work/ord5.tsv.gz" "$work/s1.cg" > "$work/cg.out" 2>/dev/null
[ "$(sed -n '1p' "$work/cg.out")" = "$(printf 'Probe_ID\tS1_M\tS1_U')" ] || { echo "FAIL: cg header"; fail=1; }
[ "$(sed -n '2p' "$work/cg.out")" = "$(printf 'cg0000001_BC21\t100\t10')" ] || { echo "FAIL: cg M/U row1"; fail=1; }

# 4) fmt3 .cg --beta -> beta = M/(M+U); probe 1 = 100/110
"$bin" describe-probe --beta --index "$work/ord5.tsv.gz" "$work/s1.cg" > "$work/beta.out" 2>/dev/null
b=$(sed -n '2p' "$work/beta.out" | cut -f2)
case "$b" in 0.9090*) : ;; *) echo "FAIL: cg beta '$b' != 0.9090..."; fail=1 ;; esac

# 5) lineage mismatch: a 4-probe ordering vs the 5-row file must fail, no stdout
printf 'Probe_ID\tM\tU\tcol\n' > "$work/ord_short.tsv"
i=1; while [ $i -le 4 ]; do printf 'cgX%06d_BC21\tNA\t1\t2\n' "$i" >> "$work/ord_short.tsv"; i=$((i+1)); done
gzip -f "$work/ord_short.tsv"
if "$bin" describe-probe --index "$work/ord_short.tsv.gz" "$work/coord.tsv" > "$work/bad.out" 2>/dev/null; then
    echo "FAIL: mismatch did not error"; fail=1
fi
[ -s "$work/bad.out" ] && { echo "FAIL: mismatch emitted output before erroring"; fail=1; } || true

if [ $fail -eq 0 ]; then
    echo "ok   describe-probe: text + 4/5-col ordering + fmt3 M/U + beta + mismatch guard"
    pass=1
else
    pass=0
fi

# --- coordinate mode ---------------------------------------------------------
# The ordering above is cg0000001_BC21 .. cg0000005_BC21, plus one replicate so
# a bare cg number has to expand. The coord table is positional over it and
# 0-BASED, the way the store's <plat>.<genome>.coord.tsv.gz is.
cfail=0
printf 'Probe_ID\tM\tU\tcol\n' > "$work/ordc.tsv"
i=1
while [ $i -le 5 ]; do
    printf 'cg%07d_BC21\tNA\t%d\t2\n' "$i" "$((1000+i))" >> "$work/ordc.tsv"
    i=$((i+1))
done
printf 'cg0000001_BC22\tNA\t9\t2\n' >> "$work/ordc.tsv"   # replicate of probe 1
gzip -f "$work/ordc.tsv"

## row 2 has no CpG in this genome, row 3 sits on an alt contig: both are
## resolvable probes that a genome-indexed store cannot address
printf 'CpG_chrm\tCpG_beg\tstrand\tmapQ\n' > "$work/coords.tsv"
printf 'chr1\t999\t+\t60\n'                  >> "$work/coords.tsv"
printf 'NA\tNA\t*\t0\n'                      >> "$work/coords.tsv"
printf 'chr1_KI270766v1_alt\t500\t+\t60\n'   >> "$work/coords.tsv"
printf 'chr7\t14\t-\t60\n'                   >> "$work/coords.tsv"
printf 'chrX\t77777\t+\t60\n'                >> "$work/coords.tsv"
printf 'chr1\t999\t+\t60\n'                  >> "$work/coords.tsv"
gzip -f "$work/coords.tsv"

describe() {
    "$bin" describe-probe --platform "$work/ordc.tsv.gz" --coords "$work/coords.tsv.gz" \
        --genome testgenome "$@"
}

## 1) a bare cg number expands to every design suffix, 0-based beg -> 1-based
printf 'cg0000001\n' > "$work/q1.txt"
describe "$work/q1.txt" > "$work/c1.out" 2>/dev/null
want=$(printf 'cg0000001_BC21\tchr1_1000\ncg0000001_BC22\tchr1_1000\n')
[ "$(cat "$work/c1.out")" = "$want" ] || { echo "FAIL: replicate expansion"; sed 's/^/    /' "$work/c1.out"; cfail=1; }

## 2) the full ID matches only itself
printf 'cg0000001_BC22\n' > "$work/q2.txt"
describe "$work/q2.txt" > "$work/c2.out" 2>/dev/null
[ "$(cat "$work/c2.out")" = "$(printf 'cg0000001_BC22\tchr1_1000')" ] || { echo "FAIL: exact ID"; cfail=1; }

## 3) input order is kept, and no header is emitted (the consumer is `cut -f2`)
printf 'cg0000005\ncg0000004\n' > "$work/q3.txt"
describe "$work/q3.txt" > "$work/c3.out" 2>/dev/null
[ "$(cat "$work/c3.out")" = "$(printf 'cg0000005_BC21\tchrX_77778\ncg0000004_BC21\tchr7_15')" ] \
    || { echo "FAIL: input order / header"; sed 's/^/    /' "$work/c3.out"; cfail=1; }

## 4) no CpG in the genome and a non-primary contig are dropped AND counted
printf 'cg0000002\ncg0000003\n' > "$work/q4.txt"
describe "$work/q4.txt" > "$work/c4.out" 2> "$work/c4.err"
[ -s "$work/c4.out" ] && { echo "FAIL: unmappable probes emitted a row"; cfail=1; }
grep -q '1 with no CpG in testgenome' "$work/c4.err" || { echo "FAIL: no-CpG not counted"; sed 's/^/    /' "$work/c4.err"; cfail=1; }
grep -q '1 on a non-primary contig' "$work/c4.err" || { echo "FAIL: alt contig not counted"; sed 's/^/    /' "$work/c4.err"; cfail=1; }

## 5) stdin, the spelling the pipeline uses
printf 'cg0000004\n' | describe - > "$work/c5.out" 2>/dev/null
[ "$(cat "$work/c5.out")" = "$(printf 'cg0000004_BC21\tchr7_15')" ] || { echo "FAIL: stdin"; cfail=1; }

## 6) an ID the platform does not carry is an error, not a short file
printf 'cg0000004\ncg9999999\n' > "$work/q6.txt"
if describe "$work/q6.txt" > "$work/c6.out" 2>/dev/null; then
    echo "FAIL: unknown probe did not error"; cfail=1
fi

## 7) the old name still works, and says so
printf 'cg0000004\n' > "$work/q7.txt"
"$bin" attach-probe --platform "$work/ordc.tsv.gz" --coords "$work/coords.tsv.gz" \
    --genome testgenome "$work/q7.txt" > "$work/c7.out" 2> "$work/c7.err"
grep -q 'attach-probe is now describe-probe' "$work/c7.err" || { echo "FAIL: alias is silent"; cfail=1; }
[ "$(cat "$work/c7.out")" = "$(printf 'cg0000004_BC21\tchr7_15')" ] || { echo "FAIL: alias output"; cfail=1; }

if [ $cfail -eq 0 ]; then
    echo "ok   describe-probe --genome: replicates + 1-based + drops counted + stdin + alias"
    pass=$((pass+1))
else
    fail=1
fi

echo
echo "passed $pass, failed $(( (fail!=0) + 0 ))"
[ $fail -eq 0 ] && [ $cfail -eq 0 ]
