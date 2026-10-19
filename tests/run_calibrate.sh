#!/bin/sh
# sesame calibrate: the titration inversion, checked value by value against
# the formula on a hand-made table and beta.cg. Self-contained (no store, no
# IDATs, no R), so it runs in CI too.
#
#   b' = (beta - floor) / max(top - floor, 0.05)
#   f  = b' / (b' + rho^gamma (1 - b')),   clamped to [0, 1]
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
yame="$root/YAME/yame"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }
[ -x "$yame" ] || { echo "SKIP calibrate: no $yame"; exit 0; }

## one row per case: floor top rho | beta | expected f at gamma 0.15
##  1 mid-curve, rho 8          2 linear probe (rho 1)
##  3 below floor -> 0          4 above top -> 1
##  5 flat probe: amp 0.05      6 NA fit -> NA       7 NA beta -> NA
cat > "$work/tab.tsv" <<'TAB'
floor	top	rho
0.1	0.9	8
0.0	1.0	1
0.2	0.8	4
0.2	0.8	4
0.30	0.32	2
NA	NA	NA
0.1	0.9	8
TAB
gzip -c "$work/tab.tsv" > "$work/tab.tsv.gz"
printf '0.5\n0.25\n0.1\n0.95\n0.32\n0.5\nNA\n' > "$work/b.txt"
"$yame" pack -f n "$work/b.txt" "$work/in.cg" 2>/dev/null
echo s1 > "$work/names.txt"
"$yame" index -s "$work/names.txt" "$work/in.cg" 2>/dev/null

PASS=0; FAIL=0
## expected values computed here, independently of the C
exp_f() { awk -v g="$1" 'NR==FNR { if (FNR > 1) { fl[FNR-1]=$1; tp[FNR-1]=$2; rh[FNR-1]=$3 } next }
    { i = FNR; b = $1
      if (b == "NA" || rh[i] == "NA") { print "NA"; next }
      amp = tp[i] - fl[i]; if (amp < 0.05) amp = 0.05
      bp = (b - fl[i]) / amp
      if (bp <= 0) f = 0; else if (bp >= 1) f = 1
      else { rg = exp(g * log(rh[i])); f = bp / (bp + rg * (1 - bp)) }
      printf "%.5f\n", f }' "$work/tab.tsv" "$work/b.txt"; }
## the calibrated .cg as text, 5 decimals (it is float32)
got_f() { "$yame" unpack -a "$1" 2>/dev/null |
    awk '{ if ($1 == "NA" || $1 == "nan") print "NA"; else printf "%.3f\n", $1 }'; }
round3() { awk '{ if ($1 == "NA") print "NA"; else printf "%.3f\n", $1 }'; }

check() {  # label, gamma ("" = default), expected-gamma
    lbl=$1; g=$2; eg=$3
    if [ -n "$g" ]; then
        "$bin" calibrate --table "$work/tab.tsv.gz" --gamma "$g" \
            "$work/in.cg" "$work/out.cg" 2>"$work/err" || { echo "FAIL $lbl: exited non-zero"; sed 's/^/    /' "$work/err"; FAIL=$((FAIL+1)); return; }
    else
        "$bin" calibrate --table "$work/tab.tsv.gz" \
            "$work/in.cg" "$work/out.cg" 2>"$work/err" || { echo "FAIL $lbl: exited non-zero"; sed 's/^/    /' "$work/err"; FAIL=$((FAIL+1)); return; }
    fi
    exp_f "$eg" | round3 > "$work/exp"
    got_f "$work/out.cg" > "$work/got"
    if cmp -s "$work/exp" "$work/got"; then
        echo "ok   $lbl: $(tr '\n' ' ' < "$work/got")"; PASS=$((PASS+1))
    else
        echo "FAIL $lbl"; paste "$work/exp" "$work/got" | sed 's/^/    exp got  /'
        FAIL=$((FAIL+1))
    fi
}
check "default gamma 0.15" "" 0.15
check "gamma 1 (the standards' own curve)" 1 1
check "gamma 0 (linear between knots)" 0 0

## the per-sample log line: 5 of the 7 probes calibrate, 1 (case 3) clamps to 0
"$bin" calibrate --table "$work/tab.tsv.gz" "$work/in.cg" "$work/out.cg" 2>"$work/err"
if grep -qF "s1: 1 of 5 calibrated probes (20.0%) clamped to 0" "$work/err"; then
    echo "ok   logs the share clamped to 0 per sample"; PASS=$((PASS+1))
else
    echo "FAIL clamp log"; sed 's/^/    /' "$work/err"; FAIL=$((FAIL+1))
fi

## refusals: a table that does not match the beta.cg, a bad header, a bad gamma
refuse() {  # label, expected stderr fragment, then the arguments
    lbl=$1; want=$2; shift 2
    if "$bin" calibrate "$@" "$work/in.cg" "$work/x.cg" 2>"$work/err"; then
        echo "FAIL refusal $lbl: exited 0"; FAIL=$((FAIL+1)); return; fi
    if grep -qF -- "$want" "$work/err"; then echo "ok   refuses $lbl"; PASS=$((PASS+1))
    else echo "FAIL refusal $lbl: stderr lacks '$want'"; sed 's/^/    /' "$work/err"; FAIL=$((FAIL+1)); fi
}
head -4 "$work/tab.tsv" | gzip -c > "$work/short.tsv.gz"
refuse "a table of the wrong length" "lineage mismatch" --table "$work/short.tsv.gz"
printf 'a\tb\tc\n1\t2\t3\n' | gzip -c > "$work/hdr.tsv.gz"
refuse "a table without the floor/top/rho header" "header must be" --table "$work/hdr.tsv.gz"
refuse "gamma outside [0, 1]" "--gamma must be" --table "$work/tab.tsv.gz" --gamma 2
refuse "no table and no platform" "needs --platform or --table"

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
