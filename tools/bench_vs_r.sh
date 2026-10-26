#!/bin/bash
## Time the same job three ways: the C binary, R's openSesame measured inside a
## live session, and R measured as a user actually pays for it -- one shell
## invocation, R startup and library load included.
##
## The middle number is the honest algorithmic comparison; the third is the one
## that matters in a pipeline, where every sample is a fresh process. Reporting
## only one of them would be a half-truth in either direction.
##
##   tools/bench_vs_r.sh <IDAT-prefix> <platform> [reps] > bench.tsv
##
## Needs: ./sesame built, an R with sesame (set RSCRIPT), a populated store.
set -euo pipefail

pfx=${1:?IDAT prefix}
plat=${2:?platform}
reps=${3:-3}
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
RSCRIPT=${RSCRIPT:-Rscript}

## median of the numbers on stdin, so one slow run cannot set the result.
## One line, no ternary: BSD awk (macOS) rejects a statement split over lines.
med() {
    sort -n |
        awk '{v[NR]=$1} END{if (NR%2) print v[(NR+1)/2]; else print (v[NR/2]+v[NR/2+1])/2}'
}

work=$(mktemp -d); trap 'rm -rf "$work"' EXIT

## 1. C: whole process, since that is all there is
for i in $(seq "$reps"); do
    s=$(date +%s.%N)
    "$root/sesame" preprocess --output beta --out "$work/c$i" "$pfx" >/dev/null 2>&1
    e=$(date +%s.%N); echo "$e - $s" | bc
done | med > "$work/c.txt"

## 2+3. R: the in-session time and the full wall clock of the same invocation.
##      Printing the inner time to stderr keeps it out of the outer measurement.
for i in $(seq "$reps"); do
    s=$(date +%s.%N)
    ## ONE openSesame call per invocation. The outer clock is then exactly
    ## what a pipeline pays per sample; the inner number, printed to stderr so
    ## it stays outside that clock, is the same call without R getting up.
    "$RSCRIPT" --vanilla -e '
        a <- commandArgs(TRUE)
        suppressMessages(library(sesame))
        t <- system.time(b <- openSesame(a[1], platform = a[2]))["elapsed"]
        cat(sprintf("%.3f\n", t), file = stderr())
    ' "$pfx" "$plat" 2>"$work/inner$i" >/dev/null
    e=$(date +%s.%N); echo "$e - $s" | bc
done | med > "$work/r_total.txt"
cat "$work"/inner* | med > "$work/r_inner.txt"

c=$(cat "$work/c.txt"); ri=$(cat "$work/r_inner.txt"); rt=$(cat "$work/r_total.txt")
printf 'tool\tseconds\tspeedup\n'
printf 'sesame (C)\t%.2f\t%.1f\n' "$c" 1
printf 'R openSesame\t%.2f\t%.1f\n' "$ri" "$(echo "$ri / $c" | bc -l)"
printf 'R + startup\t%.2f\t%.1f\n' "$rt" "$(echo "$rt / $c" | bc -l)"
