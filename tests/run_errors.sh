#!/bin/sh
# Negative tests: the errors a misconfigured run hits first.
#
# The golden ladder only ever feeds the binary WELL-FORMED input, so the whole
# error surface -- every sesame__fail() in the library and every "how to fix
# it" message in cache.c -- was unexercised: src/util.c sat at 0% coverage
# because nothing in the suite ever failed. These are the branches a first-time
# user reaches (empty store, wrong file, truncated download), and a wrong
# message there costs more than a wrong number.
#
# Fully self-contained: fixtures are built here, the store is an empty temp
# directory, and nothing needs IDATs, the real store or R.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
mu2cg="$root/mu2cg"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ]   || { echo "FAIL: $bin not built"; exit 1; }
[ -x "$mu2cg" ] || { echo "SKIP errors: $mu2cg not built"; exit 0; }

## The store every case reads. Empty, and OURS: a case must not find a real
## file and pass for the wrong reason.
empty="$work/emptystore"; mkdir -p "$empty"
## ... and one holding an ordering and nothing else, to get past index
## resolution and onto the missing-companion errors.
fake="$work/fakestore";  mkdir -p "$fake/EPICv2"
export YAME_DATA_HOME="$empty"
## XDG_DATA_HOME is the fallback when YAME_DATA_HOME is unset; point it
## somewhere empty too, so a case that clears the former cannot reach ~/.
export XDG_DATA_HOME="$work/xdg"

fail=0
pass=0

## expect <label> <expected-substring-of-stderr> -- <command...>
## The command must exit non-zero, print the substring on STDERR, and write
## nothing to stdout: a half-written table is worse than none.
expect() {
    label=$1; want=$2; shift 3
    if "$@" > "$work/out" 2> "$work/err"; then
        echo "FAIL $label: exited 0, expected an error"; fail=$((fail+1)); return
    fi
    if ! grep -qF -- "$want" "$work/err"; then
        echo "FAIL $label: stderr lacks '$want'"; sed 's/^/    /' "$work/err" | head -3
        fail=$((fail+1)); return
    fi
    if [ -s "$work/out" ]; then
        echo "FAIL $label: wrote to stdout before failing"; fail=$((fail+1)); return
    fi
    printf 'ok   %-34s %s\n' "$label" "$want"
    pass=$((pass+1))
}

# --- fixtures ----------------------------------------------------------------
printf 'not an IDAT at all\n' > "$work/text.txt"
printf 'IDAT'                     > "$work/magic.idat"      # header stops here
printf 'IDAT\002\0\0\0\0\0\0\0'   > "$work/v2.idat"         # version 2
printf 'IDAT\003\0\0\0\0\0\0\0'   > "$work/v3.idat"         # v3, no nFields

printf 'Probe_ID\tM\tU\tcol\n' > "$work/ord.tsv"
i=1
while [ $i -le 5 ]; do
    printf 'cg%07d_BC21\tNA\t%d\t2\n' "$i" "$((1000+i))" >> "$work/ord.tsv"
    i=$((i+1))
done
gzip -c "$work/ord.tsv" > "$fake/EPICv2/EPICv2.ordering.tsv.gz"

printf 'Probe_ID\tS1_M\tS1_U\n' > "$work/mu.tsv"
i=1
while [ $i -le 5 ]; do
    printf 'cg%07d_BC21\t%d\t%d\n' "$i" "$((100*i))" "$((10*i))" >> "$work/mu.tsv"
    i=$((i+1))
done
"$mu2cg" "$work/mu.tsv" "$work/s1.cg" >/dev/null 2>&1
printf 'cg0000001_BC21\n' > "$work/ids.txt"

# --- idat.c: what a file that is not the IDAT you think it is looks like -----
expect "idat: not an IDAT"    "invalid IDAT magic"       -- "$bin" idat-dump --tsv "$work/text.txt"
expect "idat: truncated head" "short read on version"    -- "$bin" idat-dump --tsv "$work/magic.idat"
expect "idat: version 2"      "unsupported IDAT version" -- "$bin" idat-dump --tsv "$work/v2.idat"
expect "idat: truncated body" "short read on nFields"    -- "$bin" idat-dump --tsv "$work/v3.idat"
expect "idat: missing file"   "cannot open"              -- "$bin" idat-dump --tsv "$work/nosuch.idat"

## The IDAT header is read field by field, and every one of those reads is a
## place a truncated download stops. A file that ends mid-field must say where.
printf 'IDAT\003\0\0\0\0\0\0\0\001\0\0\0' > "$work/short_table.idat"   # 1 field, no table
printf 'IDAT\003\0\0\0\0\0\0\0\377\377\0\0' > "$work/many_fields.idat"  # 65535 fields
expect "idat: short field table" "short read in field table" \
    -- "$bin" idat-dump --tsv "$work/short_table.idat"
expect "idat: implausible nFields" "implausible nFields" \
    -- "$bin" idat-dump --tsv "$work/many_fields.idat"

# --- index.c: an ordering that is not one ------------------------------------
## The ordering defines the row space of every positional file, so a malformed
## one has to stop the run: a wrong column count or a bad address here would
## otherwise shift every probe silently.
printf 'Probe_ID\tM\tU\tcol\n' | gzip -c > "$work/ord_empty.tsv.gz"
printf 'Probe_ID\tM\tU\tcol\ncg0000001_BC21\tNA\t7\n' | gzip -c > "$work/ord_narrow.tsv.gz"
printf 'Probe_ID\tM\tU\tcol\ncg0000001_BC21\tNA\tnotanumber\t2\n' | gzip -c > "$work/ord_addr.tsv.gz"
expect "index: no data rows"  "no data rows" \
    -- "$bin" describe-probe --platform "$work/ord_empty.tsv.gz" "$work/s1.cg"
expect "index: field count"   "expected 4 or 5" \
    -- "$bin" describe-probe --platform "$work/ord_narrow.tsv.gz" "$work/s1.cg"
expect "index: bad address"   "bad address" \
    -- "$bin" describe-probe --platform "$work/ord_addr.tsv.gz" "$work/s1.cg"

# --- cache.c: an empty store must say which command fills it -----------------
## a known platform quotes the registry's ordering name and the pinned tag
expect "store: no index"      "yame fetch -y EPICv2" \
    -- "$bin" cnv --platform EPICv2 "$work/s1.cg" "$work/a.tsv" "$work/b.tsv"
expect "store: no index path" "EPICv2/EPICv2.ordering.tsv.gz" \
    -- "$bin" cnv --platform EPICv2 "$work/s1.cg" "$work/a.tsv" "$work/b.tsv"
## an unknown one has no registry row, so both fall back
expect "store: unknown plat"  "<platform>.ordering.tsv.gz" \
    -- "$bin" cnv --platform NOPE "$work/s1.cg" "$work/a.tsv" "$work/b.tsv"
expect "store: unknown tag"   "this build expects the pinned tag" \
    -- "$bin" cnv --platform NOPE "$work/s1.cg" "$work/a.tsv" "$work/b.tsv"

# --- cache.c: index present, companion absent --------------------------------
YAME_DATA_HOME="$fake"
export YAME_DATA_HOME
expect "store: no cnv panel"  "no normal panel for EPICv2" \
    -- "$bin" cnv --platform EPICv2 "$work/s1.cg" "$work/a.tsv" "$work/b.tsv"
## a genome is not a platform: it has no ordering table, so the error must not
## send the user looking for hg38.ordering.tsv.gz (it did until 2026-09-24)
expect "store: no genome"     "no genome annotation for hg38" \
    -- "$bin" mliftover --platform EPICv2 --to hg38 "$work/s1.cg" "$work/o.cg"
## a platform whose ordering IS in the store but whose companion table is not:
## the error names the file and the one fetch that brings it
expect "store: no coord"      "no EPICv2.testgenome.coord.tsv.gz in the store" \
    -- "$bin" describe-probe --platform EPICv2 --genome testgenome "$work/ids.txt"
YAME_DATA_HOME="$empty"
export YAME_DATA_HOME

# --- describe.c / index.c: the file is there, the ordering is not --------------
## a mistyped path must be reported as a path, not resolved as if it named a
## platform (which printed <store>/<the whole path>/<platform>.ordering...)
expect "describe: path ordering" "cannot open" \
    -- "$bin" describe-probe --platform "$work/nosuch.tsv.gz" "$work/s1.cg"
expect "describe: path not name"  "reads as a path" \
    -- "$bin" describe-probe --platform "$work/nosuch.tsv.gz" "$work/s1.cg"
## and no input at all is a usage error, not a crash
expect "describe: no input"   "Usage:" -- "$bin" describe-probe --platform EPICv2
## coordinate mode cannot guess the array from a list of probe IDs
expect "describe: no platform" "--genome needs --platform" \
    -- "$bin" describe-probe --genome hg38 "$work/text.txt"

echo
echo "passed $pass, failed $fail"
[ $fail -eq 0 ]
