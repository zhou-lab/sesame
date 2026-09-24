#!/bin/sh
# mLiftOver, two kinds of lift through one map.
#
# Array -> array (probe-ID prefix join): C's OWN raw betas are lifted and
# compared to R's mLiftOver of the same source vector -- so the mapping (target
# set, order, first-match choice, NA pattern) is checked exactly, independent of
# any raw/prep divergence. Values must match within the float32 .cg precision.
# Needs the R oracle.
#
# Array -> genome -> array (coordinate join): R's mLiftOver is array-only, so
# there is no oracle. The reference is built here, independently, from the same
# two store files with yame + awk -- the coord table gives each probe chrm and
# CpG_beg (0-based), cpg_nocontig.cr gives the universe (beg0/end1), the key is
# chrm_(beg0+1) on both sides, replicate probes on one CpG are averaged. Gates:
# row count == the universe; covered rows identical; M within 1 of the reference
# (its betas pass through 3-decimal text); the lift back returns every
# non-replicate covered beta within 0.0051 (M/100 against a float32 beta).
#
# Needs orderings (store or testdata/), the yame + pipeline_dump binaries, IDATs;
# each part SKIPs on its own missing inputs.
set -eu

## The R oracle. Overridable because the binary is not called the same
## thing everywhere: on the lab HPC plain `Rscript` is a different R
## without a usable sesame, and the one to use is Rscript-4.6.0.
RSCRIPT=${RSCRIPT:-Rscript}

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
dump="$root/pipeline_dump"
yame="$root/YAME/yame"
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
[ -x "$yame" ] || { echo "SKIP mLiftOver: no $yame"; exit 0; }
have_r=0; command -v "$RSCRIPT" >/dev/null 2>&1 && have_r=1

find_ord() {   # echo the first existing ordering for platform $1
    for c in "$store/$1/$1.ordering.tsv.gz" "$root/testdata/$1.ordering.tsv.gz"; do
        [ -f "$c" ] && { echo "$c"; return 0; }
    done
    ## Explicitly 0: without it the last failed `[ -f ]` is this function's
    ## status, and `so=$(find_ord ...)` then kills the script under `set -e` --
    ## which is how `make test` aborted at test-liftover with no message and
    ## skipped the twelve targets after it, instead of printing SKIP like every
    ## other driver.
    return 0
}

PASS=0; FAIL=0
run_pair() {
    sp=$1; rel=$2; tp=$3
    ## SESAME_ONLY: an ERE over "<platform> <prefix> ..."; cases that do not
    ## match are skipped so a parallel gate can run one case per job.
    if [ -n "${SESAME_ONLY:-}" ] && ! echo "$*" | grep -Eq "$SESAME_ONLY"; then return; fi
    so=$(find_ord "$sp"); to=$(find_ord "$tp")
    pfx="$idats/$rel"
    [ -n "$so" ] && [ -n "$to" ] || { echo "SKIP $sp->$tp: missing ordering"; return; }
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "SKIP $sp->$tp: no IDAT $pfx"; return; fi
    [ "$have_r" = 1 ] || { echo "SKIP $sp->$tp: no Rscript"; return; }

    "$bin" preprocess --platform "$sp" --index "$so" --output beta --prep "" \
        --out "$work/pp" "$pfx" 2>/dev/null
    "$bin" mliftover --to "$tp" --platform "$sp" --index "$so" --index-to "$to" \
        "$work/pp/beta.cg" "$work/out.cg" 2>/dev/null
    "$yame" unpack "$work/out.cg" 2>/dev/null > "$work/vals.txt"
    zcat < "$to" | tail -n +2 | cut -f1 > "$work/tids.txt"
    paste "$work/tids.txt" "$work/vals.txt" > "$work/c.txt"
    "$dump" --index "$so" --prep "" --what beta "$pfx" 2>/dev/null > "$work/src.txt"

    if "$RSCRIPT" --vanilla - "$sp" "$tp" "$work/src.txt" "$work/c.txt" <<'PY' 2>"$work/r.err"
suppressMessages(library(sesame))
a <- commandArgs(TRUE); sp<-a[1]; tp<-a[2]
sv <- read.table(a[3], colClasses=c("character","numeric"))
rl <- mLiftOver(setNames(sv$V2, sv$V1), tp, sp)
cc <- read.table(a[4], colClasses=c("character","numeric")); cL<-setNames(cc$V2,cc$V1)

## Control probes are named differently by the two orderings, and only by
## name: the store's EPIC ordering calls them ctl_<address>_<TYPE>
## (ctl_10609447_NEGATIVE), sesameData's EPIC.address$ordering calls them
## ctl_<address>. Measured 2026-09-20: 635 such probes on EPIC, and with the
## type suffix stripped the two sets are IDENTICAL -- same addresses, none of
## them a real probe. Comparing the raw strings reported 1270 set differences
## for a difference in spelling, so normalise before comparing. A real probe
## missing from either side still fails, which is what this gate is for.
ctl_norm <- function(v) sub("^(ctl_[0-9]+)_.*$", "\\1", v)
names(rl) <- ctl_norm(names(rl)); names(cL) <- ctl_norm(names(cL))
rl <- rl[!duplicated(names(rl))]; cL <- cL[!duplicated(names(cL))]

ro <- length(setdiff(names(rl),names(cL))); co <- length(setdiff(names(cL),names(rl)))
ids <- intersect(names(rl), names(cL))
namis <- sum(xor(is.na(rl[ids]), is.na(cL[ids])))
d <- abs(rl[ids]-cL[ids]); d <- d[!is.na(d)]; mx <- if(length(d)) max(d) else 0
cat(sprintf("ok   %-6s -> %-6s: %d probes, set-diff=%d NAmis=%d max|diff|=%.2e\n",
            sp, tp, length(ids), ro+co, namis, mx))
if (ro+co != 0 || namis != 0 || mx > 2e-3) { cat("FAIL: mLiftOver diverges\n"); quit(status=1) }
PY
    then PASS=$((PASS+1)); else sed 's/^/    /' "$work/r.err" | head -4; FAIL=$((FAIL+1)); fi
}

## Formats other than beta: a .cm carries no numbers, so the reductions are
## different code -- OR for a mask bit (fmt0), first-wins for a state track
## (fmt2). EPICv2 -> EPIC is many-to-one (replicate probes collapse), which is
## exactly where a reduction can be wrong without any single probe looking it.
run_cm() {
    sp=$1; tp=$2; file=$3; label=$4
    if [ -n "${SESAME_ONLY:-}" ] && ! echo "$*" | grep -Eq "$SESAME_ONLY"; then return; fi
    so=$(find_ord "$sp"); to=$(find_ord "$tp")
    src="$store/$sp/$file"
    [ -n "$so" ] && [ -n "$to" ] && [ -f "$src" ] || {
        echo "SKIP $label: need $file and both orderings in the store"; return; }

    "$bin" mliftover --platform "$sp" --to "$tp" --index "$so" --index-to "$to" \
        "$src" "$work/cm.out" 2>/dev/null || {
        echo "FAIL $label: mliftover errored"; FAIL=$((FAIL+1)); return; }

    ## rows must be the TARGET universe, and the sample count must survive
    nt=$(zcat "$to" | tail -n +2 | wc -l | tr -d ' ')
    got=$("$yame" info "$work/cm.out" | awk 'NR==2{print $4}')
    nsi=$("$yame" info "$src"          | awk 'NR==2{print $3}')
    nso=$("$yame" info "$work/cm.out"  | awk 'NR==2{print $3}')
    fmi=$("$yame" info "$src"          | awk 'NR==2{print $5}')
    fmo=$("$yame" info "$work/cm.out"  | awk 'NR==2{print $5}')
    [ "$got" = "$nt" ] || { echo "FAIL $label: $got rows, target ordering has $nt"; FAIL=$((FAIL+1)); return; }
    [ "$nsi" = "$nso" ] || { echo "FAIL $label: $nsi samples in, $nso out"; FAIL=$((FAIL+1)); return; }
    [ "$fmi" = "$fmo" ] || { echo "FAIL $label: format $fmi in, $fmo out"; FAIL=$((FAIL+1)); return; }

    ## The join our C does, redone here: probes match on the cg number before
    ## the design suffix. For every target probe with at least one source, the
    ## lifted value must be what the reduction says -- OR over the sources for
    ## a mask, the first source for a state track.
    "$yame" unpack -a "$src"          > "$work/cm.src.txt" 2>/dev/null
    "$yame" unpack -a "$work/cm.out"  > "$work/cm.out.txt" 2>/dev/null
    zcat "$so" | tail -n +2 | cut -f1 > "$work/cm.sid"
    zcat "$to" | tail -n +2 | cut -f1 > "$work/cm.tid"
    bad=$(\awk -v red="$5" '
        ## the join is defined on the probe families that carry a design
        ## suffix; control names differ between the two orderings by an
        ## annotation suffix (DIVERGENCES.md, C:liftover-ctl-name) and are
        ## deliberately out of scope here
        function key(id) { if (id !~ /^(cg|ch|rs)/) return ""; sub(/_.*$/, "", id); return id }
        FNR==NR && FILENAME==ARGV[1] { sid[FNR]=key($0); next }
        FILENAME==ARGV[2] { tid[FNR]=key($0); nt=FNR; next }
        FILENAME==ARGV[3] { sv[FNR]=$1; next }
        FILENAME==ARGV[4] { tv[FNR]=$1; next }
        END {
            for (i in sv) { k=sid[i]
                if (k == "") continue
                if (!(k in first)) first[k]=sv[i]
                if (sv[i] != "0" && sv[i] != "NA") any[k]=1
            }
            bad=0
            for (j=1; j<=nt; j++) { k=tid[j]
                if (k == "") continue                # control row, see key()
                if (!(k in first)) continue          # no source -> absent, not checked here
                want = (red=="or") ? ((k in any) ? 1 : 0) : first[k]
                got  = tv[j]
                if (red=="or") got = (got!="0" && got!="NA") ? 1 : 0
                if (got != want) bad++
            }
            print bad
        }' "$work/cm.sid" "$work/cm.tid" "$work/cm.src.txt" "$work/cm.out.txt")
    if [ "$bad" = "0" ]; then
        echo "ok   $label: $nt rows, $nso sample(s), fmt $fmo, reduction exact"
        PASS=$((PASS+1))
    else
        echo "FAIL $label: $bad target probes disagree with the reduction"
        FAIL=$((FAIL+1))
    fi
}

## The refusals. A .cg carries no probe ids, so the ONLY thing standing between
## a wrong file and a silently misaligned lift is the row count -- which makes
## these messages load-bearing rather than cosmetic.
run_refusals() {
    if [ -n "${SESAME_ONLY:-}" ] && ! echo "refusals" | grep -Eq "$SESAME_ONLY"; then return; fi
    so=$(find_ord EPICv2); to=$(find_ord EPIC)
    [ -n "$so" ] && [ -n "$to" ] || { echo "SKIP refusals: need both orderings"; return; }
    mu2cg="$root/mu2cg"
    [ -x "$mu2cg" ] || { echo "SKIP refusals: $mu2cg not built"; return; }

    ## a 5-row .cg is not indexed to any array
    printf 'Probe_ID\tS1_M\tS1_U\n' > "$work/r.mu.tsv"
    i=1
    while [ $i -le 5 ]; do
        printf 'cg%07d_BC21\t%d\t%d\n' "$i" "$((100*i))" "$((10*i))" >> "$work/r.mu.tsv"
        i=$((i+1))
    done
    "$mu2cg" "$work/r.mu.tsv" "$work/r.cg" >/dev/null 2>&1

    ## a coord table that does not match the ordering it is paired with
    printf 'CpG_chrm\tCpG_beg\nchr1\t5\n' > "$work/r.coord.tsv"
    gzip -f "$work/r.coord.tsv"

    n=0
    check() {  # label, expected stderr fragment, then the command
        lbl=$1; want=$2; shift 2
        if "$@" > "$work/r.out" 2> "$work/r.err"; then
            echo "FAIL refusal $lbl: exited 0"; FAIL=$((FAIL+1)); return
        fi
        if grep -qF -- "$want" "$work/r.err"; then n=$((n+1)); else
            echo "FAIL refusal $lbl: stderr lacks '$want'"
            sed 's/^/    /' "$work/r.err" | head -3; FAIL=$((FAIL+1))
        fi
    }
    check "row space" "not indexed to it" \
        "$bin" mliftover --platform EPICv2 --to EPIC --index "$so" --index-to "$to" \
               "$work/r.cg" "$work/r.lift.cg"
    check "coord lineage" "lineage mismatch" \
        "$bin" mliftover --platform EPICv2 --to "$gen_for_refusal" --index "$so" \
               --coords "$work/r.coord.tsv.gz" "$work/r.cg" "$work/r.lift2.cg"
    [ "$n" -eq 2 ] && { echo "ok   refusals: wrong row space and wrong coord lineage both named"; PASS=$((PASS+1)); }
}

## array -> genome -> array, against an in-test reference (no R)
run_genome() {
    sp=$1; rel=$2; gen=$3
    ## SESAME_ONLY: an ERE over "<platform> <prefix> ..."; cases that do not
    ## match are skipped so a parallel gate can run one case per job.
    if [ -n "${SESAME_ONLY:-}" ] && ! echo "$*" | grep -Eq "$SESAME_ONLY"; then return; fi
    so=$(find_ord "$sp"); pfx="$idats/$rel"
    cr="$store/$gen/cpg_nocontig.cr"; co="$store/$sp/$sp.$gen.coord.tsv.gz"
    [ -n "$so" ] && [ -f "$cr" ] && [ -f "$co" ] || {
        echo "SKIP $sp->$gen: need $gen/cpg_nocontig.cr and $sp.$gen.coord.tsv.gz in the store"; return; }
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "SKIP $sp->$gen: no IDAT $pfx"; return; fi

    "$bin" preprocess --platform "$sp" --index "$so" --output beta --prep "" \
        --out "$work/gp" "$pfx" 2>/dev/null
    "$bin" mliftover --to "$gen" --platform "$sp" --index "$so" \
        --simulated-depth 100 "$work/gp/beta.cg" "$work/wg.cg" 2>/dev/null
    "$bin" mliftover --to "$sp" --platform "$gen" --index-to "$so" \
        "$work/wg.cg" "$work/back.cg" 2>/dev/null
    ## the same lifts at 4 threads must be byte-identical: records are
    ## independent and written in input order, whatever thread did them
    "$bin" mliftover --to "$gen" --platform "$sp" --index "$so" --threads 4 \
        --simulated-depth 100 "$work/gp/beta.cg" "$work/wg4.cg" 2>/dev/null
    "$bin" mliftover --to "$sp" --platform "$gen" --index-to "$so" --threads 4 \
        "$work/wg.cg" "$work/back4.cg" 2>/dev/null
    cmp -s "$work/wg.cg" "$work/wg4.cg" && cmp -s "$work/back.cg" "$work/back4.cg" \
        && same4=identical || same4=DIFFERENT

    ## the reference, from the same two store files
    zcat < "$co" | tail -n +2 | awk -F'\t' '{print $1"_"($2+1)}' > "$work/pk.txt"
    zcat < "$so" | tail -n +2 | cut -f1 > "$work/ids.txt"
    "$yame" unpack "$work/gp/beta.cg" 2>/dev/null > "$work/b.txt"
    ## cg probes only, by ID: an rs/nv probe's beta is a genotype fraction,
    ## and a multi-mapping one can sit on a CpG row (6 rs + 78 nv on EPICv2)
    paste "$work/ids.txt" "$work/pk.txt" "$work/b.txt" | awk -F'\t' '
        $1 ~ /^cg/ && $3!="NA" && $3>=0 {s[$2]+=$3; n[$2]++}
        END{for(k in s) printf "%s\t%.6f\n", k, s[k]/n[k]}' > "$work/kv.txt"
    "$yame" unpack "$cr" 2>/dev/null | awk -F'\t' '{print $1"_"($2+1)}' > "$work/uk.txt"
    ## the small map is held in awk, the big universe is streamed -- never the
    ## reverse (a 29M-row awk array has frozen a node before)
    awk -F'\t' -v kv="$work/kv.txt" '
        BEGIN{while((getline l < kv)>0){split(l,p,"\t"); v[p[1]]=p[2]}}
        {if($0 in v){m=int(v[$0]*100+0.5); print m"\t"(100-m)} else print "0\t0"}' \
        "$work/uk.txt" > "$work/ref.txt"

    "$yame" unpack -f -1 "$work/wg.cg" 2>/dev/null > "$work/c.txt"
    nu=$(wc -l < "$work/uk.txt"); nc=$(wc -l < "$work/c.txt")
    r1=$(paste "$work/c.txt" "$work/ref.txt" | awk -F'\t' '
        {a=$1+$2; b=$3+$4; if((a>0)!=(b>0)) cov++; else if(a>0){cv++; d=$1-$3; if(d<0)d=-d; if(d>mx)mx=d}}
        END{printf "%d %d %d", cv+0, cov+0, mx+0}')
    set -- $r1; covered=$1; covmis=$2; mxd=$3

    ## the way back: M/100 against the original beta, non-replicate probes only
    awk -F'\t' '{c[$0]++} END{for(k in c) if(c[k]>1) print k}' "$work/pk.txt" > "$work/dup.txt"
    "$yame" unpack -f -1 "$work/back.cg" 2>/dev/null > "$work/bk.txt"
    r2=$(paste "$work/ids.txt" "$work/pk.txt" "$work/b.txt" "$work/bk.txt" | awk -F'\t' -v dupf="$work/dup.txt" '
        BEGIN{while((getline l < dupf)>0) dup[l]=1}
        $1 !~ /^cg/ && $4+$5>0 {bad++}                     # a non-cg probe must come back 0,0
        $1 ~ /^cg/ && !($2 in dup) && $4+$5>0 && $3!="NA" {n++; d=$4/100-$3; if(d<0)d=-d; if(d>mx)mx=d}
        END{printf "%d %.4f %d", n+0, mx+0, bad+0}')
    set -- $r2; nback=$1; mxback=$2; noncg=$3

    echo "ok   $sp -> $gen -> $sp: universe=$nu rows=$nc covered=$covered coverage-mismatch=$covmis max|dM|=$mxd; back: $nback cg probes max|dbeta|=$mxback, non-cg carrying a value=$noncg; 4 threads: $same4"
    if [ "$nc" -ne "$nu" ] || [ "$covmis" -ne 0 ] || [ "$mxd" -gt 1 ] || [ "$covered" -eq 0 ] \
       || [ "$nback" -eq 0 ] || [ "$noncg" -ne 0 ] || [ "$same4" != identical ] || awk "BEGIN{exit !($mxback > 0.0051)}"; then
        echo "FAIL: genome lift diverges from the reference"; FAIL=$((FAIL+1)); return; fi
    PASS=$((PASS+1))
}

gen_for_refusal=hg38
run_refusals

run_cm EPICv2 EPIC KYCG/CGI.20220904.cm        "EPICv2->EPIC mask (fmt0, OR)"        or
run_cm EPICv2 EPIC KYCG/ChromHMM.20220303.cm  "EPICv2->EPIC states (fmt2, first)"   first

run_pair EPICv2 EPICv2/206909630040_R03C01           EPIC
run_pair EPIC   EPIC/GSM2995280_201868590258_R01C01  EPICv2
run_genome EPICv2 EPICv2/206909630040_R03C01 hg38

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
