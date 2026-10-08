#!/bin/sh
# Level-3 (Q): sesame's qualityMask must mask EXACTLY the union of the recommended
# YAME mask tracks, applied to the ordering.
#
# This is a self-consistency gate, not a differential test against R: the
# published .cm is a newer mask lineage than sesameData's KYCG object (see
# NUMERICS.md, "mask lineage"), so it cannot be bit-identical to R's qualityMask.
# What it CAN pin exactly is that sesame reads and unions the .cm correctly and
# aligns it to the ordering.
#
# Every platform sesame knows is a case: HM450, EPIC, EPICv2, MSA. A case
# FAILS (not SKIPs) when its .cm or test IDAT is missing, so no platform can
# drop out of the gate unnoticed -- that is how Q went empty on HM450/EPIC
# through v2.2.0 (zhou-lab/sesame#1) with this test green on MSA alone.
# The expected set must also be non-empty: an empty union is a lineage
# mismatch between the track names and the .cm, never a clean answer.
#
# Needs: the submodule `yame` binary and a fetched store.
set -eu

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

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }
## The submodule build, not whatever `yame` is on PATH -- a stale one in a
## shared bin would silently compare against a different mask lineage.
yame="$root/YAME/yame"
[ -x "$yame" ] || make -C "$root/YAME" >/dev/null 2>&1 || true
[ -x "$yame" ] || { echo "SKIP Q: could not build $yame"; exit 0; }

# platform  sample-prefix  <recommended track names...>
run_one() {
    plat=$1; rel=$2; shift 2
    cm=$(ls "$store/$plat"/*.cm 2>/dev/null | head -1 || true)
    pfx="$idats/$rel"
    if [ -z "$cm" ]; then
        echo "FAIL $plat Q: no .cm in $store/$plat"; FAIL=$((FAIL+1)); return; fi
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "FAIL $plat Q: no IDAT $pfx"; FAIL=$((FAIL+1)); return; fi

    # expected: ordering Probe_IDs where the yame union of the given tracks is set
    printf '%s\n' "$@" > "$work/names.txt"
    gzip -dc "$store/$plat/$plat.ordering.tsv.gz" | tail -n +2 | cut -f1 > "$work/ids.txt"
    "$yame" unpack -l "$work/names.txt" "$cm" > "$work/tab.txt" 2>/dev/null
    python3 - "$work/ids.txt" "$work/tab.txt" "$work/expected.txt" <<'PY'
import sys
ids=[l.rstrip("\n") for l in open(sys.argv[1])]
exp=set()
for i,row in enumerate(open(sys.argv[2])):
    if i>=len(ids): break
    if any(c!="0" for c in row.rstrip("\n").split("\t")): exp.add(ids[i])
open(sys.argv[3],"w").write("\n".join(sorted(exp))+"\n")
PY

    # actual: probes sesame masks (NA) under prep=Q, less those already NA
    # with no prep at all (HM450 and EPICv2 each carry two STAINING controls
    # that never get a beta)
    for q in Q ""; do
        YAME_DATA_HOME="$yhome" "$dump" --prep "$q" --what beta "$pfx" 2>/dev/null \
          | python3 -c "import sys; print('\n'.join(sorted(l.split('\t')[0] for l in sys.stdin if l.rstrip('\n').split('\t')[1]=='NA')))" \
          > "$work/na_$q.txt"
    done
    LC_ALL=C comm -23 "$work/na_Q.txt" "$work/na_.txt" > "$work/actual.txt"

    if ! grep -q . "$work/expected.txt"; then
        echo "FAIL $plat Q: none of the $# tracks is in $cm"
        FAIL=$((FAIL+1)); return
    fi
    if cmp -s "$work/expected.txt" "$work/actual.txt"; then
        echo "ok   $plat Q  $(wc -l < "$work/expected.txt" | tr -d ' ') probes = yame union of $# tracks"
        PASS=$((PASS+1))
    else
        echo "FAIL $plat Q: masked set != yame union"
        echo "  expected $(wc -l < "$work/expected.txt"|tr -d ' '), got $(wc -l < "$work/actual.txt"|tr -d ' ')"
        FAIL=$((FAIL+1))
    fi
}

PASS=0; FAIL=0
## the recommended set, one list for all four (src/mask.c REC_HUMAN)
rec="M_1baseSwitchSNPcommon_5pt M_2extBase_SNPcommon_5pt M_mapping M_nonuniq M_SNPcommon_5pt"
run_one HM450  HM450/3999492009_R01C01             $rec
run_one EPIC   EPIC/GSM2995280_201868590258_R01C01 $rec
run_one EPICv2 EPICv2/206909630040_R03C01          $rec
run_one MSA    MSA/207760740030_R01C03             $rec

## A .cm that carries none of the recommended tracks must stop Q with an
## error, not hand back an empty mask. Built from HM450's own .cm, so only
## the track NAME is wrong: M_general is real, just not recommended.
cm=$(ls "$store/HM450"/*.cm 2>/dev/null | head -1 || true)
pfx="$idats/HM450/3999492009_R01C01"
if [ -n "$cm" ] && [ -f "$pfx"_Grn.idat ]; then
    "$yame" subset "$cm" M_general > "$work/norec.cm" 2>/dev/null
    "$yame" index "$work/norec.cm" >/dev/null 2>&1
    if YAME_DATA_HOME="$yhome" "$bin" preprocess --platform HM450 \
           --mask "$work/norec.cm" --prep Q --output beta --out "$work/nr" \
           "$pfx" >/dev/null 2>"$work/nr.err"; then
        echo "FAIL HM450 Q: a .cm with no recommended track was accepted"
        FAIL=$((FAIL+1))
    elif grep -q "none of the 5 recommended Q tracks" "$work/nr.err"; then
        echo "ok   HM450 Q  refuses a .cm with no recommended track"
        PASS=$((PASS+1))
    else
        echo "FAIL HM450 Q: wrong error for a .cm with no recommended track"
        sed 's/^/    /' "$work/nr.err" | head -3
        FAIL=$((FAIL+1))
    fi
else
    echo "FAIL HM450 Q: no .cm or IDAT for the refusal case"; FAIL=$((FAIL+1))
fi

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
