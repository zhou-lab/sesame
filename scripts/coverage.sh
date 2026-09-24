#!/bin/sh
# scripts/coverage.sh — line coverage of `make test` over sesame's own sources.
#
#   scripts/coverage.sh            measure, print the table, rewrite docs/coverage.json
#   scripts/coverage.sh --print    measure and print only (leave the badge alone)
#   scripts/coverage.sh --check    measure and fail if docs/coverage.json is stale
#                                  by more than $TOLERANCE points (default 2.0)
#   --jobs=N                       run the suite with `make -jN` (one node);
#                                  default 1
#   --sbatch (with any of the above)
#                                  run the whole thing as one SLURM job (8
#                                  cores, `make -j8`) and wait: the release
#                                  gate in ~4.5 min instead of ~35 serial.
#
# The number is gcov's own summary over src/*.c and cli/*.c, counting
# EXECUTABLE lines. It is LINE coverage, not branch coverage -- branch coverage
# is lower, and for the option parsers in cli/main.c it is the more telling
# number, so measure it before quoting it.
#
# Only sesame's own code is counted. YAME is a submodule with its own suite and
# its own badge; folding its ~20k lines in here would produce a number that
# says nothing about either.
#
# Everything happens in a scratch copy. An instrumented binary must never
# become the one in the repo: it is built -O0 and writes .gcda files beside
# itself on every run.
#
# The suite is the R-oracle golden ladder, so this needs what `make test`
# needs: an R with the sesame package (RSCRIPT=, e.g. Rscript-4.6.0 on the lab
# HPC), test IDATs at $SESAME_TEST_IDATS, and a populated $YAME_DATA_HOME.
# Targets whose inputs are missing SKIP, and a skipped target's lines are
# simply uncovered -- so measure on a machine that has all three, or the number
# understates the suite.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$here"
TOLERANCE=${TOLERANCE:-2.0}
mode=--write; sbatch=0; jobs=1
for a in "$@"; do
  case "$a" in
    --sbatch) sbatch=1 ;;
    --jobs=*) jobs=${a#--jobs=} ;;
    --print|--check|--write) mode=$a ;;
    *) echo "coverage: unknown option $a" >&2; exit 2 ;;
  esac
done

command -v gcov >/dev/null || { echo "coverage: gcov not found (it ships with gcc)" >&2; exit 1; }
[ -f YAME/libyame.a ] || {
  echo "coverage: YAME/libyame.a is missing; run make first" >&2; exit 1; }

# Not /tmp: the instrumented tree plus the suite's outputs run to hundreds of
# megabytes, and /tmp on the lab machines is a small shared partition.
scratch=${TMPDIR:-$HOME/tmp}
mkdir -p "$scratch"
work=$(mktemp -d "$scratch/sesame-cov.XXXXXX")
trap 'rm -rf "$work"' EXIT HUP INT TERM

# docs/, README.md and conda-recipe/ come too: the docs gate reads all three
# (versions, invocations, and docs/examples/*.sh), and `make test` now runs it.
# Without them the gate dies on a missing file and takes the suite with it.
cp -r src cli include tests tools docs conda-recipe Makefile README.md "$work/"
rm -f "$work"/src/*.o "$work"/src/*.gcda "$work"/src/*.gcno \
      "$work"/cli/*.o "$work"/cli/*.gcda "$work"/cli/*.gcno
# YAME is SYMLINKED, not copied: it is already built, it is not what we are
# measuring, and rebuilding htslib into the scratch tree costs minutes per run.
# Its objects carry no coverage flags, so nothing of YAME's lands in the gcov
# output even though every binary links it.
ln -s "$here/YAME" "$work/YAME"
# testdata/ (the R-exported orderings, gitignored) is what the betas/prep/
# qmask/batch drivers index by; without it they SKIP -- silently thinning both
# the number and the gate. Symlinked, like YAME: it is input, not code.
[ -d "$here/testdata" ] && ln -s "$here/testdata" "$work/testdata"

# EXTRA_CFLAGS/EXTRA_LDFLAGS are appended to every compile and link (that is
# what `make asan` uses), so the coverage flags reach src/cache.o too -- it has
# its own rule with a hard-coded -O2. -O0 comes last and wins: optimisation
# merges and elides lines, and the annotation then stops meaning anything.
cov_c="-O0 -g --coverage"
cov_l="--coverage"
( cd "$work" && make EXTRA_CFLAGS="$cov_c" EXTRA_LDFLAGS="$cov_l" \
    && make EXTRA_CFLAGS="$cov_c" EXTRA_LDFLAGS="$cov_l" \
         normexp_test cbs_test pipeline_dump mu2cg ) >"$work/build.log" 2>&1 || true
# mu2cg is in that list because test-describe needs it: built lazily by the
# target instead, it would link instrumented objects without -lgcov and fail.
[ -x "$work/sesame" ] || { tail -20 "$work/build.log" >&2; echo "coverage: build failed" >&2; exit 1; }

# Counters accumulate across processes, so clear them and let ONLY the suite
# run: a stray `sesame version` adds main.c's version path and moves the total.
rm -f "$work"/src/*.gcda "$work"/cli/*.gcda

# BOTH suites: `make test` is the R-oracle ladder, `make test-docs` runs the
# documented examples. They cover different code -- region.c and most of
# describe.c are reached only by the docs gate -- and both are ours, so a number
# from one alone understates what is actually tested.
suite=0
if [ "$sbatch" -eq 1 ]; then
  # Re-run THIS script as one SLURM job and wait. The suite is ~4.5 min of
  # 8-core work, which is past the head node's budget, and a release should
  # not depend on a login session staying up.
  #
  # It used to fan out into 41 jobs, one per target and per R-oracle case.
  # That was 234 s against 267 s for the single node (measured 2026-09-24,
  # both 87.5%), because the critical path is one target (test-prep on
  # EPICv2, ~130 s at -O0) and no amount of jobs shortens it. Thirty seconds
  # did not pay for a hand-kept shard list, a marker protocol and a dead-job
  # guard, so the fan-out is gone; `make -j` inside the one job does the same
  # work. -j16 measured no faster than -j8 for the same reason.
  log="$work/sbatch.log"
  job=$(sbatch -p "${SBATCH_PARTITION:-defq}" -A "${SBATCH_ACCOUNT:-hpcusers}" \
      -c "${SBATCH_CPUS:-8}" --mem="${SBATCH_MEM:-48G}" -t "${SBATCH_TIME:-02:00:00}" \
      -J sesame-gate -o "$log" --parsable \
      --wrap "bash -lc 'export YAME_DATA_HOME=\"${YAME_DATA_HOME:-}\" \
              RSCRIPT=\"${RSCRIPT:-Rscript}\" \
              SESAME_TEST_IDATS=\"${SESAME_TEST_IDATS:-}\"; \
              \"$here/scripts/coverage.sh\" $mode --jobs=${SBATCH_CPUS:-8}; \
              echo \"===GATE=== \$?\"'") \
      || { echo "coverage: sbatch failed" >&2; exit 1; }
  echo "coverage: job $job; waiting" >&2
  # wait on the marker, not on squeue -- a hung slurmctld returns nothing --
  # but ask sacct too, so a job that dies without printing cannot hang this
  while ! grep -q '===GATE===' "$log" 2>/dev/null; do
    sleep 10
    case $(sacct -j "$job" -n -X -o State 2>/dev/null | tr -d ' ') in
      COMPLETED|FAILED|CANCELLED*|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL)
        grep -q '===GATE===' "$log" 2>/dev/null && break
        echo "coverage: job $job ended without finishing; see $log" >&2
        cat "$log" >&2 2>/dev/null || true
        exit 1 ;;
    esac
  done
  sed '/^===GATE===/d' "$log" >&2
  rc=$(sed -n 's/^===GATE=== \([0-9]*\)$/\1/p' "$log" | tail -1)
  exit "${rc:-1}"
fi

suite=0
# -j: the targets are independent and each driver works in its own mktemp
# dir, so they parallelise on one node. Their .gcda files are shared, which
# libgcov merges under its own file lock -- measured identical (87.5%) at
# -j1, -j8, -j16 and across 41 separate jobs.
( cd "$work" && make -j"$jobs" -k test RSCRIPT="${RSCRIPT:-Rscript}" ) >"$work/suite.log" 2>&1 || suite=$?
( cd "$work" && make test-docs ) >>"$work/suite.log" 2>&1 || suite=$((suite + 100))
# the totals and every SKIP by name: in a release run a SKIP is a hole (an
# input the driver did not find), not a pass, so it belongs in the output
grep -ohE 'passed [0-9]+, failed [0-9]+' "$work/suite.log" \
  | awk '{p+=$2; f+=$4} END{printf "coverage: suite total: passed %d, failed %d\n", p, f}' >&2
grep -h 'SKIP' "$work/suite.log" | sort | uniq -c | sed 's/^/  /' >&2 || true
skipped=$(grep -c 'SKIP' "$work/suite.log" || true)
if [ "$suite" -ne 0 ]; then
  echo "coverage: the suite is RED (make test exited $suite) -- see the tail below." >&2
  # -j interleaves targets, so name the failing ones and show their output
  for t in $(sed -n 's/^make: \*\*\* \[[^]]*: \([a-z-]*\)\] Error.*/\1/p' "$work/suite.log" | sort -u); do
    echo "--- $t ---" >&2
    sed -n "/^make\[.*$t/,/passed [0-9]*, failed/p" "$work/suite.log" | tail -15 >&2
  done
  tail -5 "$work/suite.log" >&2
fi
[ "$skipped" -gt 0 ] && echo "coverage: $skipped target(s) SKIPPED; their lines count as uncovered" >&2

( cd "$work" && gcov -n -o src src/*.c && gcov -n -o cli cli/*.c ) >"$work/gcov.txt" 2>/dev/null || true

pct=$(python3 - "$work/gcov.txt" <<'PY'
import re, sys
t = open(sys.argv[1]).read()
rows = re.findall(r"File '([^']+)'\nLines executed:([0-9.]+)% of (\d+)", t)
rows = [(f, float(p), int(n)) for f, p, n in rows
        if f.startswith("src/") or f.startswith("cli/")]
if not rows:
    sys.exit("coverage: gcov produced no rows for src/ or cli/")
tot = sum(n for _, _, n in rows)
cov = sum(p / 100 * n for _, p, n in rows)
for f, p, n in sorted(rows, key=lambda r: r[1]):
    print(f"  {f:<24}{n:>7}{p:>8.1f}%", file=sys.stderr)
print(f"  {'TOTAL':<24}{tot:>7}{100*cov/tot:>8.1f}%", file=sys.stderr)
print(f"{100*cov/tot:.1f}")
PY
)

echo "coverage: ${pct}% of executable lines (gcov, line coverage: make test + make test-docs)"

badge="$here/docs/coverage.json"
red=0
if [ "$suite" -ne 0 ]; then
    failed=$(sed -n 's/^make: \*\*\* \[[^]]*: \([a-z-]*\)\] Error.*/\1/p' "$work/suite.log" | sort -u)
    for t in $failed; do
        case " ${KNOWN_FAIL:-} " in (*" $t "*) ;; (*) red=1 ;; esac
    done
    [ "$red" -eq 0 ] || echo "coverage: the gate is RED on:$(echo $failed | sed 's/^/ /')" >&2
fi
case "$mode" in
  --print) exit $red ;;
  --check)
      [ "$red" -eq 0 ] || exit 1
      [ -f "$badge" ] || { echo "coverage: $badge is missing; run scripts/coverage.sh" >&2; exit 1; }
      old=$(sed -n 's/.*"message" *: *"\([0-9.]*\)%".*/\1/p' "$badge")
      awk -v a="$old" -v b="$pct" -v t="$TOLERANCE" 'BEGIN{d=a-b; if(d<0)d=-d; exit !(d>t)}' \
          && { echo "coverage: badge says ${old}% but the suite measures ${pct}% (> ${TOLERANCE} points); run scripts/coverage.sh and commit docs/coverage.json" >&2; exit 1; }
      echo "coverage: badge (${old}%) is within ${TOLERANCE} points"
      exit 0 ;;
esac

# The badge records a measurement, so refuse to write one from a red suite: a
# number that fell because a target died is not a coverage change. KNOWN_FAIL
# is the exception, and it is a list of target names, not a switch -- set it to
# exactly the failures already documented in the org, and a NEW one still stops
# the write. As of 2026-09-19 that is:
#   KNOWN_FAIL="test-pneg test-liftover test-vcf" scripts/coverage.sh
# Re-measure that list per release; it moves with the R environment.
if [ "$suite" -ne 0 ]; then
    failed=$(sed -n 's/^make: \*\*\* \[[^]]*: \([a-z-]*\)\] Error.*/\1/p' "$work/suite.log" | sort -u)
    unexpected=""
    for t in $failed; do
        case " ${KNOWN_FAIL:-} " in (*" $t "*) ;; (*) unexpected="$unexpected $t" ;; esac
    done
    [ -z "$unexpected" ] || {
        echo "coverage: refusing to write the badge; unexpected failure(s):$unexpected" >&2
        echo "coverage: set KNOWN_FAIL to the documented ones if they are expected" >&2
        exit 1; }
    echo "coverage: suite red only on the documented failures ($(echo $failed | tr '\n' ' '))" >&2
fi

colour=$(awk -v p="$pct" 'BEGIN{
  print (p>=90)?"brightgreen":(p>=80)?"green":(p>=70)?"yellowgreen":(p>=60)?"yellow":(p>=50)?"orange":"red"}')
cat > "$badge" <<JSON
{
  "schemaVersion": 1,
  "label": "coverage",
  "message": "${pct}%",
  "color": "${colour}"
}
JSON
echo "coverage: wrote docs/coverage.json (${pct}%, ${colour})"
