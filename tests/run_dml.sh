#!/bin/sh
# Level (DML): per-probe differential methylation vs R's DML. No data-lineage
# caveat -- the betas matrix is the input -- so the OLS estimates, t/F p-values,
# effect sizes, and BH adjustment must match R to ~1e-8. The R oracle generates a
# small dataset, runs `sesame dml`, and compares (see compare_dml.R).
set -eu

## The R oracle. Overridable because the binary is not called the same
## thing everywhere: on the lab HPC plain `Rscript` is a different R
## without a usable sesame, and the one to use is Rscript-4.6.0.
RSCRIPT=${RSCRIPT:-Rscript}

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }

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
