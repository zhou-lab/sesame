#!/bin/sh
# GCT (bisConversionControl) in qc.tsv. GCT uses the store's
# <plat>.typeI_ext.tsv.gz extension lists, which the annotation pipeline derives
# uniformly across platforms; that set is deliberately slightly broader than
# sesameData's probeInfo, so GCT is NOT bit-identical to R's bisConversionControl
# (~1e-2 on EPIC/HM450 -- a curation choice, see NUMERICS.md). We validate the
# ARITHMETIC exactly instead: recompute GCT in R from the SAME ext lists + raw
# signal and require C to match. Works on every platform carrying a typeI_ext
# asset (EPIC/EPICv2/HM450/MSA at annotation v8.1+). Needs the store, R, IDATs.
set -eu

## The R oracle. Overridable because the binary is not called the same
## thing everywhere: on the lab HPC plain `Rscript` is a different R
## without a usable sesame, and the one to use is Rscript-4.6.0.
RSCRIPT=${RSCRIPT:-Rscript}

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
## The shared store yame fetch fills. It is keyed on the browser
## hierarchy -- species/platform -- so a platform's assets sit directly
## at $store/$plat/ and a genome build's at $store/<build>/.
yhome=${YAME_DATA_HOME:-${XDG_DATA_HOME:-$HOME/.local/share}/yame}
store=$yhome
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "SKIP gct: no Rscript"; exit 0; }

PASS=0; FAIL=0
run_one() {
    plat=$1; rel=$2
    ## SESAME_ONLY: an ERE over "<platform> <prefix> ..."; cases that do not
    ## match are skipped so a parallel gate can run one case per job.
    if [ -n "${SESAME_ONLY:-}" ] && ! echo "$*" | grep -Eq "$SESAME_ONLY"; then return; fi
    ord="$store/$plat/$plat.ordering.tsv.gz"
    ext="$store/$plat/$plat.typeI_ext.tsv.gz"
    pfx="$idats/$rel"
    [ -f "$ord" ] || { echo "SKIP $plat GCT: no $ord"; return; }
    [ -f "$ext" ] || { echo "SKIP $plat GCT: no $ext (annotation v8.1+)"; return; }
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "SKIP $plat GCT: no IDAT $pfx"; return; fi

    YAME_DATA_HOME="$yhome" "$bin" preprocess --platform "$plat" --output qc \
        --out "$work/pp" "$pfx" 2>/dev/null
    cgct=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="GCT")c=i} NR==2{print $c}' "$work/pp/qc.tsv")

    if "$RSCRIPT" --vanilla - "$plat" "$pfx" "$ord" "$ext" "$cgct" <<'PY' 2>"$work/r.err"
suppressMessages(library(sesame))
a <- commandArgs(TRUE); plat<-a[1]; pfx<-a[2]; ordp<-a[3]; extp<-a[4]; cg<-as.numeric(a[5])
ids <- read.table(gzfile(ordp), header=TRUE, sep="\t", stringsAsFactors=FALSE)$Probe_ID
ext <- readLines(gzfile(extp))[-1]                      # positional, drop header
stopifnot(length(ids) == length(ext))
extC <- ids[ext == "C"]; extT <- ids[ext == "T"]
sdf <- readIDATpair(pfx, platform=plat)                 # raw, as qc computes GCT
gm <- setNames(sdf$MG, sdf$Probe_ID); gu <- setNames(sdf$UG, sdf$Probe_ID)
ref <- mean(c(gm[extC], gu[extC]), na.rm=TRUE) / mean(c(gm[extT], gu[extT]), na.rm=TRUE)
d <- abs(ref - cg)
cat(sprintf("ok   %-6s GCT: ref=%.7f C=%.7f |diff|=%.2e  (%d C / %d T probes)\n",
            plat, ref, cg, d, length(extC), length(extT)))
if (!is.finite(cg) || d > 1e-4) { cat("FAIL: GCT arithmetic diverges\n"); quit(status=1) }
PY
    then PASS=$((PASS+1)); else sed 's/^/    /' "$work/r.err" | head -3; FAIL=$((FAIL+1)); fi
}

run_one EPIC   EPIC/GSM2995280_201868590258_R01C01
run_one HM450  HM450/3999492009_R01C01
run_one EPICv2 EPICv2/206909630040_R03C01
run_one MSA    MSA/207760740030_R01C03

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
