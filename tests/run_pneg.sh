#!/bin/sh
# detectionPnegEcdf: the negative-control ECDF detection p-value. Compared on the
# RAW signal (prep ""), where it is channel-agnostic (pmax over each colour), so
# a pure Infinium-I G<->R channel flip between the two manifest lineages does not
# enter. Two lineage differences DO enter and are corrected for here, both
# written up in DIVERGENCES.md:
#   C:pneg-design-type   a probe the two orderings call Infinium-I on one side
#                        and Infinium-II on the other is read from different
#                        addresses, so its query is not comparable -- excluded,
#                        on a per-platform count that the harness pins.
#   C:pneg-negctl-type   sesameData's HM450 control table has no RESTORATION
#                        category and 614 NEGATIVE where Illumina's manifest has
#                        613 + 1 RESTORATION -- the extra control is dropped from
#                        R's pool so both sides run the ECDF over 613.
# Gate: bit-identical (<=1 ULP) to R's detectionPnegEcdf(sdf, return.pval=TRUE)
# over every remaining probe R reports. C additionally emits p=1 for all-NA
# control probes that R's SigDF omits (positional ordering carries them) --
# those are allowed.
#
# Needs a store with the platform ordering, the R oracle, and IDATs; skips clean.
set -eu

## The R oracle. Overridable because the binary is not called the same
## thing everywhere: on the lab HPC plain `Rscript` is a different R
## without a usable sesame, and the one to use is Rscript-4.6.0.
RSCRIPT=${RSCRIPT:-Rscript}

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
dump="$root/pipeline_dump"
## The shared store yame fetch fills. It is keyed on the browser
## hierarchy -- species/platform -- so a platform's assets sit directly
## at $store/$plat/ and a genome build's at $store/<build>/.
yhome=${YAME_DATA_HOME:-${XDG_DATA_HOME:-$HOME/.local/share}/yame}
store=$yhome
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$dump" ] || { echo "FAIL: $dump not built (make pipeline_dump)"; exit 1; }
command -v "$RSCRIPT" >/dev/null 2>&1 || { echo "SKIP pneg: no Rscript"; exit 0; }

PASS=0; FAIL=0
run_one() {
    plat=$1; rel=$2; n_design=$3
    ord="$store/$plat/$plat.ordering.tsv.gz"
    pfx="$idats/$rel"
    [ -f "$ord" ] || { echo "SKIP $plat pneg: no $ord"; return; }
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "SKIP $plat pneg: no IDAT $pfx"; return; fi

    YAME_DATA_HOME="$yhome" "$dump" --prep "" --what pval --detection pneg "$pfx" \
        2>/dev/null > "$work/c_pneg.txt"

    if "$RSCRIPT" --vanilla - "$plat" "$pfx" "$work/c_pneg.txt" "$ord" "$n_design" \
         <<'PY' 2>"$work/r.err"
suppressMessages(library(sesame))
a <- commandArgs(TRUE)
plat <- a[1]; pfx <- a[2]; cf <- a[3]; ordf <- a[4]
n_design <- as.integer(a[5])
sdf <- readIDATpair(pfx, platform=plat)

## C:pneg-negctl-type. sesameData's HM450 control table carries no
## RESTORATION category and 614 NEGATIVE; Illumina's manifest -- and the
## store's v8.1 ordering -- has 613 NEGATIVE plus one RESTORATION, at
## address 41636384, which sesameData lists as "Negative 604". R would
## run the ECDF over a pool one control larger than C's, moving ~4500
## HM450 probes by ~1.4e-3 (82/614 vs 81/613). Drop that one row so both
## sides see the same pool.
##
## Key it on the ROW NAME, not the row index: readControls() drops
## controls whose signal is NA, so attr(sdf,"controls") is shorter than
## the sesameData control table (848 vs 850 on this chip) and indices do
## not carry over -- indexing by the table's position silently removes a
## different negative control and leaves a 1/613 residual. The rownames
## are make.names(Name), so "Negative 604" arrives as "Negative.604".
## Cross-check its signal against the SigDF row for the address, so a
## rename in a future sesameData fails here rather than passing quietly.
if (plat == "HM450") {
    at <- attr(sdf, "controls")
    k  <- which(rownames(at) == "Negative.604")
    sig <- sdf[sub("^(ctl_[0-9]+).*$", "\\1", sdf$Probe_ID) == "ctl_41636384",
               c("UG","UR")]
    if (length(k) != 1 || nrow(sig) != 1 ||
        at$G[k] != sig$UG[1] || at$R[k] != sig$UR[1]) {
        cat(paste0("FAIL: cannot locate sesameData's mislabelled HM450 ",
            "RESTORATION control (address 41636384, expected row name ",
            "\"Negative.604\"). sesameData's control table moved -- ",
            "re-check C:pneg-negctl-type in DIVERGENCES.md.\n"))
        quit(status=1)
    }
    attr(sdf, "controls") <- at[-k, , drop=FALSE]
}

rp <- detectionPnegEcdf(sdf, return.pval=TRUE)
cb <- read.table(cf, colClasses=c("character","numeric"))
cP <- setNames(cb$V2, cb$V1)

## C:pneg-design-type. Probes the two orderings give a different DESIGN
## TYPE (Infinium-I one side, Infinium-II the other) are read from
## different addresses, so their ECDF queries are not comparable at all.
## Pure G<->R flips need no exclusion -- the raw pmax() over each colour
## is channel-agnostic. sesameData encodes Infinium-II as col=NA (the
## factor has levels G/R only), so a plain col.R != col.C comparison
## silently drops exactly these rows; normalise NA to "2" first.
norm <- function(v) { v <- as.character(v); v[is.na(v)] <- "2"; v }
ann <- sesameDataGet(sprintf("%s.address", plat))$ordering
ord <- read.table(gzfile(ordf), header=TRUE, sep="\t", stringsAsFactors=FALSE)
m <- merge(data.frame(id=ann$Probe_ID, col.R=norm(ann$col),
                      stringsAsFactors=FALSE),
           data.frame(id=ord$Probe_ID, col.C=norm(ord$col),
                      stringsAsFactors=FALSE), by="id")
dis    <- m[m$col.R != m$col.C, ]
design <- dis$id[xor(dis$col.R == "2", dis$col.C == "2")]
if (length(design) != n_design) {
    cat(sprintf(paste0("FAIL: %s has %d design-type disagreements between the ",
        "store ordering and sesameData, expected %d. The manifests moved; ",
        "re-check which is right and update the count with a note in ",
        "DIVERGENCES.md.\n"), plat, length(design), n_design))
    quit(status=1)
}

ids   <- intersect(names(rp), names(cP))
both  <- setdiff(ids[!is.na(rp[ids]) & !is.na(cP[ids])], design)
mx    <- if (length(both)) max(abs(rp[both]-cP[both])) else 0
## any probe R reports (non-NA) that C left NA is a real mismatch:
conly <- ids[!is.na(rp[ids]) & is.na(cP[ids])]
cat(sprintf(paste0("ok   %-8s pneg: max|diff|=%.2e over %d probes ",
                   "(C-extra NA=%d, %d design-excluded)\n"),
            plat, mx, length(both), length(conly), length(design)))
if (mx > 1e-9 || length(conly) != 0) { cat("FAIL: pneg diverges\n"); quit(status=1) }
PY
    then PASS=$((PASS+1)); else sed 's/^/    /' "$work/r.err" | head -4; FAIL=$((FAIL+1)); fi
}

## The third argument pins how many probes the two manifest lineages give a
## different design type. Measured 2026-09-22 against InfiniumAnnotation v8.1
## and sesameData 1.29.10; a change means a manifest moved -- see DIVERGENCES.md.
run_one EPICv2 EPICv2/206909630040_R03C01      0
run_one MSA    MSA/207760740030_R01C03         0
run_one EPIC   EPIC/GSM2995280_201868590258_R01C01 2
run_one HM450  HM450/3999492009_R01C01         2

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
