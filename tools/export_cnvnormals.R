#!/usr/bin/env Rscript
# Export a platform's copy-number NORMAL reference from sesameData -- the same
# normal samples cnSegmentation() uses by default for EPIC/EPICv2, and a normal
# subset of the public HM450/MSA sets -- as per-probe M and U aligned to our
# ordering. Writes a TSV (Probe_ID + <name>_M/<name>_U columns); pipe it
# through `mu2cg` to get the format-3 .cg the CNV command will read.
#
#   Rscript tools/export_cnvnormals.R <platform> <ordering.tsv.gz> <out.mu.tsv.gz>

suppressMessages({ library(sesame); library(sesameData) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: export_cnvnormals.R <platform> <ordering> <out>")
platform <- args[1]; ordering <- args[2]; outfile <- args[3]

## the default normal set (mirrors sesame's cnv_normal_default)
normals <- if (platform == "EPICv2") {
    sdfs <- sesameDataGet("EPICv2.8.SigDF")
    sdfs[c("GM12878_206909630042_R08C01", "GM12878_206909630040_R03C01")]
} else if (platform == "EPIC") {
    sesameDataGet("EPIC.5.SigDF.normal")
} else if (platform == "HM450") {
    ## sesameData has no HM450 .normal set; all ten HM450.10.SigDF samples are
    ## TCGA-BLCA solid-tissue normals (barcode sample type 11A)
    sdfs <- sesameDataGet("HM450.10.SigDF")
    sdfs[grepl("^TCGA-..-....-11A-", names(sdfs))]
} else if (platform == "MSA") {
    ## as for EPICv2: the GM12878 replicates only (HCT116 and LNCaP are
    ## cancer lines)
    sdfs <- sesameDataGet("MSA.9.SigDF")
    sdfs[grepl("^GM12878_", names(sdfs))]
} else {
    stop(sprintf("no default CNV normal set for %s in sesameData", platform))
}

## our ordering's probe set/order (the .cg is positional, so align to it)
ord <- read.table(gzfile(ordering), header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, colClasses = "character")
pid <- ord$Probe_ID

signalMU <- get("signalMU", envir = asNamespace("sesame"))
cols <- list(Probe_ID = pid)
for (nm in names(normals)) {
    s <- signalMU(normals[[nm]], mask = FALSE)      # per-probe M, U (raw)
    idx <- match(pid, s$Probe_ID)                   # align to our ordering
    cols[[paste0(nm, "_M")]] <- s$M[idx]
    cols[[paste0(nm, "_U")]] <- s$U[idx]
}
df <- as.data.frame(cols, check.names = FALSE)
gz <- gzfile(outfile, "w")
write.table(df, gz, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
close(gz)
cat(sprintf("wrote %s: %d probes x %d normals (%s)\n",
            outfile, length(pid), length(normals),
            paste(names(normals), collapse = ", ")))
