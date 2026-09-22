#!/usr/bin/env Rscript
# Oracle for `sesame vcf`: R's formatVCF genotypes the SNP probes into a VCF. Run
# it on the same raw IDAT + SNP annotation the C side uses, and dump the parsed
# per-probe genotype for comparison (Probe_ID, GT, GS, PVF).
#
#   Rscript compare_vcf.R <idat_prefix> <platform> <snp.tsv.gz> <out.tsv>

suppressMessages({ library(sesame) })

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]; platform <- args[2]; snpf <- args[3]
outfile <- if (length(args) >= 4) args[4] else "vcf_oracle.tsv"

sdf  <- readIDATpair(prefix, platform = platform)          # raw
anno <- read.table(gzfile(snpf), header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
v <- formatVCF(sdf, anno)                                  # data.frame, INFO col carries GT/GS/PVF

info <- as.character(v$INFO)
pull <- function(tag) sub(paste0(".*", tag, "=([^;]*).*"), "\\1", info)
## The SigDF's colour channel comes from sesameData's manifest; the C side
## reads the store's ordering. Where the two disagree, the out-of-band
## allele fraction of a Type-I probe is computed from the OTHER channel and
## comes out as the exact complement -- an annotation difference, not a
## divergence in the port. Emit it so the comparison can separate the two.
cols <- as.data.frame(sdf)[, c("Probe_ID", "col")]
df <- data.frame(Probe_ID = pull("Probe_ID"),
                 GT = pull("GT"),
                 GS = as.integer(v$QUAL),
                 PVF = as.numeric(pull("PVF")),
                 stringsAsFactors = FALSE)
df$Rcol <- as.character(cols$col[match(df$Probe_ID, cols$Probe_ID)])
write.table(df, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %s: %d genotyped probes\n", outfile, nrow(df)))
