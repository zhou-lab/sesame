/* main.c -- sesame command line interface
 * SPDX-License-Identifier: LicenseRef-CHOP-Academic-BSD-2-Clause
 *
 * Copyright (C) 2026-present The Children's Hospital of Philadelphia
 *
 * Use of this software is available to academic and non-profit institutions
 * for research purposes under the 2-Clause BSD License; for use or transfers
 * to commercial entities, inquire with Dr. Wanding Zhou at zhouw3@chop.edu.
 * See the LICENSE file at the root of the repository for the full terms.
 */
#include "sesame.h"
#include "../src/internal.h"
#include "../src/registry.h"   /* the asset tags this build expects */
#include "yame_ui.h"           /* YAME: shared colour + usage rendering */

/* Set by the Makefile from `yame-config --version`. A build that bypasses the
 * Makefile still compiles; it just cannot name the yame it linked. */
#ifndef SESAME_YAME_VERSION
#define SESAME_YAME_VERSION "unknown"
#endif

#include <dirent.h>
#include <fcntl.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

/* The help is styled, but only when someone is looking: every colour below
 * comes from YAME's ui.c, which returns an empty string off a TTY, under
 * NO_COLOR, or on a dumb terminal. So `sesame 2>&1 | less` stays readable and
 * the bytes are what they were before any of this. Sharing the helpers with
 * yame and kycg is also what keeps the three from looking like different
 * programs. */
#define H_TITLE  yame_ui_bold()
#define H_KEY    yame_ui_cyan()
#define H_NOTE   yame_ui_dim()
#define H_OFF    yame_ui_reset()

/* --index is the old spelling of --platform <ordering>. Kept working because it
 * is in published docs and in scripts, but named once on stderr so it does not
 * quietly outlive the release it was deprecated in. "index" was always a poor
 * word for the probe ordering -- it reads like a search structure. */
static void sesame__index_deprecated(void)
{
    static int said = 0;
    if (said) return;
    said = 1;
    fprintf(stderr, "sesame: --index is deprecated; pass the ordering to "
                    "--platform instead\n");
}

/* The two facts worth stating right under the title: the backend this build is
 * coupled to, and where its store is. sesame links libyame.a from a pinned
 * submodule, so the YAME version is fixed at build time -- the store layout and
 * the annotation tags all come from it, and yame is what fills the store that
 * every ordering, mask, coord and SNP lookup resolves against. */
static void print_build_info(FILE *out)
{
    const char *env = getenv("YAME_DATA_HOME");
    char store[4096];

    sesame_store_dir(store, sizeof store);
    fprintf(out, "    %sbuilt against%s  YAME %s\n", H_NOTE, H_OFF,
            SESAME_YAME_VERSION);
    fprintf(out, "    %sstore%s          %s   %s%s%s\n", H_NOTE, H_OFF, store,
            H_NOTE, env && *env ? "(from $YAME_DATA_HOME)"
                                : "($YAME_DATA_HOME unset; yame fetch fills it)",
            H_OFF);
}

/* One row of the command index. The padding sits OUTSIDE the colour span, so
 * the reset lands right after the name rather than after a run of spaces. */
static void cmd_row(const char *name, const char *desc)
{
    int pad = 13 - (int)strlen(name);
    if (pad < 1) pad = 1;
    fprintf(stderr, "    %s%s%s%*s%s\n", H_KEY, name, H_OFF, pad, "", desc);
}

/* Top-level command index -- kept brief on purpose; each command's own -h
 * (e.g. `sesame preprocess -h`) carries the full option detail. */
static int usage(void)
{
    fprintf(stderr, "\n  %ssesame " SESAME_VERSION "%s  %s-- standalone Infinium "
            "DNA methylation toolkit%s\n\n", H_TITLE, H_OFF, H_NOTE, H_OFF);
    print_build_info(stderr);

    fprintf(stderr, "\n%sUsage%s\n    sesame <command> [options]\n",
            H_TITLE, H_OFF);

    fprintf(stderr, "\n%sPreprocessing%s\n", H_TITLE, H_OFF);
    cmd_row("preprocess", "IDATs -> betas / signal / detection p-value / QC (QCDPB)");

    fprintf(stderr, "\n%sAnalysis%s\n", H_TITLE, H_OFF);
    cmd_row("dml",       "differential methylation: per-probe linear models");
    cmd_row("cnv",       "copy number: log2 ratio vs a normal panel + CBS");
    cmd_row("vcf",       "genotype the SNP probes into a VCF (formatVCF)");
    cmd_row("region",    "a region's betas as long-form TSV, for plotting");
    cmd_row("mliftover", "lift a .cg/.cm across arrays, or to/from a genome (hg38)");
    cmd_row("impute",    "fill missing betas (matrix mean or genomic neighbours)");

    fprintf(stderr, "\n%sInspect / convert%s\n", H_TITLE, H_OFF);
    cmd_row("attach-probe", "prepend Probe_IDs to a positional .cg / .tsv");
    cmd_row("idat-dump",    "dump raw IDAT records, or a summary header");
    cmd_row("deidentify",   "remove the genetic fingerprint (SNP probes)");

    fprintf(stderr, "\n%sOther%s\n", H_TITLE, H_OFF);
    cmd_row("version", "version, the yame it links, the store");
    cmd_row("help",    "this help");

    fprintf(stderr, "\n%sData%s\n", H_TITLE, H_OFF);
    fprintf(stderr, "    %ssesame downloads nothing. %syame fetch%s%s fills the shared "
            "store\n    every tool in the suite reads.%s\n",
            H_NOTE, H_KEY, H_OFF, H_NOTE, H_OFF);

    fprintf(stderr, "\n    %ssesame <command> -h%s %sfor a command's options%s\n\n",
            H_KEY, H_OFF, H_NOTE, H_OFF);
    return 1;
}

static int usage_preprocess(void)
{
    yame_usage_head("sesame preprocess [options] <prefix|dir> [<prefix|dir> ...]");
    fputs("\n", stderr);
    yame_usage_text("Preprocess Infinium IDATs to methylation levels. Each argument is an");
    yame_usage_text("IDAT prefix (<prefix>_Grn.idat / _Red.idat, .gz ok) OR a directory,");
    yame_usage_text("searched recursively for pairs (found prefixes are sorted). Applies");
    yame_usage_text("--prep to every sample and writes one indexed YAME .cg per requested");
    yame_usage_text("output, over the whole cohort. A failed sample becomes an NA column");
    yame_usage_text("and the exit status is 1.");

    yame_usage_sec("Options:");
    yame_usage_opt("--output LIST", "comma list (default beta,intensity,pval,qc) from");
    yame_usage_cont("beta, intensity, total_intensity, pval, qc");
    yame_usage_opt("--prep STEPS", "prep code (default QCDPB): Q qualityMask,");
    yame_usage_cont("C inferChannel, D dyeBias(NL), E dyeBias(linear),");
    yame_usage_cont("P pOOBAH, B noob. Empty string = raw betas.");
    yame_usage_opt("--raw-signal", "take intensity/total from the raw signal, not the prep");
    yame_usage_opt("--detection M", "pval source: poobah (default) or pneg (negative-");
    yame_usage_cont("control ECDF). Masking always uses pOOBAH.");
    yame_usage_opt("--collapse", "average replicate probes to their cg-prefix in the");
    yame_usage_cont("beta output; axis -> beta.cg.probes");
    yame_usage_opt("--platform P", "EPIC | EPICv2 | HM450 | MSA (else from the bead");
    yame_usage_cont("count), or the path to an ordering file, which is");
    yame_usage_cont("how a custom array is named.");
    yame_usage_opt("--mask FILE", "quality/background mask .cm. Inferred for a named");
    yame_usage_cont("platform; required for a custom array using Q/P/B.");
    yame_usage_opt("--min-beads N", "mask probes with < N beads (default 0)");
    yame_usage_opt("--out DIR", "output directory (default: current directory)");
    yame_usage_opt("--tmp DIR", "scratch for the sample-major matrix (default $TMPDIR)");
    yame_usage_opt("--threads N, -t", "worker threads (default: number of CPUs)");

    yame_usage_sec("Outputs (under --out):");
    yame_usage_text("beta.cg  fmt4 betas     intensity.cg  fmt3 M/U");
    yame_usage_text("pval.cg  fmt4 det. p    total_intensity.cg  fmt4");
    yame_usage_text("qc.tsv   per-sample QC");

    yame_usage_sec("Notes:");
    yame_usage_text("Without --index the platform comes from the bead count and the");
    yame_usage_text("ordering from the shared store ($YAME_DATA_HOME), then ./.");
    return 1;
}

static int usage_dml(void)
{
    yame_usage_head("sesame dml --betas <beta.cg|matrix.tsv> [options]");
    yame_usage_text("(--formula '~a+b' --meta <s.tsv> | --design <X.tsv>)");
    fputs("\n", stderr);
    yame_usage_text("Per-probe differential methylation: OLS of each probe's betas on the");
    yame_usage_text("design, with per-coefficient t-tests, a holdout F-test per");
    yame_usage_text("categorical variable, effect sizes and BH-adjusted p-values. TSV on");
    yame_usage_text("stdout, one row per probe.");

    yame_usage_sec("Options:");
    yame_usage_opt("--betas FILE", "a preprocess beta.cg (needs --index for Probe_IDs)");
    yame_usage_cont("or a Probe_ID x sample matrix TSV");
    yame_usage_opt("--formula '~a+b'", "model over --meta columns (categoricals dummied)");
    yame_usage_opt("--meta FILE", "sample table; first column matches sample names");
    yame_usage_opt("--design FILE", "numeric design matrix (for interactions), in");
    yame_usage_cont("place of --formula/--meta");

    yame_usage_opt("--threads N, -t", "worker threads (default: number of CPUs)");
    return 1;
}

static int usage_cnv(void)
{
    yame_usage_head("sesame cnv [options] <target.cg> <segments.tsv> <bins.tsv>");
    yame_usage_text("sesame cnv --probes [options] <target.cg> <segments.tsv> <probes.tsv>");
    fputs("\n", stderr);
    yame_usage_text("Copy number: regress the target's per-probe total intensity on a");
    yame_usage_text("panel of normals (OLS), take log2(target/fitted) per probe, bin");
    yame_usage_text("along the genome (median per bin), and segment by circular binary");
    yame_usage_text("segmentation. The target may hold several samples -- one profile");
    yame_usage_text("each. Writes the CBS <segments.tsv> plus a detail file that is");
    yame_usage_text("per-bin by default, or per-probe with --probes; cinderplot reads the");
    yame_usage_text("two as separate layers (segment step + points).");

    yame_usage_sec("Arguments (all required):");
    yame_usage_opt("<target.cg>", "total-intensity .cg, from `preprocess --raw-signal");
    yame_usage_cont("--output total_intensity`; may hold several samples");
    yame_usage_opt("<segments.tsv>", "output: sample chrom start end nbin seg.mean");
    yame_usage_opt("<bins.tsv>", "output: sample chrom start end nprobes log2ratio");
    yame_usage_cont("(--probes: sample Probe_ID chrom pos log2ratio)");

    yame_usage_sec("Options:");
    yame_usage_opt("--probes", "write per-probe rows to the detail file, not per-bin");
    yame_usage_opt("--exclude LIST", "comma list of chromosomes to drop from both outputs");
    yame_usage_cont("(e.g. chrY for a female sample; unreliable on arrays)");
    yame_usage_opt("--normals FILE", "normal panel .cg -- NOT shipped with the annotation;");
    yame_usage_cont("build one with `make cnv-normals` or supply your own");
    yame_usage_opt("--platform P", "EPIC | EPICv2 | HM450 | MSA; supplies the");
    yame_usage_cont("defaults. A path to an ordering file also works,");
    yame_usage_cont("for a custom array.");
    yame_usage_opt("--index FILE", "ordering .tsv.gz (overrides --platform)");
    yame_usage_opt("--coords FILE", "probe coordinate .tsv.gz (default: the store's)");
    yame_usage_opt("--genome BUILD", "genome build (default hg38)");
    yame_usage_opt("--tilewidth BP", "genome bin width (default 50000)");
    yame_usage_opt("--min-probes N", "drop bins with < N probes (default 20)");

    yame_usage_sec("Notes:");
    yame_usage_text("coords/genome default to the store for --platform; the normal panel");
    yame_usage_text("does not, since it is not published with the annotation.");
    return 1;
}

static int usage_vcf(void)
{
    yame_usage_head("sesame vcf <prefix> [options]");
    fputs("\n", stderr);
    yame_usage_text("Genotype the SNP probes into a VCF (sesame formatVCF). Runs on the");
    yame_usage_text("raw signal -- channel inference would erase the Type-I genotype.");
    yame_usage_text("<prefix> is an IDAT pair. Sites-only VCF on stdout.");

    yame_usage_sec("Options:");
    yame_usage_opt("--platform P", "EPIC | EPICv2 | HM450 | MSA (else from the bead");
    yame_usage_cont("count), or the path to an ordering file, which is");
    yame_usage_cont("how a custom array is named.");
    yame_usage_opt("--index FILE", "ordering .tsv.gz (overrides --platform)");
    yame_usage_opt("--snp FILE", "SNP table (default: the store's");
    yame_usage_cont("<P>.<genome>.snp.tsv.gz)");
    yame_usage_opt("--genome BUILD", "genome build (default hg38)");
    yame_usage_opt("--variants", "keep only rs and channel-switching probes, dropping");
    yame_usage_cont("the non-switching Infinium-I bulk (unlikely to vary)");
    yame_usage_opt("--min-beads N", "mask probes with < N beads (default 0)");
    yame_usage_opt("--out FILE, -o", "write here instead of stdout");
    return 1;
}

static int usage_region(void)
{
    yame_usage_head("sesame region <chr:beg-end> --betas <beta.cg> [options]");
    yame_usage_text("sesame region --gene <NAME> --betas <beta.cg> [options]");
    fputs("\n", stderr);
    yame_usage_text("Extract a region's betas from a multi-sample beta.cg as long-form");
    yame_usage_text("TSV: chrom beg end Probe_ID beta sample_name -- one row per (probe");
    yame_usage_text("in region) x sample, the plot-ready feed for a locus view. Give an");
    yame_usage_text("explicit chr:beg-end (commas ok) or a gene symbol via --gene, whose");
    yame_usage_text("span is looked up in the GENCODE gene models.");

    yame_usage_sec("Options:");
    yame_usage_opt("--betas FILE", "multi-sample beta.cg (required)");
    yame_usage_opt("--gene NAME", "resolve the region from a gene symbol (e.g. ADA)");
    yame_usage_cont("rather than a chr:beg-end");
    yame_usage_opt("--gene-models F", "GENCODE genes.bed.gz for --gene (default: store)");
    yame_usage_opt("--pad BP", "extend the region by BP on each side (default 0)");
    yame_usage_opt("--platform P", "EPIC | EPICv2 | HM450 | MSA, or the path to an");
    yame_usage_cont("ordering file. Supplies Probe_IDs and coords.");
    yame_usage_opt("--index FILE", "ordering .tsv.gz (overrides --platform)");
    yame_usage_opt("--coords FILE", "probe coordinate .tsv.gz (default: the store's)");
    yame_usage_opt("--genome BUILD", "genome build (default hg38)");
    yame_usage_opt("--out FILE, -o", "write here instead of stdout");
    return 1;
}

static int usage_liftover(void)
{
    yame_usage_head("sesame mliftover --to <platform|genome> [options] <in.cx> <out.cx>");
    fputs("\n", stderr);
    yame_usage_text("Lift a .cg between row spaces (mLiftOver). Array to array joins on");
    yame_usage_text("probe-ID prefix: EPICv2/MSA carry a _suffix that EPIC/HM450/HM27");
    yame_usage_text("lack, stripped on the modern side when crossing families; each");
    yame_usage_text("target probe takes the first prefix-matched source beta (NA if");
    yame_usage_text("none). Array to GENOME (--to hg38) joins on coordinate, into the");
    yame_usage_text("genome's CpG universe (one row per CpG): replicate probes on one");
    yame_usage_text("CpG are averaged, rs/ch/unmapped probes are dropped and counted.");
    yame_usage_text("Genome to array (--platform hg38 --to MSA) takes any genome-indexed");
    yame_usage_text(".cg or .cm mask down to one row per probe; that direction is");
    yame_usage_text("lossy (CpGs with no probe are gone). Output is always positional");
    yame_usage_text("to the TARGET ordering.");

    yame_usage_sec("Options:");
    yame_usage_opt("--to P", "target: EPIC | EPICv2 | HM450 | MSA, or a genome");
    yame_usage_cont("such as hg38 whose cpg_nocontig.cr is in the store (required)");
    yame_usage_opt("--platform P", "source platform or genome, setting the direction");
    yame_usage_cont("(inferred from --index's filename if omitted)");
    yame_usage_opt("--simulated-depth N", "write betas as format 3 with M+U=N on every covered");
    yame_usage_cont("row and 0,0 elsewhere -- what a whole-genome model");
    yame_usage_cont("(methscope classify) reads. Default: keep the input format.");
    yame_usage_opt("--coords FILE", "the array side's <plat>.<genome>.coord.tsv.gz");
    yame_usage_cont("(default: the store's)");
    yame_usage_opt("--index FILE", "source ordering .tsv.gz (default: the store's)");
    yame_usage_opt("--index-to FILE", "target ordering .tsv.gz (default: the store's)");

    yame_usage_sec("Notes:");
    yame_usage_text("The empirical liftOver quality mappings and impute=TRUE are not");
    yame_usage_text("ported.");
    return 1;
}

static int usage_impute(void)
{
    yame_usage_head("sesame impute --method <mean|neighbors> [options] <in.cg> <out.cg>");
    fputs("\n", stderr);
    yame_usage_text("Fill missing (NA) beta values in a beta.cg (imputeBetas*). Output is");
    yame_usage_text("positional to the same ordering.");

    yame_usage_sec("--method mean:");
    yame_usage_text("imputeBetasMatrixByMean -- replace each NA with a mean.");
    yame_usage_opt("--axis probe", "mean of the probe across samples (default; R axis=1)");
    yame_usage_opt("--axis sample", "mean of the sample across probes (R axis=2)");

    yame_usage_sec("--method neighbors:");
    yame_usage_text("imputeBetasByGenomicNeighbors -- replace each missing probe with the");
    yame_usage_text("mean of its nearest non-missing genomic neighbours.");
    yame_usage_opt("--platform P", "EPIC | EPICv2 | HM450 | MSA (for the coord table)");
    yame_usage_opt("--coords FILE", "probe coordinate .tsv.gz (default: the store's).");
    yame_usage_cont("Matches R for mapQ>=1 probes; differs only on");
    yame_usage_cont("mapQ=0 (multi-mapping) ones.");
    yame_usage_opt("--genome BUILD", "genome build (default hg38)");
    yame_usage_opt("--max-neighbors N", "up to N nearest neighbours (default 3)");
    yame_usage_opt("--max-dist BP", "search window in bp (default 10000)");
    return 1;
}

static int usage_attach(void)
{
    yame_usage_head("sesame attach-probe [options] <file>");
    fputs("\n", stderr);
    yame_usage_text("Prepend the ordering's Probe_ID to a positional file's rows, as TSV");
    yame_usage_text("on stdout. <file> is a YAME .cg/.cm/.cx (fmt0 mask, fmt3 M/U or");
    yame_usage_text("--beta, fmt4 float) or a text .tsv[.gz], e.g. a .hg38.coord.tsv.gz.");
    yame_usage_text("The row count must match the ordering -- same platform and tag that");
    yame_usage_text("made the file.");

    yame_usage_sec("Options:");
    yame_usage_opt("--platform P", "EPIC | EPICv2 | HM450 | MSA. Omit it and the");
    yame_usage_cont("ordering comes from the filename, else from the");
    yame_usage_cont("file's row count -- each ordering has a distinct");
    yame_usage_cont("length, so a positional .cg identifies itself.");
    yame_usage_opt("--index FILE", "ordering .tsv.gz (overrides --platform)");
    yame_usage_opt("--with LIST", "also emit these ordering columns after the ID:");
    yame_usage_cont("M, U, col, mask (or all). The ordering is already");
    yame_usage_cont("open -- this saves joining it back on.");
    yame_usage_opt("--all, -a", "emit every sample column (default: first only)");
    yame_usage_opt("--beta", "read fmt3 M/U as beta");
    yame_usage_opt("--no-header", "omit the header row");
    return 1;
}

static int usage_idat(void)
{
    yame_usage_head("sesame idat-dump [options] <file.idat|file.idat.gz>");
    fputs("\n", stderr);
    yame_usage_text("Print IDAT contents. Default prints a summary header.");

    yame_usage_sec("Options:");
    yame_usage_opt("--tsv", "emit addr<TAB>mean<TAB>sd<TAB>nbeads");
    yame_usage_opt("--head N", "limit to the first N records");
    return 1;
}

static int usage_deidentify(void)
{
    yame_usage_head("sesame deidentify [options] <prefix|idat> [<out.idat>]");
    fputs("\n", stderr);
    yame_usage_text("Strip the genetic fingerprint from an IDAT: the SNP (rs) probes'");
    yame_usage_text("bead Mean intensities are zeroed (default), or reversibly scrambled");
    yame_usage_text("with --randomize --seed. Everything else is copied byte-for-byte, so");
    yame_usage_text("the file reads identically apart from the SNP beads (ports sesame's");
    yame_usage_text("deIdentify). A <prefix> processes both Grn and Red; a single");
    yame_usage_text(".idat[.gz] just that file. Output defaults to <name>_noid_*.idat,");
    yame_usage_text("or _reid_ with -r.");

    yame_usage_sec("Options:");
    yame_usage_opt("--randomize", "scramble the SNP means instead of zeroing");
    yame_usage_cont("(reversible)");
    yame_usage_opt("-r, --reidentify", "reverse a --randomize; needs the same --seed");
    yame_usage_opt("--seed N", "seed for --randomize / -r");
    yame_usage_opt("--platform P", "EPIC | EPICv2 | HM450 | MSA (else from the bead");
    yame_usage_cont("count), or the path to an ordering file, which is");
    yame_usage_cont("how a custom array is named.");
    yame_usage_opt("--index FILE", "ordering .tsv.gz (overrides --platform)");
    yame_usage_opt("--out FILE", "output path (single input only)");

    yame_usage_sec("Notes:");
    yame_usage_text("--randomize uses this tool's own PRNG, so a scrambled file");
    yame_usage_text("round-trips with `deidentify -r --seed <same N>`, not with R.");
    yame_usage_text("Zeroing is R-equivalent and irreversible.");
    return 1;
}

/* Route `sesame help <cmd>` / `sesame <cmd> -h` to the command's own help. */
static int route_usage(const char *cmd)
{
    if (!strcmp(cmd, "preprocess"))   return usage_preprocess();
    if (!strcmp(cmd, "dml"))          return usage_dml();
    if (!strcmp(cmd, "cnv"))          return usage_cnv();
    if (!strcmp(cmd, "vcf"))          return usage_vcf();
    if (!strcmp(cmd, "region"))       return usage_region();
    if (!strcmp(cmd, "mliftover"))    return usage_liftover();
    if (!strcmp(cmd, "impute"))       return usage_impute();
    if (!strcmp(cmd, "deidentify")) return usage_deidentify();
    if (!strcmp(cmd, "attach-probe")) return usage_attach();
    if (!strcmp(cmd, "idat-dump"))    return usage_idat();
    return usage();
}

/* Resolve <prefix>_Grn.idat, trying .gz as R does (R/sesame.R:331-345). */
static int resolve_idat(const char *prefix, const char *chan,
                        char *out, size_t outsz)
{
    FILE *f;
    snprintf(out, outsz, "%s_%s.idat", prefix, chan);
    if ((f = fopen(out, "rb"))) { fclose(f); return 0; }
    snprintf(out, outsz, "%s_%s.idat.gz", prefix, chan);
    if ((f = fopen(out, "rb"))) { fclose(f); return 0; }
    return -1;
}

/* --------------------------------------------------------------- betas ---
 *
 * Batch-capable: many <prefix> in one process share the opened index and the
 * loaded masks, and independent per-sample work runs across a pthread pool. */

/* Basename of a path, without modifying it (basename(3) may). */
static const char *path_basename(const char *p)
{
    const char *s = strrchr(p, '/');
    return s ? s + 1 : p;
}

/* Read one prefix's IDAT pair and build its SigDF against the shared index. If
 * expect_plat is non-NULL the sample's detected platform must match it -- this
 * guards a mixed-platform batch from silently producing a misaligned column. */
static int build_sigdf_for(const char *prefix, const sesame_index_t *ix,
                           const char *expect_plat, int min_beads,
                           sesame_sigdf_t **out, sesame_err_t *err)
{
    char gpath[4096], rpath[4096];
    sesame_idat_t *g = NULL, *r = NULL;
    sesame_sigdf_t *s = NULL;
    int rc;

    if (resolve_idat(prefix, "Grn", gpath, sizeof gpath) ||
        resolve_idat(prefix, "Red", rpath, sizeof rpath))
        return sesame__fail(err, SESAME_ERR_IO, "IDAT pair not found for %s", prefix);

    if ((rc = sesame_idat_read(gpath, &g, err)) != SESAME_OK ||
        (rc = sesame_idat_read(rpath, &r, err)) != SESAME_OK)
        goto done;

    if (expect_plat) {
        const char *p = sesame_platform_from_beads(g->n);
        if (!p || strcmp(p, expect_plat) != 0) {
            rc = sesame__fail(err, SESAME_ERR_FORMAT,
                "%s is %s (%d beads), not the batch platform %s",
                prefix, p ? p : "unknown", g->n, expect_plat);
            goto done;
        }
    }

    if (!(s = sesame_sigdf_from_idats(g, r, ix, min_beads, err))) {
        rc = err ? err->code : SESAME_ERR_IO; goto done;
    }
    *out = s; s = NULL; rc = SESAME_OK;

done:
    sesame_idat_free(g);
    sesame_idat_free(r);
    sesame_sigdf_free(s);
    return rc;
}

/* A samples x probes result store, sample-major (column j at offset j*n). Backed
 * by an unlinked temp file so the OS can page a matrix larger than RAM to disk;
 * falls back to malloc. Either backing is just a flat double * to callers. */
typedef struct {
    double *data;
    size_t  bytes;
    int     fd;         /* -1 when malloc-backed */
} betas_matrix;

static int betas_matrix_init(betas_matrix *m, int32_t n, int32_t nsamp)
{
    size_t bytes = (size_t)n * (size_t)nsamp * sizeof(double);
    const char *tmp = getenv("TMPDIR");
    char path[4096];
    int fd;

    m->data = NULL; m->bytes = bytes; m->fd = -1;
    if (bytes == 0) return -1;

    snprintf(path, sizeof path, "%s/sesame-betas-XXXXXX",
             tmp && *tmp ? tmp : "/tmp");
    if ((fd = mkstemp(path)) >= 0) {
        unlink(path);                        /* auto-clean on close/crash */
        if (ftruncate(fd, (off_t)bytes) == 0) {
            void *p = mmap(NULL, bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            if (p != MAP_FAILED) { m->data = (double *)p; m->fd = fd; return 0; }
        }
        close(fd);
    }
    m->data = (double *)malloc(bytes);       /* fallback: anonymous RAM */
    return m->data ? 0 : -1;
}

static void betas_matrix_free(betas_matrix *m)
{
    if (!m->data) return;
    if (m->fd >= 0) { munmap(m->data, m->bytes); close(m->fd); }
    else free(m->data);
    m->data = NULL;
}

/* ---------------------------------------------------------- preprocess ---
 *
 * The one pipeline command: apply --prep (default QCDPB) per sample and emit any
 * of {intensity, total_intensity, beta, pval} as YAME .cg (one indexed file per
 * output over the whole cohort) plus a qc.tsv. Signal-derived outputs reflect the
 * prep by default; --raw-signal takes intensity/total from the raw SigDF. */

typedef struct {
    char           **prefixes;
    int32_t          nsamp, n;
    const sesame_index_t *ix;
    const char      *plat, *prep;
    const uint8_t   *qmask, *bgmask, *ext;
    int32_t          qn, bgn, extn;
    int              min_beads, raw_signal, pneg_detection;
    int              want_beta, want_int, want_tot, want_pval, want_qc;
    double          *matBeta, *matM, *matU, *matTot, *matPval;   /* NULL if unwanted */
    sesame_qc_t     *qcres;
    pthread_mutex_t  lock;
    int32_t          next, nfail;
    uint32_t         status_or;
} pp_ctx;

static void fill_na(double *col, int32_t n) { int32_t k; if (col) for (k = 0; k < n; k++) col[k] = NAN; }

static void *pp_worker(void *arg)
{
    pp_ctx *c = (pp_ctx *)arg;
    int32_t n = c->n;
    for (;;) {
        int32_t j;
        sesame_sigdf_t *raw = NULL;
        sesame_err_t e;
        uint32_t st = 0;
        int rc;
        double *beta, *M, *U, *tot, *pval;

        pthread_mutex_lock(&c->lock);
        j = c->next < c->nsamp ? c->next++ : -1;
        pthread_mutex_unlock(&c->lock);
        if (j < 0) break;

        beta = c->matBeta ? c->matBeta + (size_t)j*(size_t)n : NULL;
        M    = c->matM    ? c->matM    + (size_t)j*(size_t)n : NULL;
        U    = c->matU    ? c->matU    + (size_t)j*(size_t)n : NULL;
        tot  = c->matTot  ? c->matTot  + (size_t)j*(size_t)n : NULL;
        pval = c->matPval ? c->matPval + (size_t)j*(size_t)n : NULL;

        rc = build_sigdf_for(c->prefixes[j], c->ix, c->plat, c->min_beads, &raw, &e);
        if (rc == SESAME_OK) st = raw->status;
        if (rc == SESAME_OK && c->want_qc)
            rc = sesame_qc_calc(raw, c->bgmask, c->bgn, c->ext, c->extn, &c->qcres[j], &e);
        if (rc == SESAME_OK)
            rc = sesame_pipeline(raw, c->prep, c->qmask, c->qn, c->bgmask, c->bgn,
                                 c->raw_signal, c->pneg_detection,
                                 beta, M, U, tot, pval, NULL, &e);
        sesame_sigdf_free(raw);

        if (rc != SESAME_OK) {
            fill_na(beta, n); fill_na(M, n); fill_na(U, n); fill_na(tot, n); fill_na(pval, n);
        }
        pthread_mutex_lock(&c->lock);
        if (rc != SESAME_OK) {
            c->nfail++;
            fprintf(stderr, "sesame: %s: %s\n", path_basename(c->prefixes[j]), e.msg);
        } else c->status_or |= st;
        pthread_mutex_unlock(&c->lock);
    }
    return NULL;
}

static int pp_parse_output(const char *s, pp_ctx *c)
{
    char buf[256], *t;
    snprintf(buf, sizeof buf, "%s", s);
    c->want_beta = c->want_int = c->want_tot = c->want_pval = c->want_qc = 0;
    for (t = strtok(buf, ","); t; t = strtok(NULL, ",")) {
        if      (!strcmp(t, "beta"))            c->want_beta = 1;
        else if (!strcmp(t, "intensity"))       c->want_int = 1;
        else if (!strcmp(t, "total_intensity")) c->want_tot = 1;
        else if (!strcmp(t, "pval"))            c->want_pval = 1;
        else if (!strcmp(t, "qc"))              c->want_qc = 1;
        else { fprintf(stderr, "sesame: unknown --output '%s'\n", t); return -1; }
    }
    return 0;
}

static int pp_cmp(const void *a, const void *b)
{ return strcmp(*(char *const *)a, *(char *const *)b); }

static void pp_free_prefixes(char **list, int32_t n)
{ int32_t i; if (!list) return; for (i = 0; i < n; i++) free(list[i]); free(list); }

/* Append a strdup'd copy of `p` to a growable list. -1 on OOM. */
static int pp_add(char ***list, int32_t *n, int32_t *cap, const char *p)
{
    if (*n >= *cap) {
        int32_t nc = *cap ? *cap * 2 : 16;
        char **t = (char **)realloc(*list, (size_t)nc * sizeof(char *));
        if (!t) return -1;
        *list = t; *cap = nc;
    }
    if (!((*list)[*n] = strdup(p))) return -1;
    (*n)++;
    return 0;
}

/* Recursively add every IDAT prefix under `dir` -- a file path carrying a
 * _Grn.idat[.gz] suffix -- to the list, as the path minus that suffix. */
static int pp_scan_dir(const char *dir, char ***list, int32_t *n, int32_t *cap)
{
    DIR *d = opendir(dir);
    struct dirent *de;
    if (!d) return -1;
    while ((de = readdir(d)) != NULL) {
        char full[4096]; struct stat st; size_t fl;
        if (!strcmp(de->d_name, ".") || !strcmp(de->d_name, "..")) continue;
        snprintf(full, sizeof full, "%s/%s", dir, de->d_name);
        if (stat(full, &st) != 0) continue;
        if (S_ISDIR(st.st_mode)) { pp_scan_dir(full, list, n, cap); continue; }
        fl = strlen(full);
        if (fl >= 9 && !strcmp(full + fl - 9, "_Grn.idat"))
            { full[fl - 9] = '\0'; pp_add(list, n, cap, full); }
        else if (fl >= 12 && !strcmp(full + fl - 12, "_Grn.idat.gz"))
            { full[fl - 12] = '\0'; pp_add(list, n, cap, full); }
    }
    closedir(d);
    return 0;
}

static int cmd_preprocess(int argc, char **argv)
{
    const char *idxpath = NULL, *platform = NULL, *plat = NULL, *prep = "QCDPB";
    const char *maskpath = NULL;
    const char *outlist = "beta,intensity,pval,qc", *outdir = ".", *tmpdir = NULL;
    int min_beads = 0, nthreads = 0, raw_signal = 0, pneg_detection = 0, collapse = 0, i, rc = 1;
    char **prefixes = NULL, **names = NULL;
    int32_t nsamp = 0, pcap = 0, n = 0, qn = 0, bgn = 0;
    char resolved[4096], path[4096];
    sesame_index_t *ix = NULL;
    uint8_t *qmask = NULL, *bgmask = NULL, *ext = NULL;
    int32_t extn = 0;
    sesame_qc_t *qcres = NULL;
    betas_matrix mB={NULL,0,-1}, mM={NULL,0,-1}, mU={NULL,0,-1}, mT={NULL,0,-1}, mP={NULL,0,-1};
    pp_ctx ctx; memset(&ctx, 0, sizeof ctx);
    sesame_err_t e;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--index") == 0 && i+1 < argc) {
            idxpath = argv[++i];
            sesame__index_deprecated();
        }
        else if (strcmp(argv[i], "--platform") == 0 && i+1 < argc) platform = argv[++i];
        else if (strcmp(argv[i], "--mask") == 0 && i+1 < argc) maskpath = argv[++i];
        else if (strcmp(argv[i], "--prep") == 0 && i+1 < argc) prep = argv[++i];
        else if (strcmp(argv[i], "--output") == 0 && i+1 < argc) outlist = argv[++i];
        else if (strcmp(argv[i], "--out") == 0 && i+1 < argc) outdir = argv[++i];
        else if (strcmp(argv[i], "--tmp") == 0 && i+1 < argc) tmpdir = argv[++i];
        else if (strcmp(argv[i], "--min-beads") == 0 && i+1 < argc) min_beads = (int)strtol(argv[++i],NULL,10);
        else if ((strcmp(argv[i],"--threads")==0||strcmp(argv[i],"-t")==0) && i+1<argc) nthreads = (int)strtol(argv[++i],NULL,10);
        else if (strcmp(argv[i], "--raw-signal") == 0) raw_signal = 1;
        else if (strcmp(argv[i], "--collapse") == 0) collapse = 1;
        else if (strcmp(argv[i], "--detection") == 0 && i+1 < argc) {
            const char *m = argv[++i];
            if      (!strcmp(m, "poobah")) pneg_detection = 0;
            else if (!strcmp(m, "pneg"))   pneg_detection = 1;
            else { fprintf(stderr, "sesame: --detection must be poobah or pneg\n");
                   pp_free_prefixes(prefixes, nsamp); return usage_preprocess(); }
        }
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { pp_free_prefixes(prefixes, nsamp); usage_preprocess(); return 0; }
        else if (argv[i][0] == '-' && argv[i][1] != '\0') {
            fprintf(stderr, "sesame: unknown option %s\n", argv[i]); pp_free_prefixes(prefixes, nsamp); return usage_preprocess();
        } else {
            struct stat st;
            if (stat(argv[i], &st) == 0 && S_ISDIR(st.st_mode)) {   /* a directory -> recurse for IDAT pairs */
                int32_t s0 = nsamp;
                pp_scan_dir(argv[i], &prefixes, &nsamp, &pcap);
                if (nsamp - s0 > 1) qsort(prefixes + s0, (size_t)(nsamp - s0), sizeof(char *), pp_cmp);
            } else if (pp_add(&prefixes, &nsamp, &pcap, argv[i]) != 0) {
                fprintf(stderr, "sesame: out of memory\n"); pp_free_prefixes(prefixes, nsamp); return 1;
            }

    /* --platform may name an ordering file instead of a platform (YAME's -R
     * rule: a path is a path). Then it IS the index, and there is no platform
     * name under which to look up companion assets. */
    if (sesame_platform_is_path(platform)) {
        if (!idxpath) idxpath = platform;
        platform = NULL;
    }
        }
    }
    if (nsamp == 0) {
        fprintf(stderr, "sesame: no IDAT prefixes given (a prefix, or a directory to search)\n");
        pp_free_prefixes(prefixes, nsamp); return usage_preprocess();
    }
    if (pp_parse_output(outlist, &ctx) != 0) { pp_free_prefixes(prefixes, nsamp); return 1; }
    if (tmpdir) setenv("TMPDIR", tmpdir, 1);
    mkdir(outdir, 0777);                              /* ok if it already exists */

    if (!idxpath) {
        if (!platform) {
            char gp[4096]; sesame_idat_t *g0 = NULL; int32_t beads;
            if (resolve_idat(prefixes[0], "Grn", gp, sizeof gp)) {
                fprintf(stderr, "sesame: Grn IDAT does not exist for %s\n", prefixes[0]); goto out; }
            if (sesame_idat_read(gp, &g0, &e) != SESAME_OK) { fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
            beads = g0->n; platform = sesame_platform_from_beads(beads); sesame_idat_free(g0);
            if (!platform) { fprintf(stderr, "sesame: cannot identify platform from %d beads.\n"
                "  pass --platform <P> or --index <path>.\n", beads); goto out; }
            fprintf(stderr, "sesame: detected %s (%d beads)\n", platform, beads);
        }
        if (sesame_index_locate(platform, resolved, sizeof resolved) != 0) {
            char help[1024]; sesame_index_missing_help(platform, help, sizeof help);
            fprintf(stderr, "sesame: %s\n", help); goto out;
        }
        idxpath = resolved;
    }
    if (!(ix = sesame_index_open(idxpath, &e))) { fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
    n = sesame_index_nprobes(ix); plat = platform;

    /* masks: Q needs qmask; P/B/pval/qc need bgmask */
    if (strchr(prep, 'Q')) {
        if (!plat && !maskpath) { fprintf(stderr, "sesame: Q needs a mask -- name a\n"            "  commercial platform, or pass --mask <file.cm> for a custom array\n"); goto out; }
        if (sesame_quality_mask(plat, maskpath, &qmask, &qn, &e) != SESAME_OK) { fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
    }
    if (strchr(prep,'P') || strchr(prep,'B') || ctx.want_pval || ctx.want_qc) {
        if (!plat && !maskpath) { fprintf(stderr, "sesame: P/B/pval/qc need a mask -- name a\n"            "  commercial platform, or pass --mask <file.cm> for a custom array\n"); goto out; }
        if (sesame_background_mask(plat, maskpath, &bgmask, &bgn, &e) != SESAME_OK) { fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
    }
    if (ctx.want_qc && plat) {                 /* GCT ext codes (optional) */
        char efile[256], epath[4096];
        snprintf(efile, sizeof efile, "%s.typeI_ext.tsv.gz", plat);
        if (sesame_asset_locate(plat, efile, epath, sizeof epath) == 0 &&
            sesame_load_ext_codes(epath, n, &ext, &e) == SESAME_OK) extn = n;
        else ext = NULL;                       /* absent -> GCT NaN, no error */
    }

    if ((ctx.want_beta && betas_matrix_init(&mB, n, nsamp)) ||
        (ctx.want_int  && (betas_matrix_init(&mM, n, nsamp) || betas_matrix_init(&mU, n, nsamp))) ||
        (ctx.want_tot  && betas_matrix_init(&mT, n, nsamp)) ||
        (ctx.want_pval && betas_matrix_init(&mP, n, nsamp))) {
        fprintf(stderr, "sesame: cannot allocate output matrices\n"); goto out;
    }
    if (ctx.want_qc) {
        size_t k, nd = sizeof(sesame_qc_t)/sizeof(double); double *q;
        qcres = (sesame_qc_t *)malloc((size_t)nsamp * sizeof(sesame_qc_t));
        if (!qcres) { fprintf(stderr, "sesame: oom\n"); goto out; }
        q = (double *)qcres; for (k = 0; k < (size_t)nsamp*nd; k++) q[k] = NAN;
    }

    ctx.prefixes=prefixes; ctx.nsamp=nsamp; ctx.n=n; ctx.ix=ix; ctx.plat=plat; ctx.prep=prep;
    ctx.qmask=qmask; ctx.qn=qn; ctx.bgmask=bgmask; ctx.bgn=bgn; ctx.ext=ext; ctx.extn=extn;
    ctx.min_beads=min_beads; ctx.raw_signal=raw_signal; ctx.pneg_detection=pneg_detection;
    ctx.matBeta=mB.data; ctx.matM=mM.data; ctx.matU=mU.data; ctx.matTot=mT.data; ctx.matPval=mP.data;
    ctx.qcres=qcres; ctx.next=0; ctx.nfail=0; ctx.status_or=0;
    pthread_mutex_init(&ctx.lock, NULL);
    {
        long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
        int T = nthreads > 0 ? nthreads : (ncpu > 0 ? (int)ncpu : 1), t, spawned = 0;
        pthread_t *th = NULL;
        if (T > nsamp) T = nsamp;
        if (T < 1) T = 1;
        if (T > 1) th = (pthread_t *)malloc((size_t)(T-1)*sizeof(pthread_t));
        if (th) for (t=0;t<T-1;t++) if (pthread_create(&th[spawned],NULL,pp_worker,&ctx)==0) spawned++;
        pp_worker(&ctx);
        for (t=0;t<spawned;t++) pthread_join(th[t],NULL);
        free(th);
    }
    pthread_mutex_destroy(&ctx.lock);

    names = (char **)malloc((size_t)nsamp * sizeof(char *));
    if (!names) { fprintf(stderr, "sesame: oom\n"); goto out; }
    for (i = 0; i < nsamp; i++) names[i] = (char *)path_basename(prefixes[i]);

    if (ctx.want_beta && collapse) {
        /* betasCollapseToPfx: average replicate probes to their cg-prefix, then
         * write beta.cg over the collapsed axis + a beta.cg.probes sidecar (one
         * prefix per line, positional) since the axis is no longer the ordering. */
        double *cmat = NULL; char **cpfx = NULL; int32_t m = 0; FILE *pf;
        if (sesame_betas_collapse_prefix(ix, mB.data, nsamp, &cmat, &cpfx, &m, &e)) {
            fprintf(stderr,"sesame: %s\n",e.msg); goto out; }
        snprintf(path,sizeof path,"%s/beta.cg",outdir);
        if (sesame_write_cg(path, cmat, m, nsamp, names, &e)) {
            free(cmat); for (i=0;i<m;i++) free(cpfx[i]); free(cpfx);
            fprintf(stderr,"sesame: %s\n",e.msg); goto out; }
        snprintf(path,sizeof path,"%s/beta.cg.probes",outdir);
        if ((pf = fopen(path,"w"))) {
            for (i = 0; i < m; i++) fprintf(pf, "%s\n", cpfx[i]);
            fclose(pf);
        } else fprintf(stderr,"sesame: cannot write %s\n", path);
        free(cmat); for (i=0;i<m;i++) free(cpfx[i]); free(cpfx);
    } else if (ctx.want_beta) { snprintf(path,sizeof path,"%s/beta.cg",outdir);
        if (sesame_write_cg(path, mB.data, n, nsamp, names, &e)) { fprintf(stderr,"sesame: %s\n",e.msg); goto out; } }
    if (ctx.want_int)  { snprintf(path,sizeof path,"%s/intensity.cg",outdir);
        if (sesame_write_cg_mu(path, mM.data, mU.data, n, nsamp, names, &e)) { fprintf(stderr,"sesame: %s\n",e.msg); goto out; } }
    if (ctx.want_tot)  { snprintf(path,sizeof path,"%s/total_intensity.cg",outdir);
        if (sesame_write_cg(path, mT.data, n, nsamp, names, &e)) { fprintf(stderr,"sesame: %s\n",e.msg); goto out; } }
    if (ctx.want_pval) { snprintf(path,sizeof path,"%s/pval.cg",outdir);
        if (sesame_write_cg(path, mP.data, n, nsamp, names, &e)) { fprintf(stderr,"sesame: %s\n",e.msg); goto out; } }
    if (ctx.want_qc) {
        FILE *qf; snprintf(path,sizeof path,"%s/qc.tsv",outdir);
        if (!(qf = fopen(path, "w"))) { fprintf(stderr, "sesame: cannot write %s\n", path); goto out; }
        fprintf(qf, "sample\t%s\n", sesame_qc_header());
        for (i = 0; i < nsamp; i++) {
            char row[8192];
            if (sesame_qc_format_row(&qcres[i], row, sizeof row) < 0) { fclose(qf); fprintf(stderr,"sesame: qc row too long\n"); goto out; }
            fprintf(qf, "%s\t%s\n", names[i], row);
        }
        fclose(qf);
    }
    fprintf(stderr, "sesame: preprocess wrote %s to %s/ (%d sample%s, prep=%s)\n",
            outlist, outdir, nsamp, nsamp==1?"":"s", *prep?prep:"\"\"");
    if (ctx.nfail) fprintf(stderr, "sesame: %d of %d samples failed (NA)\n", ctx.nfail, nsamp);
    rc = ctx.nfail ? 1 : 0;

out:
    betas_matrix_free(&mB); betas_matrix_free(&mM); betas_matrix_free(&mU);
    betas_matrix_free(&mT); betas_matrix_free(&mP);
    free(qcres); free(names); free(qmask); free(bgmask); free(ext); pp_free_prefixes(prefixes, nsamp);
    sesame_index_close(ix);
    return rc;
}

/* ----------------------------------------------------------------- dml ---
 *
 * Per-probe differential methylation from a betas matrix + sample metadata. The
 * heavy per-probe OLS lives in src/dml.c; here we parse the inputs, build the
 * design (a simple main-effects formula, or an explicit matrix), fit in
 * parallel, BH-adjust, and emit a TSV. */

/* Read one \n-terminated line (newline stripped) into *buf (grown as needed).
 * Returns length, or -1 at EOF. Portable substitute for getline. */
static long read_line(FILE *f, char **buf, size_t *cap)
{
    size_t len = 0;
    if (!*buf) { *cap = 65536; *buf = (char *)malloc(*cap); if (!*buf) return -1; }
    for (;;) {
        if (len + 1 >= *cap) {
            char *nb = (char *)realloc(*buf, *cap * 2); if (!nb) return -1;
            *buf = nb; *cap *= 2;
        }
        if (!fgets(*buf + len, (int)(*cap - len), f)) break;
        len += strlen(*buf + len);
        if (len && (*buf)[len-1] == '\n') { (*buf)[--len] = '\0'; return (long)len; }
        if (feof(f)) break;
    }
    return len ? (long)len : -1;
}

/* Split s in place on tabs; store up to max field pointers in tok. Returns the
 * field count (may exceed max, but only max are stored). */
static int split_tabs(char *s, char **tok, int max)
{
    int n = 0;
    for (;;) {
        char *t = strchr(s, '\t');
        if (n < max) tok[n] = s;
        n++;
        if (!t) break;
        *t = '\0'; s = t + 1;
    }
    return n;
}

static double parse_num(const char *s)
{
    char *end;
    double v;
    if (!s || !*s) return NAN;
    if (strcmp(s, "NA") == 0 || strcmp(s, "NaN") == 0) return NAN;
    v = strtod(s, &end);
    return end == s ? NAN : v;
}

/* Load a betas matrix (Probe_ID + one column per sample). Fills probe ids, an
 * nprobe*nsamp row-major matrix, and the sample names. 0 on success. */
static int load_betas(const char *path, char ***ids_out, double **mat_out,
                      int32_t *nprobe_out, int32_t *nsamp_out, char ***samp_out)
{
    FILE *f = fopen(path, "rb");
    char *line = NULL, **tok = NULL;
    size_t cap = 0, cap_rows = 0;
    int32_t nsamp = 0, np = 0, i;
    char **ids = NULL, **samp = NULL;
    double *mat = NULL;
    int maxtok;

    if (!f) return -1;
    if (read_line(f, &line, &cap) < 0) { fclose(f); free(line); return -1; }
    maxtok = (int)(cap);                            /* generous upper bound */
    tok = (char **)malloc((size_t)maxtok * sizeof(char *));
    if (!tok) { fclose(f); free(line); return -1; }
    nsamp = (int32_t)split_tabs(line, tok, maxtok) - 1;   /* minus Probe_ID col */
    if (nsamp < 1) { fclose(f); free(line); free(tok); return -1; }
    samp = (char **)malloc((size_t)nsamp * sizeof(char *));
    for (i = 0; i < nsamp; i++) samp[i] = strdup(tok[i+1]);

    while (read_line(f, &line, &cap) >= 0) {
        int nf = split_tabs(line, tok, maxtok);
        if (nf < 1 || tok[0][0] == '\0') continue;
        if ((size_t)np >= cap_rows) {
            size_t nc = cap_rows ? cap_rows * 2 : 8192;
            char **ni = (char **)realloc(ids, nc * sizeof(char *));
            double *nm = (double *)realloc(mat, nc * (size_t)nsamp * sizeof(double));
            if (!ni || !nm) { free(ni?ni:ids); free(nm?nm:mat); fclose(f);
                              free(line); free(tok); free(samp); return -1; }
            ids = ni; mat = nm; cap_rows = nc;
        }
        ids[np] = strdup(tok[0]);
        for (i = 0; i < nsamp; i++)
            mat[(size_t)np*(size_t)nsamp + (size_t)i] = (i+1 < nf) ? parse_num(tok[i+1]) : NAN;
        np++;
    }
    fclose(f); free(line); free(tok);
    *ids_out = ids; *mat_out = mat; *nprobe_out = np; *nsamp_out = nsamp;
    *samp_out = samp;
    return 0;
}

/* Benjamini-Hochberg adjust p[0..n-1] (NaN preserved) into out; scratch holds n
 * pairs. */
typedef struct { double p; int32_t i; } bh_pair;
static int bh_cmp(const void *a, const void *b)
{
    double x = ((const bh_pair *)a)->p, y = ((const bh_pair *)b)->p;
    return x < y ? -1 : x > y ? 1 : 0;
}
static void bh_adjust(const double *p, int32_t n, double *out, bh_pair *pr)
{
    int32_t cnt = 0, r;
    double prev = INFINITY;
    for (r = 0; r < n; r++) { out[r] = NAN; if (!isnan(p[r])) { pr[cnt].p = p[r]; pr[cnt].i = r; cnt++; } }
    if (cnt == 0) return;
    qsort(pr, (size_t)cnt, sizeof(bh_pair), bh_cmp);
    for (r = cnt; r >= 1; r--) {                     /* from largest p down: cummin */
        double a = pr[r-1].p * (double)cnt / (double)r;
        if (a < prev) prev = a;
        out[pr[r-1].i] = prev > 1.0 ? 1.0 : prev;
    }
}

/* The built design plus the labels needed to write the header. */
typedef struct {
    double  *X;                     /* nsamp * p, row-major */
    int32_t  p, nvar;
    char   **coef;                  /* p coefficient names */
    char   **vname;                 /* nvar categorical-variable names */
    int32_t *var_off, *var_col;     /* grouping for F-tests / effect size */
} design_t;

static void design_free(design_t *d)
{
    int32_t j;
    if (d->coef) { for (j = 0; j < d->p; j++) free(d->coef[j]); free(d->coef); }
    if (d->vname) { for (j = 0; j < d->nvar; j++) free(d->vname[j]); free(d->vname); }
    free(d->X); free(d->var_off); free(d->var_col);
    memset(d, 0, sizeof *d);
}

static char *join2(const char *a, const char *b)
{
    size_t la = strlen(a), lb = strlen(b);
    char *s = (char *)malloc(la + lb + 1);
    if (s) { memcpy(s, a, la); memcpy(s + la, b, lb + 1); }
    return s;
}

static int cmp_str(const void *a, const void *b)
{ return strcmp(*(const char *const *)a, *(const char *const *)b); }

/* One static partition of probes for a fitting thread. */
typedef struct {
    const sesame_dml_design_t *d;
    const double *betas;
    int32_t lo, hi;
    double *est, *pval, *fpval, *eff;
    int32_t *nobs;
} dml_part_t;

static void *dml_worker(void *arg)
{
    dml_part_t *P = (dml_part_t *)arg;
    sesame_dml_work_t *w = sesame_dml_work_new(P->d->m, P->d->p);
    int32_t i, nv = P->d->nvar ? P->d->nvar : 1, p = P->d->p, m = P->d->m;
    if (!w) return NULL;
    for (i = P->lo; i < P->hi; i++)
        P->nobs[i] = sesame_dml_fit(P->d, w, P->betas + (size_t)i*(size_t)m,
            P->est + (size_t)i*(size_t)p, P->pval + (size_t)i*(size_t)p,
            P->fpval + (size_t)i*(size_t)nv, P->eff + (size_t)i*(size_t)nv);
    sesame_dml_work_free(w);
    return NULL;
}

static int cmd_dml(int argc, char **argv)
{
    const char *betapath = NULL, *metapath = NULL, *formula = NULL, *designpath = NULL;
    const char *idxpath = NULL;
    int nthreads = 0, i, j, rc = 1;
    char **ids = NULL, **samp = NULL;
    double *betas = NULL;
    int32_t nprobe = 0, nsamp = 0;
    design_t D; memset(&D, 0, sizeof D);
    double *est = NULL, *pval = NULL, *fpval = NULL, *eff = NULL, *padj = NULL, *fpadj = NULL;
    int32_t *nobs = NULL;
    sesame_err_t e;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--betas") == 0 && i+1 < argc) betapath = argv[++i];
        else if (strcmp(argv[i], "--meta") == 0 && i+1 < argc) metapath = argv[++i];
        else if (strcmp(argv[i], "--formula") == 0 && i+1 < argc) formula = argv[++i];
        else if (strcmp(argv[i], "--design") == 0 && i+1 < argc) designpath = argv[++i];
        else if (strcmp(argv[i], "--index") == 0 && i+1 < argc) {
            idxpath = argv[++i];
            sesame__index_deprecated();
        }
        else if ((strcmp(argv[i], "--threads")==0||strcmp(argv[i],"-t")==0) && i+1<argc)
            nthreads = (int)strtol(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage_dml(); return 0; }
        else { fprintf(stderr, "sesame: unknown/incomplete option %s\n", argv[i]); return usage_dml(); }
    }
    if (!betapath || (!formula && !designpath)) {
        fprintf(stderr, "sesame: dml needs --betas and either --formula (with --meta) or --design\n");
        return usage_dml();
    }
    { size_t bl = strlen(betapath);
      if (bl >= 3 && strcmp(betapath + bl - 3, ".cg") == 0) {   /* YAME beta.cg */
        double *smat = NULL; int32_t np = 0, ns = 0, s, p; sesame_index_t *ix;
        if (!idxpath) { fprintf(stderr, "sesame: dml on a .cg needs --index (for probe IDs)\n"); return 1; }
        if (sesame_read_cg(betapath, &smat, &np, &ns, &samp, &e)) { fprintf(stderr, "sesame: %s\n", e.msg); return 1; }
        if (!(ix = sesame_index_open(idxpath, &e))) { fprintf(stderr, "sesame: %s\n", e.msg); free(smat); return 1; }
        if (sesame_index_nprobes(ix) != np) {
            fprintf(stderr, "sesame: %s has %d probes, index has %d\n", betapath, np, sesame_index_nprobes(ix));
            sesame_index_close(ix); free(smat); return 1; }
        nprobe = np; nsamp = ns;
        ids = (char **)malloc((size_t)np * sizeof(char *));
        betas = (double *)malloc((size_t)np * (size_t)ns * sizeof(double));
        for (p = 0; p < np; p++) ids[p] = strdup(sesame_index_probe_id(ix, p));
        for (s = 0; s < ns; s++) for (p = 0; p < np; p++)   /* sample-major -> probe-major */
            betas[(size_t)p*(size_t)ns + (size_t)s] = smat[(size_t)s*(size_t)np + (size_t)p];
        free(smat); sesame_index_close(ix);
      } else if (load_betas(betapath, &ids, &betas, &nprobe, &nsamp, &samp) != 0) {
        fprintf(stderr, "sesame: cannot read betas matrix %s\n", betapath); return 1;
      }
    }

    /* ---- build the design (samples aligned to the betas columns) ---- */
    {
        FILE *mf = NULL; char *line = NULL, **tok = NULL; size_t cap = 0;
        char **mid = NULL, ***mcell = NULL, **mhdr = NULL;
        int32_t mrows = 0, mcols = 0, *smap = NULL, s;
        const char *srcpath = designpath ? designpath : metapath;
        int maxtok, ok = 1;

        if (!srcpath) { fprintf(stderr, "sesame: --formula needs --meta\n"); goto dclean; }
        mf = fopen(srcpath, "rb");
        if (!mf) { fprintf(stderr, "sesame: cannot open %s\n", srcpath); goto dclean; }
        if (read_line(mf, &line, &cap) < 0) { fprintf(stderr,"sesame: empty %s\n",srcpath); goto dclean; }
        maxtok = (int)cap;
        tok = (char **)malloc((size_t)maxtok * sizeof *tok);
        mcols = (int32_t)split_tabs(line, tok, maxtok);          /* col0 = sample id */
        mhdr = (char **)malloc((size_t)mcols * sizeof *mhdr);
        for (j = 0; j < mcols; j++) mhdr[j] = strdup(tok[j]);
        { size_t crow = 0;
          while (read_line(mf, &line, &cap) >= 0) {
            int nf = split_tabs(line, tok, maxtok);
            if (nf < 1 || tok[0][0]=='\0') continue;
            if ((size_t)mrows >= crow) {
                size_t nc = crow ? crow*2 : 256;
                mid = realloc(mid, nc*sizeof *mid); mcell = realloc(mcell, nc*sizeof *mcell);
                crow = nc;
            }
            mid[mrows] = strdup(tok[0]);
            mcell[mrows] = (char **)malloc((size_t)mcols * sizeof(char *));
            for (j = 0; j < mcols; j++) mcell[mrows][j] = strdup(j < nf ? tok[j] : "");
            mrows++;
          } }
        fclose(mf); mf = NULL;

        /* map each betas sample to its metadata row */
        smap = (int32_t *)malloc((size_t)nsamp * sizeof(int32_t));
        for (s = 0; s < nsamp; s++) {
            int32_t r, found = -1;
            for (r = 0; r < mrows; r++) if (strcmp(mid[r], samp[s]) == 0) { found = r; break; }
            if (found < 0) { fprintf(stderr, "sesame: sample %s not in %s\n", samp[s], srcpath); ok = 0; break; }
            smap[s] = found;
        }
        if (!ok) goto dclean;

        if (designpath) {                          /* explicit numeric design */
            D.p = mcols - 1; D.nvar = 0;
            D.X = (double *)malloc((size_t)nsamp * (size_t)D.p * sizeof(double));
            D.coef = (char **)malloc((size_t)D.p * sizeof(char *));
            for (j = 0; j < D.p; j++) D.coef[j] = strdup(mhdr[j+1]);
            for (s = 0; s < nsamp; s++)
                for (j = 0; j < D.p; j++)
                    D.X[(size_t)s*(size_t)D.p + (size_t)j] = parse_num(mcell[smap[s]][j+1]);
        } else {                                   /* main-effects formula */
            char fbuf[1024]; char *terms[64]; int nterm = 0, ok2 = 1;
            const char *fp = formula;
            while (*fp == ' ' || *fp == '~') fp++;
            if (strpbrk(fp, "*:^(|")) {
                fprintf(stderr, "sesame: --formula supports main effects only "
                    "(no * : ^ () ). Use --design for interactions/splines.\n");
                goto dclean;
            }
            snprintf(fbuf, sizeof fbuf, "%s", fp);
            { char *t = strtok(fbuf, "+");
              while (t && nterm < 64) {
                while (*t==' ') t++;
                { char *ep = t + strlen(t); while (ep>t && ep[-1]==' ') *--ep = '\0'; }
                if (*t && strcmp(t,"1")!=0) terms[nterm++] = t;
                t = strtok(NULL, "+");
              } }
            /* pass 1: classify each term, collect levels; count columns */
            {
                int *is_cat = (int *)calloc((size_t)nterm, sizeof(int));
                int *mcol   = (int *)malloc((size_t)nterm * sizeof(int));
                char ***levs = (char ***)calloc((size_t)nterm, sizeof(char **));
                int *nlev = (int *)calloc((size_t)nterm, sizeof(int));
                int32_t p = 1, nvar = 0, k;
                for (k = 0; k < nterm; k++) {
                    int c = -1; int cat = 0;
                    for (j = 1; j < mcols; j++) if (strcmp(mhdr[j], terms[k])==0) { c = j; break; }
                    if (c < 0) { fprintf(stderr,"sesame: no metadata column '%s'\n",terms[k]); ok2=0; break; }
                    mcol[k] = c;
                    for (s = 0; s < nsamp; s++) {
                        const char *v = mcell[smap[s]][c];
                        if (*v && isnan(parse_num(v))) { cat = 1; break; }
                    }
                    is_cat[k] = cat;
                    if (cat) {                     /* collect sorted unique levels */
                        char **lv = (char **)malloc((size_t)nsamp * sizeof(char *));
                        int nl = 0, u;
                        for (s = 0; s < nsamp; s++) {
                            const char *v = mcell[smap[s]][c]; int seen = 0;
                            for (u = 0; u < nl; u++) if (strcmp(lv[u],v)==0){seen=1;break;}
                            if (!seen) lv[nl++] = (char *)v;
                        }
                        qsort(lv, (size_t)nl, sizeof(char *), cmp_str);
                        levs[k] = lv; nlev[k] = nl;
                        p += nl - 1; nvar++;
                    } else p += 1;
                }
                if (!ok2) { free(is_cat);free(mcol);free(nlev);
                            for(k=0;k<nterm;k++) free(levs[k]); free(levs); goto dclean; }
                /* pass 2: fill X, coef names, variable groupings */
                D.p = p; D.nvar = nvar;
                D.X = (double *)calloc((size_t)nsamp * (size_t)p, sizeof(double));
                D.coef = (char **)malloc((size_t)p * sizeof(char *));
                D.vname = nvar ? (char **)malloc((size_t)nvar*sizeof(char*)) : NULL;
                D.var_off = (int32_t *)malloc((size_t)(nvar+1)*sizeof(int32_t));
                D.var_col = (int32_t *)malloc((size_t)(p)*sizeof(int32_t));
                { int32_t col = 0, vv = 0, vc = 0, lvl;
                  for (s = 0; s < nsamp; s++) D.X[(size_t)s*(size_t)p + 0] = 1.0;   /* intercept */
                  D.coef[0] = strdup("(Intercept)");
                  D.var_off[0] = 0;
                  col = 1;
                  for (k = 0; k < nterm; k++) {
                    int c = mcol[k];
                    if (!is_cat[k]) {
                        for (s = 0; s < nsamp; s++)
                            D.X[(size_t)s*(size_t)p + (size_t)col] = parse_num(mcell[smap[s]][c]);
                        D.coef[col] = strdup(terms[k]);
                        col++;
                    } else {
                        D.vname[vv] = strdup(terms[k]);
                        for (lvl = 1; lvl < nlev[k]; lvl++) {
                            for (s = 0; s < nsamp; s++) {
                                const char *v = mcell[smap[s]][c];
                                D.X[(size_t)s*(size_t)p + (size_t)col] = (strcmp(v, levs[k][lvl])==0) ? 1.0 : 0.0;
                            }
                            D.coef[col] = join2(terms[k], levs[k][lvl]);
                            D.var_col[vc++] = col;
                            col++;
                        }
                        vv++; D.var_off[vv] = vc;
                    }
                  }
                }
                free(is_cat); free(mcol); free(nlev);
                for (k = 0; k < nterm; k++) free(levs[k]);
                free(levs);
            }
            if (!ok2) goto dclean;
        }
    dclean:
        if (mf) fclose(mf);
        free(line); free(tok);
        if (mid) { for (j=0;j<mrows;j++) free(mid[j]); free(mid); }
        if (mcell) { int32_t r; for (r=0;r<mrows;r++){ if(mcell[r]){int q;for(q=0;q<mcols;q++)free(mcell[r][q]);free(mcell[r]);} } free(mcell); }
        if (mhdr) { for (j=0;j<mcols;j++) free(mhdr[j]); free(mhdr); }
        free(smap);
        if (!D.X) goto out;
    }

    /* ---- fit every probe in parallel ---- */
    est   = (double *)malloc((size_t)nprobe * (size_t)D.p * sizeof(double));
    pval  = (double *)malloc((size_t)nprobe * (size_t)D.p * sizeof(double));
    fpval = (double *)malloc((size_t)nprobe * (size_t)(D.nvar?D.nvar:1) * sizeof(double));
    eff   = (double *)malloc((size_t)nprobe * (size_t)(D.nvar?D.nvar:1) * sizeof(double));
    nobs  = (int32_t *)malloc((size_t)nprobe * sizeof(int32_t));
    if (!est || !pval || !fpval || !eff || !nobs) { fprintf(stderr,"sesame: oom\n"); goto out; }

    {
        sesame_dml_design_t dz = { nsamp, D.p, D.nvar, D.X, D.var_off, D.var_col };
        long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
        int T = nthreads > 0 ? nthreads : (ncpu > 0 ? (int)ncpu : 1);
        dml_part_t *parts; pthread_t *th; int t;
        if (T > nprobe) T = nprobe;
        if (T < 1) T = 1;
        parts = (dml_part_t *)malloc((size_t)T * sizeof *parts);
        th = (pthread_t *)malloc((size_t)T * sizeof *th);
        if (!parts || !th) { free(parts); free(th); fprintf(stderr,"sesame: oom\n"); goto out; }
        for (t = 0; t < T; t++) {
            parts[t].d = &dz; parts[t].betas = betas;
            parts[t].lo = (int32_t)((long)nprobe*t/T);
            parts[t].hi = (int32_t)((long)nprobe*(t+1)/T);
            parts[t].est=est; parts[t].pval=pval; parts[t].fpval=fpval;
            parts[t].eff=eff; parts[t].nobs=nobs;
        }
        for (t = 1; t < T; t++) pthread_create(&th[t], NULL, dml_worker, &parts[t]);
        dml_worker(&parts[0]);
        for (t = 1; t < T; t++) pthread_join(th[t], NULL);
        free(parts); free(th);
    }

    /* ---- BH-adjust each p-value column across probes ---- */
    padj  = (double *)malloc((size_t)nprobe * (size_t)D.p * sizeof(double));
    fpadj = (double *)malloc((size_t)nprobe * (size_t)(D.nvar?D.nvar:1) * sizeof(double));
    {
        bh_pair *pr = (bh_pair *)malloc((size_t)nprobe * sizeof(bh_pair));
        double *colb = (double *)malloc((size_t)nprobe * sizeof(double));
        double *cola = (double *)malloc((size_t)nprobe * sizeof(double));
        int32_t k;
        if (pr && colb && cola) {
            for (j = 0; j < D.p; j++) {
                for (k = 0; k < nprobe; k++) colb[k] = pval[(size_t)k*(size_t)D.p+(size_t)j];
                bh_adjust(colb, nprobe, cola, pr);
                for (k = 0; k < nprobe; k++) padj[(size_t)k*(size_t)D.p+(size_t)j] = cola[k];
            }
            for (j = 0; j < D.nvar; j++) {
                for (k = 0; k < nprobe; k++) colb[k] = fpval[(size_t)k*(size_t)D.nvar+(size_t)j];
                bh_adjust(colb, nprobe, cola, pr);
                for (k = 0; k < nprobe; k++) fpadj[(size_t)k*(size_t)D.nvar+(size_t)j] = cola[k];
            }
        }
        free(pr); free(colb); free(cola);
    }

    /* ---- emit TSV ---- */
    fputs("Probe_ID\tN", stdout);
    for (j = 0; j < D.p; j++) printf("\tEst_%s", D.coef[j]);
    for (j = 0; j < D.p; j++) printf("\tPval_%s", D.coef[j]);
    for (j = 0; j < D.nvar; j++) printf("\tFPval_%s", D.vname[j]);
    for (j = 0; j < D.nvar; j++) printf("\tEff_%s", D.vname[j]);
    for (j = 0; j < D.p; j++) printf("\tPadj_%s", D.coef[j]);
    for (j = 0; j < D.nvar; j++) printf("\tFPadj_%s", D.vname[j]);
    putchar('\n');
    for (i = 0; i < nprobe; i++) {
        printf("%s\t%d", ids[i], nobs[i]);
        for (j = 0; j < D.p; j++) { double x=est[(size_t)i*(size_t)D.p+(size_t)j];   if(isnan(x))fputs("\tNA",stdout);else printf("\t%.10g",x); }
        for (j = 0; j < D.p; j++) { double x=pval[(size_t)i*(size_t)D.p+(size_t)j];  if(isnan(x))fputs("\tNA",stdout);else printf("\t%.10g",x); }
        for (j = 0; j < D.nvar; j++){ double x=fpval[(size_t)i*(size_t)D.nvar+(size_t)j];if(isnan(x))fputs("\tNA",stdout);else printf("\t%.10g",x);}
        for (j = 0; j < D.nvar; j++){ double x=eff[(size_t)i*(size_t)D.nvar+(size_t)j]; if(isnan(x))fputs("\tNA",stdout);else printf("\t%.10g",x);}
        for (j = 0; j < D.p; j++) { double x=padj[(size_t)i*(size_t)D.p+(size_t)j];  if(isnan(x))fputs("\tNA",stdout);else printf("\t%.10g",x); }
        for (j = 0; j < D.nvar; j++){ double x=fpadj[(size_t)i*(size_t)D.nvar+(size_t)j];if(isnan(x))fputs("\tNA",stdout);else printf("\t%.10g",x);}
        putchar('\n');
    }
    rc = 0;

out:
    free(est); free(pval); free(fpval); free(eff); free(padj); free(fpadj); free(nobs);
    design_free(&D);
    if (ids) { for (i=0;i<nprobe;i++) free(ids[i]); free(ids); }
    if (samp) { for (i=0;i<nsamp;i++) free(samp[i]); free(samp); }
    free(betas);
    return rc;
}

/* Is `val` one of the comma-separated names in `csv`? (for cnv --exclude) */
static int csv_has(const char *csv, const char *val)
{
    size_t vl = strlen(val);
    const char *p = csv;
    if (!csv) return 0;
    while (*p) {
        const char *c = strchr(p, ',');
        size_t len = c ? (size_t)(c - p) : strlen(p);
        if (len == vl && strncmp(p, val, vl) == 0) return 1;
        if (!c) break;
        p = c + 1;
    }
    return 0;
}

static int cmd_cnv(int argc, char **argv)
{
    const char *normals = NULL, *idxpath = NULL, *platform = NULL;
    const char *coords = NULL, *genome = "hg38", *exclude = NULL;
    const char *pos[3] = { NULL, NULL, NULL };
    const char *target, *segout, *detout;
    int tilewidth = 50000, min_probes = 20, probes = 0, npos = 0, i, rc = 1;
    char store[4096], resolved[4096], cbuf[4096], nbuf[4096], sbuf[4096], gbuf[4096];
    sesame_index_t *ix = NULL;
    sesame_cnv_t *res = NULL;
    int32_t nres = 0, s, k;
    FILE *fseg = NULL, *fdet = NULL;
    sesame_err_t e;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--normals") == 0 && i+1 < argc) normals = argv[++i];
        else if (strcmp(argv[i], "--index") == 0 && i+1 < argc) {
            idxpath = argv[++i];
            sesame__index_deprecated();
        }
        else if (strcmp(argv[i], "--platform") == 0 && i+1 < argc) platform = argv[++i];
        else if (strcmp(argv[i], "--coords") == 0 && i+1 < argc) coords = argv[++i];
        else if (strcmp(argv[i], "--genome") == 0 && i+1 < argc) genome = argv[++i];
        else if (strcmp(argv[i], "--tilewidth") == 0 && i+1 < argc) tilewidth = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--min-probes") == 0 && i+1 < argc) min_probes = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--exclude") == 0 && i+1 < argc) exclude = argv[++i];
        else if (strcmp(argv[i], "--probes") == 0) probes = 1;
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage_cnv(); return 0; }
        else if (argv[i][0] == '-' && argv[i][1] != '\0') { fprintf(stderr, "sesame: unknown option %s\n", argv[i]); return usage_cnv(); }
        else if (npos < 3) pos[npos++] = argv[i];
        else { fprintf(stderr, "sesame: cnv takes exactly three positional args\n"); return usage_cnv(); }
    }

    /* --platform may name an ordering file instead of a platform (YAME's -R
     * rule: a path is a path). Then it IS the index, and there is no platform
     * name under which to look up companion assets. */
    if (sesame_platform_is_path(platform)) {
        if (!idxpath) idxpath = platform;
        platform = NULL;
    }
    if (npos != 3) {
        fprintf(stderr, "sesame: cnv needs <target.cg> <segments.tsv> <%s.tsv>\n",
                probes ? "probes" : "bins");
        return usage_cnv();
    }
    target = pos[0]; segout = pos[1]; detout = pos[2];
    if (tilewidth < 1 || min_probes < 1) { fprintf(stderr, "sesame: --tilewidth and --min-probes must be positive\n"); return 1; }

    /* Defaults for normals/coords/index come from the store, keyed on platform. */
    sesame_store_dir(store, sizeof store);
    if (!idxpath) {
        if (!platform) { fprintf(stderr, "sesame: cnv needs --platform (or --index)\n"); return 1; }
        if (sesame_index_locate(platform, resolved, sizeof resolved) != 0) {
            char help[1024]; sesame_index_missing_help(platform, help, sizeof help);
            fprintf(stderr, "sesame: %s\n", help); return 1;
        }
        idxpath = resolved;
    }
    if (!normals) {
        /* The normal panel is NOT part of the InfiniumAnnotation release, so
         * `yame fetch` does not supply it -- build one with `make cnv-normals`
         * or point --normals at your own. */
        char nfile[256];
        if (!platform) { fprintf(stderr, "sesame: cnv needs --normals or --platform\n"); return 1; }
        snprintf(nfile, sizeof nfile, "%s.cnvnormals.cg", platform);
        if (sesame_asset_locate(platform, nfile, nbuf, sizeof nbuf) != 0) {
            fprintf(stderr,
                "sesame: no normal panel for %s\n"
                "  searched: %s/%s/%s\n"
                "  the panel is not published with the annotation -- pass\n"
                "  --normals <file>, or build one with `make cnv-normals`\n",
                platform, store, platform, nfile);
            return 1;
        }
        normals = nbuf;
    }
    if (!coords) {
        char cfile[256];
        if (!platform) { fprintf(stderr, "sesame: cnv needs --coords or --platform\n"); return 1; }
        snprintf(cfile, sizeof cfile, "%s.%s.coord.tsv.gz", platform, genome);
        if (sesame_asset_locate(platform, cfile, cbuf, sizeof cbuf) != 0) {
            char help[1024];
            sesame_asset_missing_help(platform, cfile, help, sizeof help);
            fprintf(stderr, "sesame: %s\n", help); return 1;
        }
        coords = cbuf;
    }
    if (sesame_genome_locate(genome, "seqinfo.tsv.gz", sbuf, sizeof sbuf) != 0 ||
        sesame_genome_locate(genome, "gaps.tsv.gz", gbuf, sizeof gbuf) != 0) {
        char help[1024];
        sesame_genome_missing_help(genome, help, sizeof help);
        fprintf(stderr, "sesame: %s\n", help);
        return 1;
    }

    if (!(ix = sesame_index_open(idxpath, &e))) { fprintf(stderr, "sesame: %s\n", e.msg); return 1; }
    if (sesame_cnv_run(target, normals, ix, coords, sbuf, gbuf,
                       tilewidth, min_probes, &res, &nres, &e) != SESAME_OK) {
        fprintf(stderr, "sesame: %s\n", e.msg); goto out;
    }

    /* Always write the CBS segments; the detail file is per-bin (default) or
     * per-probe (--probes). Both are required positional outputs. cinderplot
     * reads the two files as separate layers (points + segment step). */
    if (!(fseg = fopen(segout, "w"))) { fprintf(stderr, "sesame: cannot write %s\n", segout); goto out; }
    fputs("sample\tchrom\tstart\tend\tnbin\tseg.mean\n", fseg);
    for (s = 0; s < nres; s++)
        for (k = 0; k < res[s].n_seg; k++) {
            if (csv_has(exclude, res[s].seg_chrom[k])) continue;
            fprintf(fseg, "%s\t%s\t%d\t%d\t%d\t%.6g\n", res[s].sample, res[s].seg_chrom[k],
                    res[s].seg_start[k], res[s].seg_end[k], res[s].seg_nbin[k], res[s].seg_mean[k]);
        }

    if (!(fdet = fopen(detout, "w"))) { fprintf(stderr, "sesame: cannot write %s\n", detout); goto out; }
    if (probes) {
        fputs("sample\tProbe_ID\tchrom\tpos\tlog2ratio\n", fdet);
        for (s = 0; s < nres; s++)
            for (k = 0; k < res[s].n_probe; k++) {
                if (csv_has(exclude, res[s].chrom[k])) continue;
                fprintf(fdet, "%s\t%s\t%s\t%d\t%.6g\n", res[s].sample, res[s].probe_id[k],
                        res[s].chrom[k], res[s].pos[k], res[s].log2ratio[k]);
            }
    } else {
        fputs("sample\tchrom\tstart\tend\tnprobes\tlog2ratio\n", fdet);
        for (s = 0; s < nres; s++)
            for (k = 0; k < res[s].n_bin; k++) {
                if (csv_has(exclude, res[s].bin_chrom[k])) continue;
                fprintf(fdet, "%s\t%s\t%d\t%d\t%d\t%.6g\n", res[s].sample, res[s].bin_chrom[k],
                        res[s].bin_start[k], res[s].bin_end[k], res[s].bin_nprobe[k], res[s].bin_log2ratio[k]);
            }
    }
    for (s = 0; s < nres; s++)
        fprintf(stderr, "sesame: cnv %s -- %d probes on %d normals, %d bins, %d segments\n",
                res[s].sample, res[s].n_used, res[s].n_normal, res[s].n_bin, res[s].n_seg);
    rc = 0;
out:
    if (fseg) fclose(fseg);
    if (fdet) fclose(fdet);
    sesame_cnv_free_array(res, nres);
    sesame_index_close(ix);
    return rc;
}

static int cmd_vcf(int argc, char **argv)
{
    const char *prefix = NULL, *idxpath = NULL, *platform = NULL, *snp = NULL;
    const char *genome = "hg38", *outpath = NULL;
    int min_beads = 0, variants_only = 0, i, rc = 1;
    char resolved[4096], snpbuf[4096];
    sesame_index_t *ix = NULL;
    sesame_sigdf_t *sdf = NULL;
    FILE *out = stdout;
    sesame_err_t e;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--index") == 0 && i+1 < argc) {
            idxpath = argv[++i];
            sesame__index_deprecated();
        }
        else if (strcmp(argv[i], "--platform") == 0 && i+1 < argc) platform = argv[++i];
        else if (strcmp(argv[i], "--snp") == 0 && i+1 < argc) snp = argv[++i];
        else if (strcmp(argv[i], "--genome") == 0 && i+1 < argc) genome = argv[++i];
        else if ((strcmp(argv[i], "--out") == 0 || strcmp(argv[i], "-o") == 0) && i+1 < argc) outpath = argv[++i];
        else if (strcmp(argv[i], "--min-beads") == 0 && i+1 < argc) min_beads = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--variants") == 0) variants_only = 1;
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage_vcf(); return 0; }
        else if (argv[i][0] == '-' && argv[i][1] != '\0') { fprintf(stderr, "sesame: unknown option %s\n", argv[i]); return usage_vcf(); }
        else if (!prefix) prefix = argv[i];
        else { fprintf(stderr, "sesame: vcf takes one IDAT prefix\n"); return usage_vcf(); }
    }

    /* --platform may name an ordering file instead of a platform (YAME's -R
     * rule: a path is a path). Then it IS the index, and there is no platform
     * name under which to look up companion assets. */
    if (sesame_platform_is_path(platform)) {
        if (!idxpath) idxpath = platform;
        platform = NULL;
    }
    if (!prefix) { fprintf(stderr, "sesame: vcf needs an IDAT <prefix>\n"); return usage_vcf(); }

    if (idxpath) {
        if (!(ix = sesame_index_open(idxpath, &e))) { fprintf(stderr, "sesame: %s\n", e.msg); return 1; }
    } else {
        if (!platform) {                             /* auto-detect from the Grn bead count */
            char g[4096]; sesame_idat_t *g0 = NULL;
            if (resolve_idat(prefix, "Grn", g, sizeof g) != 0) { fprintf(stderr, "sesame: cannot find %s_Grn.idat\n", prefix); return 1; }
            if (sesame_idat_read(g, &g0, &e) != SESAME_OK) { fprintf(stderr, "sesame: %s\n", e.msg); return 1; }
            platform = sesame_platform_from_beads(g0->n);
            if (!platform) { fprintf(stderr, "sesame: unknown platform (%d beads); pass --platform or --index\n", g0->n); sesame_idat_free(g0); return 1; }
            sesame_idat_free(g0);
        }
        if (sesame_index_locate(platform, resolved, sizeof resolved) != 0) {
            char help[1024]; sesame_index_missing_help(platform, help, sizeof help);
            fprintf(stderr, "sesame: %s\n", help); return 1;
        }
        if (!(ix = sesame_index_open(resolved, &e))) { fprintf(stderr, "sesame: %s\n", e.msg); return 1; }
    }

    if (!snp) {                                      /* default: the store's snp table */
        char sfile[256];
        if (!platform) { fprintf(stderr, "sesame: vcf needs --snp <file> or --platform\n"); goto out; }
        snprintf(sfile, sizeof sfile, "%s.%s.snp.tsv.gz", platform, genome);
        if (sesame_asset_locate(platform, sfile, snpbuf, sizeof snpbuf) != 0) {
            char help[1024];
            sesame_asset_missing_help(platform, sfile, help, sizeof help);
            fprintf(stderr, "sesame: %s\n", help); goto out;
        }
        snp = snpbuf;
    }

    /* raw SigDF -- genotyping must not infer channels (that erases the Type-I
     * channel-switch genotype signal). */
    if (build_sigdf_for(prefix, ix, NULL, min_beads, &sdf, &e) != SESAME_OK) {
        fprintf(stderr, "sesame: %s\n", e.msg); goto out;
    }
    if (outpath && !(out = fopen(outpath, "w"))) { fprintf(stderr, "sesame: cannot write %s\n", outpath); goto out; }
    if (sesame_format_vcf(sdf, snp, genome, variants_only, out, &e) != SESAME_OK) { fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
    rc = 0;
out:
    if (out && out != stdout) fclose(out);
    if (sdf) sesame_sigdf_free(sdf);
    sesame_index_close(ix);
    return rc;
}

static int cmd_region(int argc, char **argv)
{
    const char *regstr = NULL, *betapath = NULL, *idxpath = NULL, *platform = NULL;
    const char *coords = NULL, *genome = "hg38", *outpath = NULL;
    const char *gene = NULL, *genesbed = NULL;
    long beg = 0, end = 0, pad = 0;
    char resolved[4096], cbuf[4096], gbuf[4096], chrom[64] = "";
    sesame_index_t *ix = NULL;
    FILE *out = stdout;
    int i, rc = 1;
    sesame_err_t e;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--betas") == 0 && i+1 < argc) betapath = argv[++i];
        else if (strcmp(argv[i], "--gene") == 0 && i+1 < argc) gene = argv[++i];
        else if (strcmp(argv[i], "--gene-models") == 0 && i+1 < argc) genesbed = argv[++i];
        else if (strcmp(argv[i], "--pad") == 0 && i+1 < argc) pad = strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--index") == 0 && i+1 < argc) {
            idxpath = argv[++i];
            sesame__index_deprecated();
        }
        else if (strcmp(argv[i], "--platform") == 0 && i+1 < argc) platform = argv[++i];
        else if (strcmp(argv[i], "--coords") == 0 && i+1 < argc) coords = argv[++i];
        else if (strcmp(argv[i], "--genome") == 0 && i+1 < argc) genome = argv[++i];
        else if ((strcmp(argv[i], "--out") == 0 || strcmp(argv[i], "-o") == 0) && i+1 < argc) outpath = argv[++i];
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage_region(); return 0; }
        else if (argv[i][0] == '-' && argv[i][1] != '\0') { fprintf(stderr, "sesame: unknown option %s\n", argv[i]); return usage_region(); }
        else if (!regstr) regstr = argv[i];
        else { fprintf(stderr, "sesame: region takes one <chr:beg-end>\n"); return usage_region(); }
    }

    /* --platform may name an ordering file instead of a platform (YAME's -R
     * rule: a path is a path). Then it IS the index, and there is no platform
     * name under which to look up companion assets. */
    if (sesame_platform_is_path(platform)) {
        if (!idxpath) idxpath = platform;
        platform = NULL;
    }
    if (!betapath || (!regstr && !gene) || (regstr && gene)) {
        fprintf(stderr, "sesame: region needs --betas and exactly one of <chr:beg-end> or --gene <NAME>\n");
        return usage_region();
    }

    if (gene) {                                  /* resolve a gene symbol -> span */
        if (!genesbed) {
            if (sesame_genome_locate(genome, "genes.bed.gz", gbuf, sizeof gbuf) != 0) {
                char help[1024];
                sesame_genome_missing_help(genome, help, sizeof help);
                fprintf(stderr, "sesame: no gene models -- %s\n"
                        "  (or pass --gene-models <genes.bed.gz>)\n", help);
                return 1;
            }
            genesbed = gbuf;
        }
        if (sesame_gene_span(genesbed, gene, chrom, sizeof chrom, &beg, &end, &e) != SESAME_OK) {
            fprintf(stderr, "sesame: %s\n", e.msg); return 1;
        }
    } else {                                     /* parse chr:beg-end (commas allowed) */
        char tmp[256], *colon, *dash, *p;
        snprintf(tmp, sizeof tmp, "%s", regstr);
        colon = strchr(tmp, ':'); dash = colon ? strchr(colon + 1, '-') : NULL;
        if (!colon || !dash) { fprintf(stderr, "sesame: bad region '%s' (want chr:beg-end)\n", regstr); return 1; }
        *colon = '\0'; *dash = '\0';
        snprintf(chrom, sizeof chrom, "%s", tmp);
        for (p = colon + 1; *p; p++) if (*p >= '0' && *p <= '9') beg = beg*10 + (*p - '0');
        for (p = dash + 1;  *p; p++) if (*p >= '0' && *p <= '9') end = end*10 + (*p - '0');
    }
    if (pad) { beg -= pad; if (beg < 0) beg = 0; end += pad; }
    if (end < beg) { fprintf(stderr, "sesame: region end < beg\n"); return 1; }

    if (!idxpath) {
        if (!platform) { fprintf(stderr, "sesame: region needs --platform (or --index) for probe IDs\n"); return 1; }
        if (sesame_index_locate(platform, resolved, sizeof resolved) != 0) {
            char help[1024]; sesame_index_missing_help(platform, help, sizeof help);
            fprintf(stderr, "sesame: %s\n", help); return 1;
        }
        idxpath = resolved;
    }
    if (!coords) {
        char cfile[256];
        if (!platform) { fprintf(stderr, "sesame: region needs --coords or --platform\n"); return 1; }
        snprintf(cfile, sizeof cfile, "%s.%s.coord.tsv.gz", platform, genome);
        if (sesame_asset_locate(platform, cfile, cbuf, sizeof cbuf) != 0) {
            char help[1024];
            sesame_asset_missing_help(platform, cfile, help, sizeof help);
            fprintf(stderr, "sesame: %s\n", help); return 1;
        }
        coords = cbuf;
    }
    if (!(ix = sesame_index_open(idxpath, &e))) { fprintf(stderr, "sesame: %s\n", e.msg); return 1; }
    if (outpath && !(out = fopen(outpath, "w"))) { fprintf(stderr, "sesame: cannot write %s\n", outpath); goto out; }
    if (sesame_region_extract(betapath, ix, coords, chrom, beg, end, out, &e) != SESAME_OK) {
        fprintf(stderr, "sesame: %s\n", e.msg); goto out;
    }
    rc = 0;
out:
    if (out && out != stdout) fclose(out);
    sesame_index_close(ix);
    return rc;
}

static const char *platform_from_basename(const char *path);

static int cmd_liftover(int argc, char **argv)
{
    const char *inpath = NULL, *outpath = NULL, *src_plat = NULL, *tgt_plat = NULL;
    const char *src_idx = NULL, *tgt_idx = NULL, *coords = NULL;
    const char *how;                  /* for the report: how the map was built */
    char sres[4096], tres[4096], crbuf[4096], cobuf[4096], cfile[512], help[1024];
    sesame_index_t *six = NULL, *tix = NULL;
    sesame_rowmap_t *map = NULL, *inv = NULL;
    int depth = 0, to_genome = 0, from_genome = 0, i, rc = 1;
    sesame_err_t e;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--to") == 0 && i+1 < argc) tgt_plat = argv[++i];
        else if (strcmp(argv[i], "--platform") == 0 && i+1 < argc) src_plat = argv[++i];
        else if (strcmp(argv[i], "--simulated-depth") == 0 && i+1 < argc) depth = atoi(argv[++i]);
        else if (strcmp(argv[i], "--coords") == 0 && i+1 < argc) coords = argv[++i];
        else if (strcmp(argv[i], "--index") == 0 && i+1 < argc) {
            src_idx = argv[++i];
            sesame__index_deprecated();
        }
        else if (strcmp(argv[i], "--index-to") == 0 && i+1 < argc) {
            tgt_idx = argv[++i];
            sesame__index_deprecated();
        }
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage_liftover(); return 0; }
        else if (argv[i][0] == '-' && argv[i][1] != '\0') { fprintf(stderr, "sesame: unknown option %s\n", argv[i]); return usage_liftover(); }
        else if (!inpath) inpath = argv[i];
        else if (!outpath) outpath = argv[i];
        else { fprintf(stderr, "sesame: mliftover takes <in.cx> <out.cx>\n"); return usage_liftover(); }
    }
    if (!inpath || !outpath || !tgt_plat) {
        fprintf(stderr, "sesame: mliftover needs --to <platform|genome>, <in.cx> and <out.cx>\n");
        return usage_liftover();
    }
    if (depth < 0) { fprintf(stderr, "sesame: --simulated-depth must be >= 0\n"); return 1; }

    /* Every lift is one map applied by one routine; only the build differs.
     * A name is a genome when the store holds its cpg_nocontig.cr. */
    to_genome = sesame_asset_locate(tgt_plat, "cpg_nocontig.cr", crbuf, sizeof crbuf) == 0;
    from_genome = !to_genome && src_plat &&
        sesame_asset_locate(src_plat, "cpg_nocontig.cr", crbuf, sizeof crbuf) == 0;

    if (to_genome || from_genome) {
        /* coordinate join through the array side's coord table */
        const char *genome  = to_genome ? tgt_plat : src_plat;
        const char *plat    = to_genome ? src_plat : tgt_plat;
        const char *ordpath = to_genome ? src_idx  : tgt_idx;
        if (!plat) { fprintf(stderr, "sesame: a genome lift needs --platform for the array side\n"); return 1; }
        if (!ordpath) {
            if (sesame_index_locate(plat, sres, sizeof sres) != 0) {
                sesame_index_missing_help(plat, help, sizeof help);
                fprintf(stderr, "sesame: %s\n", help); return 1;
            }
            ordpath = sres;
        }
        if (!coords) {
            snprintf(cfile, sizeof cfile, "%s.%s.coord.tsv.gz", plat, genome);
            if (sesame_asset_locate(plat, cfile, cobuf, sizeof cobuf) != 0) {
                sesame_asset_missing_help(plat, "coord table", help, sizeof help);
                fprintf(stderr, "sesame: no %s for %s -- %s\n", cfile, plat, help); return 1;
            }
            coords = cobuf;
        }
        if (!(six = sesame_index_open(ordpath, &e))) { fprintf(stderr, "sesame: %s\n", e.msg); return 1; }
        if (sesame_rowmap_coords(coords, sesame_index_nprobes(six), crbuf, &map, &e) != SESAME_OK) {
            fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
        if (from_genome) {                       /* same pairs, columns swapped */
            if (sesame_rowmap_invert(map, &inv, &e) != SESAME_OK) { fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
            sesame_rowmap_free(map); map = inv; inv = NULL;
        }
        how = "by coordinate";
        if (!src_plat) src_plat = plat;
    } else {
        /* probe-ID prefix join between two array orderings */
        if (!src_idx) {
            if (!src_plat) { fprintf(stderr, "sesame: mliftover needs --platform or --index for the source\n"); return 1; }
            if (sesame_index_locate(src_plat, sres, sizeof sres) != 0) {
                sesame_index_missing_help(src_plat, help, sizeof help);
                fprintf(stderr, "sesame: %s\n", help); return 1;
            }
            src_idx = sres;
        } else if (!src_plat) src_plat = platform_from_basename(src_idx);
        if (!src_plat) { fprintf(stderr, "sesame: mliftover needs --platform for the source (sets the join direction)\n"); return 1; }
        if (!tgt_idx) {
            if (sesame_index_locate(tgt_plat, tres, sizeof tres) != 0) {
                sesame_index_missing_help(tgt_plat, help, sizeof help);
                fprintf(stderr, "sesame: %s\n", help); return 1;
            }
            tgt_idx = tres;
        }
        if (!(six = sesame_index_open(src_idx, &e))) { fprintf(stderr, "sesame: %s\n", e.msg); return 1; }
        if (!(tix = sesame_index_open(tgt_idx, &e))) { fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
        if (sesame_rowmap_prefix(src_plat, six, tgt_plat, tix, &map, &e) != SESAME_OK) {
            fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
        how = "by probe-ID prefix";
    }

    /* Always report the map: a wrong coordinate convention shows up here as
     * "0 of N", not as an error. */
    fprintf(stderr, "sesame: %s -> %s %s: %lld of %lld source rows map onto %lld of %lld "
            "target rows; %lld targets take more than one source\n",
            src_plat, tgt_plat, how, (long long)map->nsrc_mapped, (long long)map->nsrc,
            (long long)map->ntgt_covered, (long long)map->ntgt, (long long)map->ntgt_multi);

    if (sesame_liftover_apply(inpath, outpath, map, depth, &e) != SESAME_OK) {
        fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
    fprintf(stderr, "sesame: lifted %s -> %s to %s%s\n", inpath, tgt_plat, outpath,
            depth ? " (format 3, simulated depth on every covered row; 0,0 elsewhere)" : "");
    rc = 0;
out:
    sesame_rowmap_free(map); sesame_rowmap_free(inv);
    sesame_index_close(six); sesame_index_close(tix);
    return rc;
}

static int cmd_impute(int argc, char **argv)
{
    const char *inpath = NULL, *outpath = NULL, *method = NULL;
    const char *platform = NULL, *coords = NULL, *genome = "hg38", *axis = "probe";
    int max_neighbors = 3, i, rc = 1;
    long max_dist = 10000;
    char cbuf[4096];
    double *mat = NULL; char **names = NULL;
    int32_t nprobe = 0, nsamp = 0;
    sesame_err_t e;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--method") == 0 && i+1 < argc) method = argv[++i];
        else if (strcmp(argv[i], "--axis") == 0 && i+1 < argc) axis = argv[++i];
        else if (strcmp(argv[i], "--platform") == 0 && i+1 < argc) platform = argv[++i];
        else if (strcmp(argv[i], "--coords") == 0 && i+1 < argc) coords = argv[++i];
        else if (strcmp(argv[i], "--genome") == 0 && i+1 < argc) genome = argv[++i];
        else if (strcmp(argv[i], "--max-neighbors") == 0 && i+1 < argc) max_neighbors = (int)strtol(argv[++i],NULL,10);
        else if (strcmp(argv[i], "--max-dist") == 0 && i+1 < argc) max_dist = strtol(argv[++i],NULL,10);
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage_impute(); return 0; }
        else if (argv[i][0] == '-' && argv[i][1] != '\0') { fprintf(stderr, "sesame: unknown option %s\n", argv[i]); return usage_impute(); }
        else if (!inpath) inpath = argv[i];
        else if (!outpath) outpath = argv[i];
        else { fprintf(stderr, "sesame: impute takes <in.cg> <out.cg>\n"); return usage_impute(); }
    }
    if (!inpath || !outpath || !method) {
        fprintf(stderr, "sesame: impute needs --method <mean|neighbors>, <in.cg> and <out.cg>\n");
        return usage_impute();
    }

    if (sesame_read_cg(inpath, &mat, &nprobe, &nsamp, &names, &e) != SESAME_OK) {
        fprintf(stderr, "sesame: %s\n", e.msg); return 1; }

    if (strcmp(method, "mean") == 0) {
        int ax = !strcmp(axis, "probe") ? 1 : !strcmp(axis, "sample") ? 2 : 0;
        if (!ax) { fprintf(stderr, "sesame: --axis must be probe or sample\n"); goto out; }
        if (sesame_impute_mean(mat, nprobe, nsamp, ax, &e) != SESAME_OK) { fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
    } else if (strcmp(method, "neighbors") == 0) {
        if (!coords) {
            char cfile[256];
            if (!platform) { fprintf(stderr, "sesame: impute --method neighbors needs --coords or --platform\n"); goto out; }
            snprintf(cfile, sizeof cfile, "%s.%s.coord.tsv.gz", platform, genome);
            if (sesame_asset_locate(platform, cfile, cbuf, sizeof cbuf) != 0) {
                char help[1024];
                sesame_asset_missing_help(platform, cfile, help, sizeof help);
                fprintf(stderr, "sesame: %s\n", help); goto out;
            }
            coords = cbuf;
        }
        if (sesame_impute_neighbors(mat, nprobe, nsamp, coords, max_neighbors, max_dist, &e) != SESAME_OK) {
            fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
    } else {
        fprintf(stderr, "sesame: --method must be mean or neighbors\n"); goto out;
    }

    if (sesame_write_cg(outpath, mat, nprobe, nsamp, names, &e) != SESAME_OK) {
        fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
    fprintf(stderr, "sesame: imputed %d sample(s) x %d probes (%s) to %s\n",
            nsamp, nprobe, method, outpath);
    rc = 0;
out:
    free(mat);
    if (names) { for (i = 0; i < nsamp; i++) free(names[i]); free(names); }
    return rc;
}

/* Infer a platform from a filename like "MSA.hg38.coord.tsv.gz" -> "MSA". Returns
 * a static registry string, or NULL. */
static const char *platform_from_basename(const char *path)
{
    static const char *plats[] = { "EPICv2", "EPIC", "HM450", "MSA", NULL };
    const char *base = path_basename(path);
    for (int i = 0; plats[i]; i++) {
        size_t l = strlen(plats[i]);
        if (strncmp(base, plats[i], l) == 0 && base[l] == '.') return plats[i];
    }
    return NULL;
}

static int cmd_attach_probe(int argc, char **argv)
{
    const char *paths[64], *idxpath = NULL, *platform = NULL;
    int npath = 0;
    const char *path;
    sesame_attach_opt_t opt; memset(&opt, 0, sizeof opt);
    char resolved[4096];
    sesame_index_t *ix = NULL;
    sesame_err_t e;
    int i, rc = 1;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--index") == 0 && i+1 < argc) {
            idxpath = argv[++i];
            sesame__index_deprecated();
        }
        else if (strcmp(argv[i], "--platform") == 0 && i+1 < argc) platform = argv[++i];
        else if (strcmp(argv[i], "--with") == 0 && i+1 < argc) {
            /* comma list of ordering columns to carry through */
            char *spec = argv[++i], *t;
            for (t = strtok(spec, ","); t; t = strtok(NULL, ",")) {
                if (!strcmp(t, "col"))       opt.with |= SESAME_WITH_COL;
                else if (!strcmp(t, "M"))    opt.with |= SESAME_WITH_M;
                else if (!strcmp(t, "U"))    opt.with |= SESAME_WITH_U;
                else if (!strcmp(t, "mask")) opt.with |= SESAME_WITH_MASK;
                else if (!strcmp(t, "all"))  opt.with |= SESAME_WITH_M | SESAME_WITH_U |
                                                         SESAME_WITH_COL | SESAME_WITH_MASK;
                else {
                    fprintf(stderr, "sesame: --with: unknown column '%s' "
                        "(M, U, col, mask, all)\n", t);
                    return 1;
                }
            }
        }
        else if (strcmp(argv[i], "--all") == 0 || strcmp(argv[i], "-a") == 0) opt.all = 1;
        else if (strcmp(argv[i], "--beta") == 0) opt.beta = 1;
        else if (strcmp(argv[i], "--no-header") == 0) opt.no_header = 1;
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage_attach(); return 0; }
        else if (argv[i][0] == '-' && argv[i][1] != '\0') {
            fprintf(stderr, "sesame: unknown option %s\n", argv[i]); return usage_attach();
        } else if (npath < (int)(sizeof paths / sizeof *paths)) paths[npath++] = argv[i];
        else { fprintf(stderr, "sesame: too many input files\n"); return 1; }
    }
    if (!npath) return usage_attach();

    /* --platform may name an ordering file instead of a platform (YAME's -R
     * rule). Then it IS the index, and there are no companion assets to look
     * up under a platform name. */
    if (sesame_platform_is_path(platform)) {
        if (!idxpath) idxpath = platform;
        platform = NULL;
    }
    path = paths[0];                    /* the one we identify the ordering by */

    /* Probe IDs come from the ordering: --index wins, else --platform, else the
     * filename prefix (e.g. MSA.hg38.coord.tsv.gz), else the file's own row
     * count -- each platform's ordering has a distinct length, so a positional
     * .cg identifies itself. Filename first because it is free and already
     * documented; the row count only has to catch what it misses. */
    if (!idxpath && !platform) platform = platform_from_basename(path);

    if (!idxpath && !platform) {
        const char *pname = NULL, *fetch = NULL;
        int got = sesame_index_from_rows(path, resolved, sizeof resolved,
                                         &pname, &fetch);
        if (got == 0) {
            idxpath = resolved;
        } else if (got == 1) {
            fprintf(stderr, "sesame: %s looks like %s by row count, but that "
                "ordering is not in the store\n  fix: yame fetch %s\n",
                path, pname, fetch ? fetch : pname);
            return 1;
        } else {
            fprintf(stderr,
                "sesame: cannot tell which ordering to use for %s\n"
                "  the row count matches no published platform -- a --collapse'd\n"
                "  beta.cg is averaged to cg-prefixes and no longer lines up with\n"
                "  any ordering\n"
                "  fix: pass --index <ordering.tsv.gz>, or --platform <P>\n",
                path);
            return 1;
        }
    }

    if (!idxpath) {
        if (sesame_index_locate(platform, resolved, sizeof resolved) != 0) {
            char help[1024];
            sesame_index_missing_help(platform, help, sizeof help);
            fprintf(stderr, "sesame: %s\n", help);
            return 1;
        }
        idxpath = resolved;
    }

    if (!(ix = sesame_index_open(idxpath, &e))) {
        fprintf(stderr, "sesame: %s\n", e.msg); return 1;
    }
    if (sesame_attach_probe_n(paths, npath, ix, &opt, stdout, &e) != SESAME_OK) {
        fprintf(stderr, "sesame: %s\n", e.msg); goto out;
    }
    rc = 0;
out:
    sesame_index_close(ix);
    return rc;
}

static int cmd_idat_dump(int argc, char **argv)
{
    const char *path = NULL;
    long head = -1;
    int tsv = 0, i;
    sesame_idat_t *d = NULL;
    sesame_err_t e;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--tsv") == 0) {
            tsv = 1;
        } else if (strcmp(argv[i], "--head") == 0 && i + 1 < argc) {
            head = strtol(argv[++i], NULL, 10);
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            usage_idat(); return 0;
        } else if (argv[i][0] == '-' && argv[i][1] != '\0') {
            fprintf(stderr, "sesame: unknown option %s\n", argv[i]);
            return usage_idat();
        } else {
            path = argv[i];
        }
    }
    if (!path) return usage_idat();

    if (sesame_idat_read(path, &d, &e) != SESAME_OK) {
        fprintf(stderr, "sesame: %s: %s\n", path, e.msg);
        return 1;
    }

    if (tsv) {
        long lim = (head < 0 || head > d->n) ? d->n : head;
        for (long k = 0; k < lim; k++)
            printf("%u\t%u\t%u\t%u\n", d->addr[k], d->mean[k],
                   d->sd[k], d->nbeads[k]);
    } else {
        long lim = (head < 0) ? 5 : (head > d->n ? d->n : head);
        printf("file       %s\n", path);
        printf("version    %lld\n", (long long)d->version);
        printf("nFields    %d\n", d->n_fields);
        printf("nSNPsRead  %d\n", d->n);
        printf("\n%-12s %8s %8s %7s\n", "addr", "mean", "sd", "nbeads");
        for (long k = 0; k < lim; k++)
            printf("%-12u %8u %8u %7u\n", d->addr[k], d->mean[k],
                   d->sd[k], d->nbeads[k]);
        if (lim < d->n) printf("... (%d more)\n", d->n - (int)lim);
    }

    sesame_idat_free(d);
    return 0;
}

/* Default output path: insert `_<tag>` before _Grn/_Red (or the extension),
 * mirroring R's deIdentify naming (<name>_noid_Grn.idat). */
static void deid_outname(const char *in, const char *tag, char *out, size_t n)
{
    char stem[4096]; size_t l;
    snprintf(stem, sizeof stem, "%s", in);
    l = strlen(stem);
    if (l >= 8 && !strcmp(stem + l - 8, ".idat.gz")) stem[l - 8] = '\0';
    else if (l >= 5 && !strcmp(stem + l - 5, ".idat")) stem[l - 5] = '\0';
    l = strlen(stem);
    if (l >= 4 && !strcmp(stem + l - 4, "_Grn")) { stem[l - 4] = '\0'; snprintf(out, n, "%s_%s_Grn.idat", stem, tag); }
    else if (l >= 4 && !strcmp(stem + l - 4, "_Red")) { stem[l - 4] = '\0'; snprintf(out, n, "%s_%s_Red.idat", stem, tag); }
    else snprintf(out, n, "%s_%s.idat", stem, tag);
}

static int cmd_deidentify(int argc, char **argv)
{
    const char *input = NULL, *outpath = NULL, *idxpath = NULL, *platform = NULL, *tag;
    const char *files[2];
    char grn[4096], red[4096], resolved[4096], out[4096];
    int randomize = 0, reidentify = 0, nf = 0, i, rc = 1;
    long long seed = 0;
    sesame_index_t *ix = NULL;
    sesame_err_t e;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--index") == 0 && i+1 < argc) {
            idxpath = argv[++i];
            sesame__index_deprecated();
        }
        else if (strcmp(argv[i], "--platform") == 0 && i+1 < argc) platform = argv[++i];
        else if (strcmp(argv[i], "--out") == 0 && i+1 < argc) outpath = argv[++i];
        else if (strcmp(argv[i], "--seed") == 0 && i+1 < argc) seed = strtoll(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--randomize") == 0) randomize = 1;
        else if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--reidentify") == 0) reidentify = 1;
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage_deidentify(); return 0; }
        else if (argv[i][0] == '-' && argv[i][1] != '\0') { fprintf(stderr, "sesame: unknown option %s\n", argv[i]); return usage_deidentify(); }
        else if (!input) input = argv[i];
        else if (!outpath) outpath = argv[i];
        else { fprintf(stderr, "sesame: too many arguments\n"); return usage_deidentify(); }
    }
    tag = reidentify ? "reid" : "noid";
    if (!input) { fprintf(stderr, "sesame: deidentify needs an IDAT prefix or file\n"); return usage_deidentify(); }
    if (reidentify && seed == 0) { fprintf(stderr, "sesame: -r (reidentify) needs --seed <N> (the deidentify seed)\n"); return usage_deidentify(); }

    {   /* resolve the IDAT file(s): a .idat[.gz] is one file, else a prefix pair */
        size_t il = strlen(input);
        int is_idat = (il >= 5 && !strcmp(input + il - 5, ".idat")) ||
                      (il >= 8 && !strcmp(input + il - 8, ".idat.gz"));
        if (is_idat) files[nf++] = input;
        else {
            if (resolve_idat(input, "Grn", grn, sizeof grn) == 0) files[nf++] = grn;
            if (resolve_idat(input, "Red", red, sizeof red) == 0) files[nf++] = red;
            if (nf == 0) { fprintf(stderr, "sesame: no IDATs found for prefix %s\n", input); return 1; }
        }
    }
    if (outpath && nf != 1) { fprintf(stderr, "sesame: --out/<out> only with a single .idat input\n"); return 1; }

    if (!idxpath) {                                   /* platform -> ordering (rs M/U) */
        if (!platform) {
            sesame_idat_t *g0 = NULL; int32_t beads;
            if (sesame_idat_read(files[0], &g0, &e) != SESAME_OK) { fprintf(stderr, "sesame: %s\n", e.msg); return 1; }
            beads = g0->n; platform = sesame_platform_from_beads(beads); sesame_idat_free(g0);
            if (!platform) { fprintf(stderr, "sesame: cannot identify platform from %d beads; pass --platform or --index\n", beads); return 1; }
            fprintf(stderr, "sesame: detected %s\n", platform);
        }
        if (sesame_index_locate(platform, resolved, sizeof resolved) != 0) {
            char help[1024]; sesame_index_missing_help(platform, help, sizeof help);
            fprintf(stderr, "sesame: %s\n", help); return 1;
        }
        idxpath = resolved;
    }
    if (!(ix = sesame_index_open(idxpath, &e))) { fprintf(stderr, "sesame: %s\n", e.msg); return 1; }

    for (i = 0; i < nf; i++) {
        int r;
        if (outpath) snprintf(out, sizeof out, "%s", outpath);
        else deid_outname(files[i], tag, out, sizeof out);
        r = reidentify ? sesame_reidentify(files[i], out, ix, (uint64_t)seed, &e)
                       : sesame_deidentify(files[i], out, ix, randomize, (uint64_t)seed, &e);
        if (r != SESAME_OK) { fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
    }
    rc = 0;
out:
    sesame_index_close(ix);
    return rc;
}

int main(int argc, char **argv)
{
    if (argc < 2) return usage();
    /* `sesame help`, `sesame -h`, `sesame --help`, `sesame help <cmd>` */
    if (!strcmp(argv[1], "help") || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        if (argc >= 3) { route_usage(argv[2]); return 0; }
        usage(); return 0;
    }
    if (strcmp(argv[1], "idat-dump") == 0)
        return cmd_idat_dump(argc - 2, argv + 2);
    if (strcmp(argv[1], "preprocess") == 0)
        return cmd_preprocess(argc - 2, argv + 2);
    if (strcmp(argv[1], "dml") == 0)
        return cmd_dml(argc - 2, argv + 2);
    if (strcmp(argv[1], "cnv") == 0)
        return cmd_cnv(argc - 2, argv + 2);
    if (strcmp(argv[1], "vcf") == 0)
        return cmd_vcf(argc - 2, argv + 2);
    if (strcmp(argv[1], "deidentify") == 0)
        return cmd_deidentify(argc - 2, argv + 2);
    if (strcmp(argv[1], "region") == 0)
        return cmd_region(argc - 2, argv + 2);
    if (strcmp(argv[1], "mliftover") == 0)
        return cmd_liftover(argc - 2, argv + 2);
    if (strcmp(argv[1], "impute") == 0)
        return cmd_impute(argc - 2, argv + 2);
    if (strcmp(argv[1], "attach-probe") == 0)
        return cmd_attach_probe(argc - 2, argv + 2);
    if (strcmp(argv[1], "version") == 0 || strcmp(argv[1], "--version") == 0) {
        printf("sesame %s\n", SESAME_VERSION);
        print_build_info(stdout);
        return 0;
    }
    return usage();
}
