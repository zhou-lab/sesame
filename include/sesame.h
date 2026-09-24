/* sesame.h -- public API for libsesame
 *
 * A standalone C implementation of sesame's basic Infinium preprocessing.
 * See NUMERICS.md for documented divergences from the R implementation.
 *
 * SPDX-License-Identifier: LicenseRef-CHOP-Academic-BSD-2-Clause
 *
 * Copyright (C) 2026-present The Children's Hospital of Philadelphia
 *
 * Use of this software is available to academic and non-profit institutions
 * for research purposes under the 2-Clause BSD License; for use or transfers
 * to commercial entities, inquire with Dr. Wanding Zhou at zhouw3@chop.edu.
 * See the LICENSE file at the root of the repository for the full terms.
 */
#ifndef SESAME_H
#define SESAME_H

/* Keep in sync with conda-recipe/meta.yaml and the git release tag. */
#define SESAME_VERSION "2.1.0"

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SESAME_OK              0
#define SESAME_ERR_IO          1
#define SESAME_ERR_FORMAT      2
#define SESAME_ERR_UNSUPPORTED 3
#define SESAME_ERR_NOMEM       4

/* Errors are returned by status code and described in a caller-owned struct.
 * The library never longjmps and never writes to stderr. */
typedef struct {
    int  code;
    char msg[256];
} sesame_err_t;

/* ---------------------------------------------------------------- IDAT ---
 *
 * Illumina IDAT v3 (non-encrypted). Only the four numeric sections are read
 * (IlluminaID/Mean/SD/NBeads); string and RunInfo sections are deliberately
 * skipped because they are unreliable in gzipped and de-identified files.
 * This mirrors R/readIDAT.R:244-246.
 */
typedef struct {
    int32_t   n;        /* nSNPsRead */
    uint32_t *addr;     /* field 102, IlluminaID (bead address) */
    uint16_t *mean;     /* field 104 */
    uint16_t *sd;       /* field 103 */
    uint8_t  *nbeads;   /* field 107 */
    int64_t   version;
    int32_t   n_fields;
    int64_t   off_mean; /* byte offset of the Mean section (for de-identify rewrite) */
} sesame_idat_t;

/* Reads a .idat or .idat.gz. On success *out is a heap object owned by the
 * caller; free with sesame_idat_free. On failure returns non-zero, *out is
 * untouched, and err (if non-NULL) is filled in. */
int  sesame_idat_read(const char *path, sesame_idat_t **out, sesame_err_t *err);
void sesame_idat_free(sesame_idat_t *idat);

/* --------------------------------------------------------------- index ---
 *
 * The per-platform ordering table: Probe_ID, the M/U bead addresses, the
 * Infinium-I colour channel, and the default design mask.
 *
 * Opaque by design -- the on-disk format is an implementation detail. Today it
 * parses a TSV; a mmap-able binary can be swapped in without touching callers.
 *
 * Row order is load-bearing: it defines the canonical SigDF row order and
 * mirrors R/sesame.R:504. Never sort it.
 */
typedef struct sesame_index_t sesame_index_t;

/* Infinium design type / in-band channel. Mirrors sesame's col factor
 * levels c("G","R","2"). */
#define SESAME_COL_II 0  /* Infinium-II  (no M address; UG=meth, UR=unmeth) */
#define SESAME_COL_G  1  /* Infinium-I, green in-band */
#define SESAME_COL_R  2  /* Infinium-I, red in-band   */

/* Platform auto-detection from the IDAT bead count (nSNPsRead). Returns NULL
 * when the count is not one we have verified -- callers must then be told the
 * platform explicitly rather than have it guessed. Avoids depending on
 * sesameData's idatSignature. */
const char *sesame_platform_from_beads(int32_t beads);

/* The shared asset store, read-only as far as sesame is concerned -- `yame
 * fetch` fills it. One variable for the whole suite, because a store only
 * dedupes if every tool agrees where it is:
 *   $YAME_DATA_HOME | ${XDG_DATA_HOME:-~/.local/share}/yame
 * Returns out. */
const char *sesame_store_dir(char *out, size_t n);

/* Finds a per-platform asset at <store>/<platform>/<file>,
 * then ./<file>. 0 and fills out on success, -1 if absent. Never downloads,
 * never prompts. */
/* Warn on stderr when a store directory holds files from an EARLIER
 * annotation release than this build pins. The check and the wording are
 * yame's (yame_store_state); sesame only supplies the compiled file list.
 * Called automatically when an asset resolves out of the store. */
void sesame_warn_if_stale(const char *store_subdir);

int sesame_asset_locate(const char *platform, const char *file,
                        char *out, size_t n);

/* The directory those assets live in, whether or not it exists yet -- for
 * callers that scan it rather than name one file. */
void sesame_asset_dir(const char *platform, char *out, size_t n);

/* True when a --platform argument is really a path to an ordering file, so a
 * custom array needs no separate flag. A path is a path; anything else is a
 * platform name. */
int sesame_platform_is_path(const char *arg);

/* Finds an existing ordering table for platform. 0 and fills out on success,
 * -1 if absent. */
int sesame_index_locate(const char *platform, char *out, size_t n);

/* Identify a positional file's platform from its row count -- each ordering
 * has a distinct length -- and return that ordering's path. 0 = resolved into
 * out; 1 = identified but not in the store (*platform and *fetch are set, for
 * the error message); -1 = the count matches no known ordering. Resolves only,
 * never downloads. */
int sesame_index_from_rows(const char *file, char *out, size_t n,
                           const char **platform, const char **fetch);

/* Finds a genome-level file (e.g. "seqinfo.tsv.gz") at
 * <store>/<genome>/. These drive CNV binning (seqinfo+gaps), the CNV
 * ideogram (cytoband) and `region --gene` (genes.bed.gz). 0 on success. */
int sesame_genome_locate(const char *genome, const char *file,
                         char *out, size_t n);

/* Fill msg with actionable "not in the store, here is the yame fetch to run"
 * text, naming the asset tag this build expects. */
void sesame_index_missing_help(const char *platform, char *msg, size_t n);
void sesame_asset_missing_help(const char *platform, const char *file,
                               char *msg, size_t n);
/* 1 when the name is a genome build in the registry (hg38, mm10, ...), which
 * decides whether a missing-asset error talks about an ordering or a genome. */
int sesame_is_genome(const char *name);

void sesame_genome_missing_help(const char *genome, char *msg, size_t n);

sesame_index_t *sesame_index_open(const char *path, sesame_err_t *err);
void            sesame_index_close(sesame_index_t *ix);
int32_t         sesame_index_nprobes(const sesame_index_t *ix);
const char     *sesame_index_probe_id(const sesame_index_t *ix, int32_t i);

/* ----------------------------------------------------------- describe ---
 *
 * Attach the ordering's Probe_IDs to a positional data file, writing a labeled
 * TSV to `out`. The file's rows are positionally aligned to the ordering (the
 * lineage MUST match -- use the same platform/tag that produced the file), so
 * row i is labeled with ix's i-th Probe_ID. Handles YAME .cg/.cm/.cx (any
 * format: fmt0 mask bit, fmt3 M/U or beta, fmt4 float, ...) and plain text
 * .tsv[.gz] (e.g. <platform>.hg38.coord.tsv.gz), whose own header is kept and
 * prefixed with "Probe_ID". Errors if the row count does not match ix. */
/* Ordering columns describe-probe can emit beside the Probe_ID. The ordering is
 * already open to supply the IDs; these let a caller have the rest of the row
 * without re-reading the file and joining it back on. */
#define SESAME_WITH_M     (1u << 0)   /* methylated bead address (NA -> "NA") */
#define SESAME_WITH_U     (1u << 1)   /* unmethylated bead address            */
#define SESAME_WITH_COL   (1u << 2)   /* design/in-band channel: G, R or 2    */
#define SESAME_WITH_MASK  (1u << 3)   /* the ordering's mask bit, when it has one */

typedef struct {
    int all;        /* YAME: emit every sample/record, not just the first     */
    int beta;       /* format 3: print beta instead of M<TAB>U                 */
    int no_header;  /* suppress the header line (text: treat line 1 as data)   */
    unsigned with;  /* SESAME_WITH_* bits, emitted right after Probe_ID        */
} sesame_describe_opt_t;

int sesame_describe_probe(const char *path, const sesame_index_t *ix,
                        const sesame_describe_opt_t *opt, FILE *out,
                        sesame_err_t *err);

/* Several YAME files side by side. They are positional to the same ordering,
 * so this is a column concatenation, not a join -- row i is probe i in every
 * one of them. Row counts must agree. Text input takes a single file only. */
/* describe-probe's coordinate mode: probe IDs in, "<Probe_ID>\t<chrm>_<beg1>"
 * out, in input order, one line per resolved probe (replicates give several).
 * The coordinate is 1-based, which is what `yame rowsub -L` reads.
 *
 * A bare cg number matches every design suffix (EPICv2/MSA); an ID the
 * platform does not carry is an error, since it is a mistake in the query. A
 * probe with no CpG in this genome, or one on a non-primary contig that a
 * genome-indexed store cannot address, is dropped and counted. */
typedef struct {
    int64_t n_query;      /* IDs read */
    int64_t n_out;        /* lines written */
    int64_t n_unmapped;   /* resolved, but no coordinate in this genome */
    int64_t n_altcontig;  /* resolved to an alt/random/fix contig */
} sesame_describe_stat_t;

int sesame_describe_coords(const char *ids_path, const sesame_index_t *ix,
                           const char *coords_path, FILE *out,
                           sesame_describe_stat_t *st, sesame_err_t *err);

int sesame_describe_probe_n(const char *const *paths, int npath,
                          const sesame_index_t *ix,
                          const sesame_describe_opt_t *opt, FILE *out,
                          sesame_err_t *err);

/* --------------------------------------------------------------- sigdf ---
 *
 * Signals, one row per index probe, in index order. NaN encodes NA.
 */

/* Status bits. R signals these conditions silently or not at all; we surface
 * them so a pipeline can count fallbacks instead of finding them in the betas. */
#define SESAME_STAT_ADDR_MISSING   (1u << 0) /* manifest address absent from IDAT */
#define SESAME_STAT_PAIR_MISMATCH  (1u << 1) /* Grn/Red address vectors disagree  */
#define SESAME_STAT_DYEBIAS_FAILED (1u << 2) /* D gave up: green channel failed   */
#define SESAME_STAT_NOOB_SKIPPED   (1u << 3) /* B skipped: <100 background probes  */
#define SESAME_STAT_NOOB_MAD0      (1u << 4) /* B: huber MAD==0 fallback (D6)       */

typedef struct {
    const sesame_index_t *ix;
    int32_t   n;
    double   *MG, *MR, *UG, *UR;
    uint8_t  *col;    /* SESAME_COL_* */
    uint8_t  *mask;
    uint32_t  status; /* SESAME_STAT_* */
    int32_t   n_addr_missing;
} sesame_sigdf_t;

/* Assembles a SigDF from a Grn/Red IDAT pair against the index.
 *
 * Divergence D1: R merges the two IDATs by *position* with no join on
 * IlluminaID (R/sesame.R:293-298). That assumption holds for well-formed files
 * -- verified across every test array -- but fails silently and array-wide when
 * mismatched files are paired. We verify it and error instead.
 *
 * min_beads: 0 disables bead filtering (R's default is NULL/off). Otherwise a
 * probe is masked when any relevant bead count is missing or < min_beads. */
sesame_sigdf_t *sesame_sigdf_from_idats(const sesame_idat_t *grn,
                                        const sesame_idat_t *red,
                                        const sesame_index_t *ix,
                                        int min_beads,
                                        sesame_err_t *err);
void sesame_sigdf_free(sesame_sigdf_t *sdf);

/* Deep copy, sharing the read-only index (NULL on OOM). */
sesame_sigdf_t *sesame_sigdf_dup(const sesame_sigdf_t *sdf);

/* --------------------------------------------------------------- prep ---
 *
 * prepSesame steps. Order matters (R/open.R:30-35). Q, C, D, E, P, B are all
 * implemented; detectionPnegEcdf is exposed separately (preprocess --detection).
 */

/* Q -- load the recommended quality mask for a platform as a 0/1 vector aligned
 * to the ordering (1 = masked). Shells out to the `yame` binary to read the .cm
 * mask in the store; *out is malloc'd (caller frees), *out_n is the probe count.
 * yame runs as a separate process here; other paths link libyame directly. */
int sesame_quality_mask(const char *platform, const char *maskpath,
                        uint8_t **out, int32_t *out_n, sesame_err_t *err);

/* Q -- apply a mask vector from sesame_quality_mask to a SigDF (mask |= qmask).
 * qn must equal the probe count. */
int sesame_prep_quality_mask(sesame_sigdf_t *sdf, const uint8_t *qmask,
                             int32_t qn, sesame_err_t *err);

/* P -- the background mask (backgroundMask names) as a 0/1 vector aligned to the
 * ordering (1 = excluded from background estimation). Reads the same .cm. */
int sesame_background_mask(const char *platform, const char *maskpath,
                           uint8_t **out, int32_t *out_n, sesame_err_t *err);

/* P -- pOOBAH (R/detection.R:141-170). Masks probes whose detection p-value
 * exceeds pval_threshold, from the ecdf of out-of-band + negative-control
 * background (probes in bgmask excluded). combine_neg adds negative controls.
 * Implements D2 (NA-channel probes are masked, not silently kept). */
int sesame_prep_poobah(sesame_sigdf_t *sdf, const uint8_t *bgmask, int32_t bn,
                       double pval_threshold, int combine_neg, sesame_err_t *err);

/* The pOOBAH detection p-value per probe into pout[nprobes] (no masking); shared
 * by the P step and the `pval` output. D2: no-signal probes get p=1. */
int sesame_poobah_pvals(const sesame_sigdf_t *sdf, const uint8_t *bgmask, int32_t bn,
                        int combine_neg, double *pout, sesame_err_t *err);

/* detectionPnegEcdf (R/detection.R:100-118): detection p per probe from the
 * empirical CDF of the negative-control probes (funcG over their UG, funcR over
 * their UR); p = pmin of the two 1-F channel values (na.rm), no-signal -> 1.
 * An alternative to pOOBAH needing only negative controls. Into pout[nprobes]. */
int sesame_detection_pneg_pvals(const sesame_sigdf_t *sdf, double *pout,
                                sesame_err_t *err);

/* Per-sample pipeline: apply prep to a copy of raw, filling whichever of
 * beta/M/U/total/pval are non-NULL (each nprobes). Signal is from the
 * preprocessed copy unless raw_signal. Shared by `preprocess` and validation. */
int sesame_pipeline(const sesame_sigdf_t *raw, const char *prep,
                    const uint8_t *qmask, int32_t qn,
                    const uint8_t *bgmask, int32_t bgn, int raw_signal,
                    int pneg_detection,
                    double *beta, double *M, double *U, double *total,
                    double *pval, uint8_t *col, sesame_err_t *err);

/* C -- inferInfiniumIChannel (R/channel_inference.R:20-55). Reassigns each
 * Infinium-I probe to its brighter channel, judged against the 95th percentile
 * of the pooled out-of-band signal under the new assignment. Defaults
 * (switch_failed=0, mask_failed=0) revert failed probes to the manifest channel
 * and do not mask them. */
int sesame_prep_infer_channel(sesame_sigdf_t *sdf, int switch_failed,
                              int mask_failed, sesame_err_t *err);

/* D -- dyeBiasNL (R/dye_bias.R:118-167). Pulls the Red and Green Infinium-I
 * distributions to their common midpoint. If the green channel has failed
 * (RGdistort NA or > 10) it masks all Inf-I green probes, sets
 * SESAME_STAT_DYEBIAS_FAILED and returns -- R does this silently. */
int sesame_prep_dye_bias_nl(sesame_sigdf_t *sdf, sesame_err_t *err);

/* E -- dyeBiasL (R/dye_bias.R:191-208). Linear dye-bias correction: scale each
 * colour channel so the median in-band Inf-I intensity meets a common reference
 * (the sample's mean in-band intensity). Only Inf-II betas change (per-channel
 * scaling leaves Inf-I M/(M+U) invariant). Both medians and the reference use
 * unmasked probes (signalMU/meanIntensity mask=TRUE). */
int sesame_prep_dye_bias_l(sesame_sigdf_t *sdf, sesame_err_t *err);

/* B -- noob (R/background.R:85-123). Normal-exponential background subtraction:
 * background (out-of-band + negative controls, probes in bgmask excluded) is
 * modelled Normal, true signal Exponential, and each channel's signal is
 * deconvolved via the inverse Mills ratio (D5), then shifted by offset (R's
 * default 15). combine_neg adds negative controls to the background, as in P.
 * If either channel has <100 positive background values, R returns the SigDF
 * unchanged; we do the same and set SESAME_STAT_NOOB_SKIPPED. D6: a huber MAD of
 * 0 sets SESAME_STAT_NOOB_MAD0 instead of aborting. */
int sesame_prep_noob(sesame_sigdf_t *sdf, const uint8_t *bgmask, int32_t bn,
                     int combine_neg, double offset, sesame_err_t *err);

/* Beta values, one per probe, in index order. NaN for NA and for masked probes
 * when apply_mask is non-zero. out must hold sesame_index_nprobes() doubles. */
int sesame_get_betas(const sesame_sigdf_t *sdf, int apply_mask,
                     double *out, sesame_err_t *err);

/* betasCollapseToPfx (R/sesame.R:157): average betas that share a probe-ID
 * prefix (the part before the first '_'), collapsing EPICv2/MSA replicate probes
 * to one value per prefix. mat_in is sample-major [nsamp*nprobes]; *mat_out is a
 * fresh sample-major [nsamp*(*m_out)] over the *m_out sorted unique prefixes
 * returned in *pfx_out. NaN-skipping mean; an all-NaN group is NaN. Caller frees
 * *mat_out, each *pfx_out[k], and *pfx_out. */
int sesame_betas_collapse_prefix(const sesame_index_t *ix,
                                 const double *mat_in, int32_t nsamp,
                                 double **mat_out, char ***pfx_out, int32_t *m_out,
                                 sesame_err_t *err);

/* mLiftOver (R/mLiftOver.R, mapping=NULL/impute=FALSE path): lift betas from the
 * source platform/ordering to the target platform/ordering by a probe-ID prefix
 * join (the modern side -- EPICv2/MSA -- has its _suffix stripped when crossing
 * to/from a legacy platform; same-family joins on the full ID). Each target
 * probe takes the first prefix-matched source beta, else NaN. mat_in is
 * sample-major [nsamp*src_nprobes]; *mat_out is a fresh sample-major
 * [nsamp*tgt_nprobes], positional to the target ordering. Caller frees. */
int sesame_liftover_betas(const char *src_platform, const sesame_index_t *src_ix,
                          const char *tgt_platform, const sesame_index_t *tgt_ix,
                          const double *mat_in, int32_t nsamp,
                          double **mat_out, sesame_err_t *err);

/* A row map: an index-to-index mapping between a SOURCE row space and a
 * TARGET row space, as (src, tgt) pairs sorted by target. One representation
 * for every lift -- one-to-one, many-to-one (array replicates on one CpG),
 * one-to-many (one CpG under several probes), and misses (a row in no pair).
 * Built three ways: sesame_rowmap_prefix (array -> array, the probe-ID prefix
 * join of R's mLiftOver), sesame_rowmap_coords (array -> genome, joined on
 * coordinate: <plat>.<genome>.coord.tsv.gz against <genome>/cpg_nocontig.cr,
 * through YAME's row index; ~2 s, in memory, nothing persisted), and
 * sesame_rowmap_invert (the same pairs with the columns swapped, which is the
 * genome -> array map for free). The counts are for reporting: a wrong
 * coordinate convention presents as "0 of N mapped", never as an error.
 *
 * The coordinate join lifts cg probes ONLY, decided by Probe_ID, never by
 * where the probe maps: an rs probe's beta is an allele fraction and an nv
 * probe's a variant fraction, not methylation, and a multi-mapping rs probe
 * whose coordinate lands on a CpG row (329 on MSA, 6 on EPICv2, plus 78 nv)
 * would otherwise write that fraction into the CpG -- or, lifting a mask
 * the other way, take a pattern state it has no business in. Those probes
 * are counted in nsrc_skipped, apart from the cg probes with no CpG row.
 *
 * sesame_liftover_apply streams any format 0/1/3/4/5/6 record (.cg or .cm)
 * through a map: each run of equal targets is reduced onto that target --
 * betas by the mean of non-NA sources, real M/U counts pooled, bits OR'd,
 * 2-bit codes by max, other formats first-wins -- and a target in no pair
 * keeps the format's "absent" (NA / 0,0 / bit 0 / code 0). Output is
 * positional to the target space. depth > 0 writes format 3 with M+U == depth
 * on every covered row (betas re-encoded as pseudo-counts, what a whole-genome
 * model such as methscope reads); depth == 0 keeps the input's format.
 * nthreads workers share the one map in-process (0 = the machine's CPUs);
 * records are independent, so the output is byte-identical at any count. */
typedef struct {
    int64_t nsrc, ntgt;               /* rows in the source / target space   */
    int64_t npair;
    int64_t *src, *tgt;               /* the pairs, sorted by (tgt, src)     */
    int64_t nsrc_mapped;              /* distinct sources in some pair       */
    int64_t nsrc_skipped;             /* sources not lifted by TYPE: a coord */
                                      /* join lifts cg probes only, by ID    */
    int64_t ntgt_covered;             /* distinct targets in some pair       */
    int64_t ntgt_multi;               /* targets with >1 source (reduced)    */
} sesame_rowmap_t;
int sesame_rowmap_prefix(const char *src_platform, const sesame_index_t *src_ix,
                         const char *tgt_platform, const sesame_index_t *tgt_ix,
                         sesame_rowmap_t **out, sesame_err_t *err);
int sesame_rowmap_coords(const char *coords_path, const sesame_index_t *ix,
                         const char *cr_path, sesame_rowmap_t **out,
                         sesame_err_t *err);
int sesame_rowmap_invert(const sesame_rowmap_t *in, sesame_rowmap_t **out,
                         sesame_err_t *err);
void sesame_rowmap_free(sesame_rowmap_t *m);
int sesame_liftover_apply(const char *in_cx, const char *out_cx,
                          const sesame_rowmap_t *m, int depth, int nthreads,
                          sesame_err_t *err);

/* imputeBetasMatrixByMean (R/impute.R): fill NaN in a sample-major matrix
 * [nsamp*nprobe] in place with axis=1 (per probe, mean across samples) or axis=2
 * (per sample, mean across probes). An all-NaN row/column is left NaN. */
int sesame_impute_mean(double *mat, int32_t nprobe, int32_t nsamp, int axis,
                       sesame_err_t *err);

/* imputeBetasByGenomicNeighbors (R/impute.R): fill each missing probe's NaN with
 * the mean of up to max_neighbors nearest non-missing probes within max_dist bp,
 * searched 5'->3' along the probe's strand (per R's start-anchored strand-aware
 * resize + findOverlaps). coords_path is the positional <plat>.<genome>.coord
 * table (CpG_chrm, CpG_beg, strand). In place, sample-major [nsamp*nprobe]. */
int sesame_impute_neighbors(double *mat, int32_t nprobe, int32_t nsamp,
                            const char *coords_path, int max_neighbors,
                            long max_dist, sesame_err_t *err);

/* Per-probe methylated/unmethylated signal (R's signalMU), index order: col G ->
 * (MG,UG), col R -> (MR,UR), Inf-II -> (UG,UR). NaN preserved. Either out NULL. */
int sesame_signal_mu(const sesame_sigdf_t *sdf, double *M, double *U,
                     sesame_err_t *err);

/* Total intensity (M+U) per probe, in index order (R's totalIntensities,
 * mask=FALSE). NaN if either allele is NA. The CNV signal input. */
int sesame_total_intensities(const sesame_sigdf_t *sdf, double *out,
                             sesame_err_t *err);

/* Genotype the SNP (rs) probes and write a VCF (sesame's formatVCF). sdf should
 * be the RAW SigDF -- channel inference would erase the Type-I channel-switch
 * genotype signal. snp_path is the platform's <platform>.<genome>.snp.tsv.gz
 * (columns chrm beg end strand rs designType U REF ALT Probe_ID). Sites-only
 * VCF with per-probe genotype/score in INFO, sorted by chrom then position.
 * variants_only drops the non-switching Infinium-I probes (REF_InfI with no rs)
 * -- the bulk that carry no known variant -- keeping rs and channel-switch probes. */
int sesame_format_vcf(const sesame_sigdf_t *sdf, const char *snp_path,
                      const char *genome, int variants_only, FILE *out,
                      sesame_err_t *err);

/* Extract a genomic region's betas from a multi-sample beta .cg as long-form
 * (tidy) TSV -- one row per (probe in [chrom:beg-end]) x (sample), columns
 * chrom/beg/end/Probe_ID/beta/sample_name. beta_cg, ix, and coords_path are all
 * positional in the same ordering. The plot-ready feed for a region/locus view. */
int sesame_region_extract(const char *beta_cg, const sesame_index_t *ix,
                          const char *coords_path, const char *chrom,
                          long beg, long end, FILE *out, sesame_err_t *err);

/* De-identify an IDAT by removing the genetic fingerprint: the SNP (rs) probes'
 * mean intensities are zeroed (randomize=0) or reversibly scrambled with `seed`
 * (randomize=1), and a new IDAT is written to out_path (everything else copied
 * byte-for-byte). Ports sesame's deIdentify. `ix` supplies the rs-probe M/U
 * addresses. reIdentify restores a scrambled file given the same seed. NOTE: the
 * scramble PRNG is this tool's own, so a randomized file round-trips with sesame
 * reidentify --seed, not with R (zeroing is fully R-equivalent and irreversible). */
int sesame_deidentify(const char *in_path, const char *out_path,
                      const sesame_index_t *ix, int randomize, uint64_t seed,
                      sesame_err_t *err);
int sesame_reidentify(const char *in_path, const char *out_path,
                      const sesame_index_t *ix, uint64_t seed, sesame_err_t *err);

/* Resolve a gene symbol to its span [beg,end) on chrom from a GENCODE BED12+
 * gene-model file (genes.bed.gz; field 14 = gene_name). A symbol on several
 * chromosomes resolves to the one with the most transcripts. Feeds
 * sesame_region_extract so `region --gene ADA` works without a chr:beg-end. */
int sesame_gene_span(const char *genes_bed, const char *gene,
                     char *chrom, size_t chrom_n, long *beg, long *end,
                     sesame_err_t *err);

/* Write sample-major matrices as a YAME .cg (+ <path>.idx of sample names),
 * consumable by the `yame` toolchain. Mat is sample-major: sample j, probe i at
 * mat[j*nprobe+i].
 *
 * _mu writes format 3 (M/U counts, the default) -- yame derives beta (MU2beta)
 * and total intensity (MU2cov); M,U round to integers. The plain writer writes
 * format 4 (one exact float32 per probe, NA=-1.0), on request. */
int sesame_write_cg_mu(const char *path, const double *matM, const double *matU,
                       int32_t nprobe, int32_t nsamp, char *const *names,
                       sesame_err_t *err);
int sesame_write_cg(const char *path, const double *mat, int32_t nprobe,
                    int32_t nsamp, char *const *names, sesame_err_t *err);

/* Read a format-4 .cg (+ .idx names) into a sample-major matrix (sample j, probe
 * i at mat[j*nprobe+i]); NA float -> NaN. Caller frees mat, each name, names. */
int sesame_read_cg(const char *path, double **mat_out, int32_t *nprobe_out,
                   int32_t *nsamp_out, char ***names_out, sesame_err_t *err);

/* Read a .cg's per-probe TOTAL intensity, format-agnostic: format 4 -> the
 * stored float, format 3 -> M+U. NA -> NaN. Sample-major, like sesame_read_cg.
 * Used by CNV for both the target (total_intensity.cg) and the normal reference
 * (cnvnormals.cg, format 3). Caller frees mat, each name, names. */
int sesame_read_cg_total(const char *path, double **mat_out, int32_t *nprobe_out,
                         int32_t *nsamp_out, char ***names_out, sesame_err_t *err);

/* ----------------------------------------------------------------- QC ---
 *
 * The sesameQC panel (R sesameQC_calcStats), computed from a *raw* SigDF. Each
 * metric is one field, tagged I (integer count) or D (real) for formatting. The
 * groups mirror R/QC.R: detection, numProbes, intensity, channel, dyeBias (incl.
 * GCT, bisConversionControl), betas.
 *
 * Two faithful-but-noted differences (see NUMERICS.md):
 *  - the detection group runs pOOBAH internally, so it carries the same mask
 *    lineage as `P`, and `num_dtna` counts probes with no signal in *either*
 *    channel (the D2 fix makes R's pOOBAH no longer return NA there);
 *  - the betas group is computed on a D->B->P-processed copy, exactly as
 *    sesameQC_calcStats_betas does (getBetas(pOOBAH(noob(dyeBiasNL(sdf))))).
 */
#define SESAME_QC_FIELDS(_) \
    _(I, num_dtna)      _(D, frac_dtna) \
    _(I, num_dt)        _(D, frac_dt) \
    _(I, num_dt_mk)     _(D, frac_dt_mk) \
    _(I, num_dt_cg)     _(D, frac_dt_cg) \
    _(I, num_dt_ch)     _(D, frac_dt_ch) \
    _(I, num_dt_rs)     _(D, frac_dt_rs) \
    _(I, num_probes)    _(I, num_probes_II) _(I, num_probes_IR) _(I, num_probes_IG) \
    _(I, num_probes_cg) _(I, num_probes_ch) _(I, num_probes_rs) \
    _(D, mean_intensity)   _(D, mean_intensity_MU) _(D, mean_ii) \
    _(D, mean_inb_grn)     _(D, mean_inb_red) \
    _(D, mean_oob_grn)     _(D, mean_oob_red) \
    _(I, na_intensity_M)   _(I, na_intensity_U) \
    _(I, na_intensity_ig)  _(I, na_intensity_ir)  _(I, na_intensity_ii) \
    _(I, InfI_switch_R2R)  _(I, InfI_switch_G2G) \
    _(I, InfI_switch_R2G)  _(I, InfI_switch_G2R) \
    _(D, medR) _(D, medG) _(D, topR) _(D, topG) _(D, RGratio) _(D, RGdistort) _(D, GCT) \
    _(D, mean_beta)    _(D, median_beta)    _(D, frac_unmeth)    _(D, frac_meth)    _(I, num_na)    _(D, frac_na) \
    _(D, mean_beta_cg) _(D, median_beta_cg) _(D, frac_unmeth_cg) _(D, frac_meth_cg) _(I, num_na_cg) _(D, frac_na_cg) \
    _(D, mean_beta_ch) _(D, median_beta_ch) _(D, frac_unmeth_ch) _(D, frac_meth_ch) _(I, num_na_ch) _(D, frac_na_ch) \
    _(D, mean_beta_rs) _(D, median_beta_rs) _(D, frac_unmeth_rs) _(D, frac_meth_rs) _(I, num_na_rs) _(D, frac_na_rs)

typedef struct {
#define SESAME_QC_DECL(t, nm) double nm;
    SESAME_QC_FIELDS(SESAME_QC_DECL)
#undef SESAME_QC_DECL
} sesame_qc_t;

/* Compute the panel from a raw SigDF. bgmask is the background mask (as for P/B);
 * required, since the detection and beta groups run pOOBAH internally. `ext`
 * (optional, positional to the ordering, length extn) carries the Type-I
 * extension-base class per probe -- 1 = extension-C, 2 = extension-A/T, 0 = none
 * -- for the GCT bisulfite-conversion metric; pass NULL (extn 0) to leave GCT
 * NaN (no ext asset). The lists come from the annotation's <plat>.typeI_ext, a
 * uniform superset of R's probeInfo curation (see NUMERICS.md). Does not modify
 * sdf. */
int sesame_qc_calc(const sesame_sigdf_t *sdf, const uint8_t *bgmask, int32_t bn,
                   const uint8_t *ext, int32_t extn,
                   sesame_qc_t *out, sesame_err_t *err);

/* Load positional Type-I extension-base codes (<plat>.typeI_ext.tsv.gz) for GCT:
 * C/T/. -> 1/2/0, one per probe in ordering order. *out is malloc'd (caller
 * frees). Errors (incl. absent file) leave GCT NaN at the call site. */
int sesame_load_ext_codes(const char *path, int32_t np, uint8_t **out,
                          sesame_err_t *err);

/* Tab-separated metric-name header, no leading sample column and no trailing
 * tab. */
const char *sesame_qc_header(void);

/* Format the metric values tab-separated into buf (integer metrics as integers,
 * reals as %.10g, NaN as "NA"). Returns the length written, or -1 if truncated. */
int sesame_qc_format_row(const sesame_qc_t *q, char *buf, size_t n);

/* ---------------------------------------------------------------- DML ---
 *
 * Per-probe differential methylation: OLS of a probe's betas across samples on a
 * design matrix, with per-coefficient t-tests and a holdout F-test per
 * categorical variable (sesame's DML). No genomic coordinates needed; region
 * calling (DMR) is separate.
 */

/* The design. X is m*p, row-major (sample i, coef j at X[i*p+j]); the caller
 * builds it from a formula or supplies it directly. Each of the nvar categorical
 * variables owns the coefficient columns var_col[var_off[v] .. var_off[v+1]) --
 * used for the holdout F-test and effect size; continuous terms and the
 * intercept belong to no variable. */
typedef struct {
    int32_t         m, p, nvar;
    const double   *X;
    const int32_t  *var_off;   /* length nvar+1 */
    const int32_t  *var_col;   /* length var_off[nvar] */
} sesame_dml_design_t;

/* Per-thread scratch, sized for one design. */
typedef struct sesame_dml_work_t sesame_dml_work_t;
sesame_dml_work_t *sesame_dml_work_new(int32_t m, int32_t p);
void               sesame_dml_work_free(sesame_dml_work_t *w);

/* Fit one probe. y[m] is the probe's betas across samples; NaN samples are
 * dropped. Fills est[p], pval[p] (per coefficient) and fpval[nvar], eff[nvar]
 * (per categorical variable) -- NaN where not estimable. Returns the number of
 * samples used, or 0 if unfittable (too few observations, or rank-deficient). */
int32_t sesame_dml_fit(const sesame_dml_design_t *d, sesame_dml_work_t *w,
                       const double *y, double *est, double *pval,
                       double *fpval, double *eff);

/* ---------------------------------------------------------------- CNV ---
 *
 * Copy-number: for one target sample, regress its per-probe total intensity on a
 * panel of normal samples' total intensities (OLS, with intercept), then take
 * log2(target / max(fitted, 1)) per probe -- the copy-number signal, mirroring
 * sesame's cnSegmentation. Probes are then binned along the genome (tile,
 * subtract assembly gaps, left/right-merge to >= min_probes) and each bin gets
 * the median probe signal. CBS segmentation (DNAcopy) is NOT included -- the bin
 * log2 ratios are the profile; feed them to a segmenter downstream if needed. */
typedef struct {
    char     *sample;          /* target sample name */
    int32_t   n_normal;        /* normals in the panel */
    int32_t   n_used;          /* probes used in the fit (signal + coord) */
    /* per-probe (length n_probe) */
    int32_t   n_probe;
    char    **probe_id;
    char    **chrom;
    int32_t  *pos;
    double   *log2ratio;
    /* per-bin (length n_bin) */
    int32_t   n_bin;
    char    **bin_chrom;
    int32_t  *bin_start;
    int32_t  *bin_end;
    int32_t  *bin_nprobe;
    double   *bin_log2ratio;
    /* per-segment (length n_seg): circular binary segmentation over the bins */
    int32_t   n_seg;
    char    **seg_chrom;
    int32_t  *seg_start;
    int32_t  *seg_end;
    int32_t  *seg_nbin;
    double   *seg_mean;
} sesame_cnv_t;

/* Run CNV for every sample in target_cg against normals_cg. target_cg is a
 * total_intensity.cg (format 4) or any .cg whose value is total intensity;
 * normals_cg is the format-3 normal reference (M+U). Both are positional in the
 * ordering `ix`. coords_path is the per-probe <platform>.hg38.coord.tsv.gz;
 * seqinfo_path/gaps_path are the fetched genome tables. Returns an array of
 * n_out results (one per target sample); free with sesame_cnv_free_array. */
int sesame_cnv_run(const char *target_cg, const char *normals_cg,
                   const sesame_index_t *ix, const char *coords_path,
                   const char *seqinfo_path, const char *gaps_path,
                   int32_t tilewidth, int32_t min_probes,
                   sesame_cnv_t **out, int32_t *n_out, sesame_err_t *err);
void sesame_cnv_free_array(sesame_cnv_t *a, int32_t n);

#ifdef __cplusplus
}
#endif

#endif /* SESAME_H */
