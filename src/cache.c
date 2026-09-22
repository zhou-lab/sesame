/* cache.c -- locating annotation in the shared store.
 *
 * sesame does not download anything. YAME owns fetching for the whole suite --
 * `yame fetch` writes the store, every tool that links libyame reads it -- so
 * this file is only path resolution plus the errors that name the command to
 * run. See YAME/src/assets.h.
 *
 * The store is keyed on the browser hierarchy -- species, then platform or
 * genome build -- which is the address a user navigates to, so a file one tool
 * downloads is found by the next:
 *
 *   <store>/MSA/SHA256SUMS         byte-identical copy of the remote
 *   <store>/MSA/MSA.ordering.tsv.gz
 *   <store>/MSA/MSA.hg38.mask.cm
 *   <store>/MSA/KYCG/...           kycg's per-platform sets
 *   <store>/hg38/seqinfo.tsv.gz    genome-level annotation
 *   <store>/hg38/KYCG/...          kycg's sequencing knowledgebase
 *
 * Note the store path is NOT the fetch address: `yame fetch` names a unit as
 * <source>/<target> (InfiniumAnnotation/MSA, genomes/hg38), because one
 * directory can be filled from more than one upstream -- <store>/hg38 holds
 * both the genomes files and, under KYCG/, the KYCGKB sets.
 *
 * Before this, sesame kept its own cache under ~/.cache/sesame and kycg kept one
 * under ~/.cache/kycg, and the two downloaded the identical InfiniumAnnotation
 * files separately.
 *
 * Mirroring means the store is self-describing and hand-verifiable per platform
 * (cd <store>/MSA && shasum -a 256 -c SHA256SUMS), with no
 * local bookkeeping -- and that stored SHA256SUMS is how a directory says which
 * upstream tag filled it. registry.h records the tag THIS build expects, which
 * is why the errors below quote it.
 *
 * Store location -- one rule, yame_assets_root():
 *
 *   $YAME_DATA_HOME, else ${XDG_DATA_HOME:-~/.local/share}/yame
 *
 * One variable for the whole suite, because a store only dedupes if every tool
 * agrees where it is. The data tier rather than a cache tier is deliberate:
 * these are multi-GB references that should survive somebody clearing ~/.cache.
 * An existing ~/.cache/sesame is detected and named in a one-line notice by
 * libyame, never read from and never moved.
 *
 * Resolution order (explicit always wins; same rule for library and CLI):
 *
 *   1. --index / --coords / --normals <path>          (caller-supplied)
 *   2. <store>/<platform>/<file>
 *   3. ./<file>
 *
 * sesame NEVER prompts and NEVER touches the network. If an asset is missing,
 * the caller gets an error naming the exact `yame fetch` to run. An interactive
 * prompt would hang forever in a Nextflow job or a Docker build -- silently, on
 * question one, across a whole run.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present Wanding Zhou
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "sesame.h"
#include "internal.h"
#include "registry.h"
#include "assets.h"                  /* YAME: the shared store */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* `yame fetch` addresses a unit as <source>/<target> (InfiniumAnnotation/MSA),
 * but the STORE is keyed on the target alone -- the browser hierarchy the
 * user navigates, species then platform. So the fetch address and the store
 * path are not the same string, and only the latter belongs in a path. */
/* From YAME v1.50 the fetch address is the store path, and a directory
 * address is that directory's own files: `yame fetch -y EPICv2` is the six
 * platform files, `yame fetch -y hg38` the genome directory (which also
 * carries the suite's coordinate stream, ~30 MB sesame never opens). The
 * pre-1.50 `InfiniumAnnotation/EPICv2` spelling resolves as an alias for one
 * release. -y is not optional: without it yame refuses off a terminal, which
 * is exactly where a printed fix has to work. */

static int is_file(const char *p)
{
    return p && *p && yame_assets_is_file(p);
}

const sesame_reg_t *sesame__reg_for_platform(const char *platform)
{
    if (!platform) return NULL;
    for (const sesame_reg_t *r = SESAME_REGISTRY; r->platform; r++)
        if (strcmp(r->platform, platform) == 0) return r;
    return NULL;
}

/* The catalog tag this build expects for a genome build; "the pinned tag"
 * when the build is not one we carry. */
static const char *sesame__genome_tag(const char *genome)
{
    if (!genome) return "the pinned tag";
    for (const sesame_genome_reg_t *g = SESAME_GENOMES; g->genome; g++)
        if (strcmp(g->genome, genome) == 0) return g->tag;
    return "the pinned tag";
}

/* Is the store holding an EARLIER annotation release than this build pins?
 *
 * sesame does no verification of its own: yame_store_state() compares each
 * compiled digest against the store directory's own SHA256SUMS line and
 * composes the advice, so the sentence a sesame user reads is the one a yame
 * user reads, and it changes when yame changes it. The comparison is a string
 * compare over one small text file -- no hashing -- so this is cheap enough to
 * run on every lookup.
 *
 * Why it matters here: our only lineage guard is a row-count check, so an
 * ordering re-released with the SAME probe count and a different order passes
 * it and then mislabels every probe, silently. That is the failure this
 * catches.
 */
static const yame_fetch_cfg_t *sesame__fetch_cfg(void)
{
    static const yame_fetch_cfg_t cfg = {
        .files = SESAME_FILES, .n_files = SESAME_N_FILES,
        .tool = "yame",          /* the command that fixes it is yame's */
        .tool_env = NULL,        /* sesame reads $YAME_DATA_HOME, the shared one */
        /* sesame documents that it NEVER prompts -- a prompt hangs a Nextflow
         * job or a Docker build with no indication why. From YAME v1.52 the
         * resolver will ask on a terminal unless this is set, and a
         * designated initializer zero-fills it, i.e. defaults to prompting. */
        .no_prompt = 1,
    };
    return &cfg;
}

void sesame_warn_if_stale(const char *store_subdir)
{
    char advice[1024];
    if (!store_subdir) return;
    if (yame_store_state(sesame__fetch_cfg(), store_subdir, advice, sizeof advice)
            == YAME_STORE_STALE && advice[0])
        fprintf(stderr, "sesame: %s\n", advice);
}

const char *sesame_platform_from_beads(int32_t beads)
{
    for (const sesame_reg_t *r = SESAME_REGISTRY; r->platform; r++)
        if (r->beads && r->beads == beads) return r->platform;
    return NULL;
}

/* The one store every tool in the suite shares. */
const char *sesame_store_dir(char *out, size_t n)
{
    yame_assets_root(NULL, NULL, out, n);
    return out;
}

/* A --platform argument may be a NAME (EPICv2) or a PATH to an ordering file.
 *
 * YAME settled this rule for its own -R/-m and it is worth matching: an
 * existing path is used as given -- the ordinary spelling costs nothing and
 * cannot change meaning -- and anything else is a name resolved against the
 * store. So a custom or pre-release array needs no new flag; --platform simply
 * takes the ordering that describes it.
 *
 * --index is still the explicit spelling, and still separable: --platform names
 * the companion assets (mask, coord, SNP) while --index names the probe order,
 * and passing both is how a custom ordering borrows a known platform's masks. */
int sesame_platform_is_path(const char *arg)
{
    return arg && *arg && yame_assets_is_file(arg);
}

/* <store>/<platform>/<file>, else ./<file>. 0 on success.
 * Every per-platform asset resolves through here -- ordering, .cm mask, coord
 * table, SNP table, typeI_ext -- so there is exactly one place that knows the
 * layout, and it cannot drift from what `yame fetch` writes. */
int sesame_asset_locate(const char *platform, const char *file,
                        char *out, size_t n)
{
    char dir[4096];

    if (!platform || !file) return -1;
    sesame_store_dir(dir, sizeof dir);
    snprintf(out, n, "%s/%s/%s", dir, platform, file);
    if (is_file(out)) { sesame_warn_if_stale(platform); return 0; }

    snprintf(out, n, "./%s", file);   /* cwd convenience */
    if (is_file(out)) return 0;

    out[0] = '\0';
    return -1;
}

/* Directory holding a platform's assets, whether or not it exists yet. Callers
 * that scan it (the .cm mask lookup) need the directory, not one file. */
void sesame_asset_dir(const char *platform, char *out, size_t n)
{
    char dir[4096];
    sesame_store_dir(dir, sizeof dir);
    snprintf(out, n, "%s/%s", dir, platform ? platform : "");
}

int sesame_index_locate(const char *platform, char *out, size_t n)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    char file[512];

    if (!platform) return -1;
    snprintf(file, sizeof file, "%s",
             reg && reg->ordering ? reg->ordering : "");
    if (!file[0]) snprintf(file, sizeof file, "%s.ordering.tsv.gz", platform);
    return sesame_asset_locate(platform, file, out, n);
}

int sesame_genome_locate(const char *genome, const char *file, char *out, size_t n)
{
    char dir[4096];
    if (!genome || !file) return -1;
    sesame_store_dir(dir, sizeof dir);
    /* <store>/<genome>/<file> -- the same species-keyed slot the browser
     * shows. The genomes files sit at its root; kycg's sequencing sets sit
     * under KYCG/ in the same directory, from a different upstream. */
    snprintf(out, n, "%s/%s/%s", dir, genome, file);
    if (is_file(out)) return 0;

    snprintf(out, n, "./%s", file);   /* cwd convenience */
    if (is_file(out)) return 0;

    out[0] = '\0';
    return -1;
}

/* Identify a positional file's platform from its ROW COUNT, and return the
 * ordering that goes with it.
 *
 * A .cg is positional: row i is probe i of some ordering, and nothing in the
 * file says which. But each platform's ordering has a distinct length (EPIC
 * 866553, EPICv2 937690, HM450 486427, MSA 284309, ...), so the count alone
 * identifies it. YAME keeps that row-space table for its own -R/-m inference;
 * this asks the same question, and like YAME it RESOLVES ONLY -- a command
 * that finds the ordering missing should say which fetch supplies it, not
 * start a multi-megabyte transfer in the middle of printing a table.
 *
 * Returns 0 and fills out with the ordering path; 1 if the platform is
 * identified but not in the store (platform/fetch are set, out is not); -1 if
 * the count matches nothing -- which is the honest answer for a --collapse'd
 * beta.cg, whose rows were averaged to cg-prefixes and no longer line up with
 * any ordering. Text input returns -1 too: there is no header to read, and
 * those files carry the platform in their name anyway. */
int sesame_index_from_rows(const char *file, char *out, size_t n,
                           const char **platform, const char **fetch)
{
    uint64_t rows;
    int rc;

    if (platform) *platform = NULL;
    if (fetch) *fetch = NULL;
    if (out && n) out[0] = '\0';
    if (!file) return -1;

    rows = yame_ref_file_rows(file);
    if (!rows) return -1;                 /* not a YAME file, or unreadable */

    rc = yame_ref_for_rows(rows, NULL, "array", out, n, platform, fetch);
    if (rc == YAME_REF_OK) return 0;
    if (rc == YAME_REF_MISSING) return 1;
    return -1;
}

/* Fills msg with the "here is how to fix it" text used when no index is found.
 * Deliberately verbose: this is the error a first-time user will hit. The
 * searched path must match sesame_asset_locate(). */
void sesame_index_missing_help(const char *platform, char *msg, size_t n)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    char dir[4096];
    sesame_store_dir(dir, sizeof dir);
    snprintf(msg, n,
        "no index found for platform %s\n"
        "  searched:\n"
        "    %s/%s/%s\n"
        "    ./%s.ordering.tsv.gz\n"
        "  fix, any of:\n"
        "    yame fetch -y %s   (this build expects %s)\n"
        "    export YAME_DATA_HOME=<dir>",
        platform,
        dir, platform, reg ? reg->ordering : "<platform>.ordering.tsv.gz",
        platform,
        platform, reg ? reg->tag : "the pinned tag");
}

/* Same, for a per-platform companion asset (coord table, SNP table, .cm mask):
 * they all arrive in one `yame fetch` of the platform. */
void sesame_asset_missing_help(const char *platform, const char *file,
                               char *msg, size_t n)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    char dir[4096];
    sesame_store_dir(dir, sizeof dir);
    snprintf(msg, n,
        "no %s in the store for %s\n"
        "  searched: %s/%s/\n"
        "  fix: yame fetch -y %s   (this build expects %s)",
        file ? file : "asset", platform, dir, platform,
        platform, reg ? reg->tag : "the pinned tag");
}

/* Genome-level annotation (seqinfo/gaps/cytoband/genes) from zhou-lab/genomes,
 * which carries its own tag. */
void sesame_genome_missing_help(const char *genome, char *msg, size_t n)
{
    char dir[4096];
    sesame_store_dir(dir, sizeof dir);
    snprintf(msg, n,
        "no genome annotation for %s in the store\n"
        "  searched: %s/%s/\n"
        "  fix: yame fetch -y %s   (this build expects %s)",
        genome, dir, genome, genome, sesame__genome_tag(genome));
}
