/* mask.c -- quality masking (the Q step), reading YAME .cm mask files.
 *
 * The recommended mask sets live in a YAME-packed .cm (a family of format-0
 * bit-tracks, one bit per ordering row). We link YAME directly (both projects
 * are AGPL-3.0, same author) and read the .cm in-process: no `yame` binary at
 * runtime, no fork/exec.
 *
 * A probe is quality-masked iff it is set in ANY recommended track -- the union
 * of recommendedMaskNames(platform) (R/mask.R:197-220).
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present The Children's Hospital of Philadelphia
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "sesame.h"
#include "internal.h"

/* YAME (linked): cdata/cfile I/O, the sample index, the fmt0 bit accessor, and
 * prepare_mask (normalizes a mask block to fmt0). */
#include "cfile.h"    /* open_cfile, read_cdata1, cdata_t, cfile_t, bgzf_close */
#include "index.h"    /* loadSampleNamesFromIndex, snames_t, cleanSampleNames2 */
#include "summary.h"  /* prepare_mask */
#include "cdata.h"    /* FMT0_IN_SET, free_cdata */

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* recommendedMaskNames() (R/mask.R:197-220) -- the Q set; NULL-terminated. */
static const char *const REC_EPIC[] = {
    "mapping", "channel_switch", "snp5_GMAF1p", "extension", "sub30_copy", NULL };
static const char *const REC_MSA[] = {   /* also EPICv2 */
    "M_1baseSwitchSNPcommon_5pt", "M_2extBase_SNPcommon_5pt",
    "M_mapping", "M_nonuniq", "M_SNPcommon_5pt", NULL };

/* backgroundMask() (R/mask.R:170-185) -- the set EXCLUDED from background
 * estimation in pOOBAH/noob. Platform-agnostic name list, intersected with the
 * tracks the .cm actually has (a track not present simply never matches). */
static const char *const BG_NAMES[] = {
    "M_nonuniq", "nonunique", "sub35_copy", "multi", "design_issue",
    "M_1baseSwitchSNPcommon_5pt", "M_1baseSwitchSNPcommon_1pt", NULL };

static const char *const *recommended_names(const char *platform)
{
    if (!platform) return NULL;
    if (strcmp(platform, "EPIC") == 0 || strcmp(platform, "HM450") == 0)
        return REC_EPIC;
    if (strcmp(platform, "MSA") == 0 || strcmp(platform, "EPICv2") == 0)
        return REC_MSA;
    return NULL;
}

static int is_recommended(const char *name, const char *const *rec)
{
    for (; *rec; rec++) if (strcmp(name, *rec) == 0) return 1;
    return 0;
}

/* Locate the platform's .cm mask in the store: a file in the platform's asset
 * directory ending in ".cm" (not ".cm.idx"). Scanned rather than named because
 * the file carries its genome build (MM285.mm10.mask.cm), which the caller does
 * not have here.
 *
 * A platform may ship more than one -- MM285 has both mm10 and mm39 -- so take
 * the lexicographically first rather than whatever readdir happens to return.
 * Arbitrary, but DETERMINISTIC, which the previous first-match was not: readdir
 * order is filesystem-dependent, so the same command could silently mask
 * against a different genome build on another machine. It also lands on mm10
 * for MM285, which is the one with the full track set (110 vs 3). Name --mask
 * explicitly when the choice matters. 0 on success. */
static int find_cm(const char *platform, char *out, size_t n)
{
    char pdir[4096], best[512];
    DIR *d;
    struct dirent *de;
    int found = 0;

    sesame_asset_dir(platform, pdir, sizeof pdir);
    if (!(d = opendir(pdir))) return -1;
    best[0] = '\0';
    while ((de = readdir(d))) {
        size_t L = strlen(de->d_name);
        if (L > 3 && strcmp(de->d_name + L - 3, ".cm") == 0)
            if (!found || strcmp(de->d_name, best) < 0) {
                snprintf(best, sizeof best, "%s", de->d_name);
                found = 1;
            }
    }
    closedir(d);
    if (found) snprintf(out, n, "%s/%s", pdir, best);
    return found ? 0 : -1;
}

/* Read the platform's .cm and return the union of the named tracks as a 0/1
 * vector aligned to the ordering (out[i] == 1 -> probe i is in some track).
 * Tracks not present in the .cm simply never match. */
static int mask_union(const char *platform, const char *maskpath,
                      const char *const *names,
                      uint8_t **out, int32_t *out_n, sesame_err_t *err)
{
    char cm[4096];
    cfile_t cf;
    snames_t sn;
    uint8_t *mask = NULL;
    int32_t n = 0;
    int k = 0;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    /* No recommended track list means we do not know this platform -- a custom
     * array. An explicit --mask covers that, but only if it is unambiguous:
     * see the single-track check below. */
    if (!names && !(maskpath && *maskpath))
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "no mask set for platform '%s' -- pass --mask <file.cm>",
            platform ? platform : "(none)");
    if (maskpath && *maskpath) {
        snprintf(cm, sizeof cm, "%s", maskpath);      /* named outright */
        /* A .cm handed over directly must say exactly what to mask. Without a
         * recommended track list there is no basis for choosing among several,
         * and unioning them all is the wrong default -- a stock EPICv2 mask has
         * 25 tracks where the recommended set is 5, so quietly taking all of
         * them masks 111k probes instead of 50k. Refuse and say so. */
        if (!names) {
            snames_t probe = loadSampleNamesFromIndex(cm);
            int ntrk = probe.n;
            cleanSampleNames2(probe);
            if (ntrk > 1)
                return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
                    "--mask needs a single-track .cm here (%s has %d). Extract "
                    "the track you want, or name a published platform so its "
                    "recommended set picks the tracks.", cm, ntrk);
        }
    } else if (find_cm(platform, cm, sizeof cm) != 0) {
        char help[1024];
        sesame_asset_missing_help(platform, "mask (.cm)", help, sizeof help);
        return sesame__fail(err, SESAME_ERR_IO, "%s", help);
    }

    /* All tracks come back in index order; union the ones we want. The file is
     * small, so reading every block is cheap. */
    cf = open_cfile(cm);
    sn = loadSampleNamesFromIndex(cm);
    for (;;) {
        cdata_t c = read_cdata1(&cf);
        const char *name;
        if (c.n == 0) break;                 /* EOF (methscope pattern) */
        prepare_mask(&c);                    /* normalize to fmt0 bit form */
        if (!mask) {
            n = (int32_t)c.n;                /* bit count == probe count */
            mask = (uint8_t *)calloc((size_t)n, 1);
        }
        name = (k < sn.n) ? sn.s[k] : "";
        if (mask && (!names || is_recommended(name, names))) {
            int32_t i;
            for (i = 0; i < n; i++) if (FMT0_IN_SET(c, i)) mask[i] = 1;
        }
        free_cdata(&c);
        k++;
    }
    bgzf_close(cf.fh);
    cleanSampleNames2(sn);

    if (!mask)
        return sesame__fail(err, SESAME_ERR_FORMAT, "empty .cm mask for %s", platform);
    *out = mask;
    *out_n = n;
    return SESAME_OK;
}

/* Q: the recommended quality mask (union of recommendedMaskNames). */
int sesame_quality_mask(const char *platform, const char *maskpath,
                        uint8_t **out, int32_t *out_n,
                        sesame_err_t *err)
{
    return mask_union(platform, maskpath, recommended_names(platform), out, out_n, err);
}

/* P/B: the background mask (union of backgroundMask names), the probes EXCLUDED
 * from background estimation. */
int sesame_background_mask(const char *platform, const char *maskpath,
                           uint8_t **out, int32_t *out_n,
                           sesame_err_t *err)
{
    return mask_union(platform, maskpath, BG_NAMES, out, out_n, err);
}
