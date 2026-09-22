/* mLiftOver (R/mLiftOver.R): lift beta values from one Infinium platform to
 * another by a probe-ID prefix join. Modern platforms (EPICv2, MSA) carry a
 * `_suffix` that the legacy platforms (EPIC, HM450, HM27) lack; the join strips
 * the suffix on whichever side is modern when crossing that boundary, and joins
 * on the full ID within a family. Each TARGET probe takes the beta of the first
 * source probe sharing its prefix (R's target_uniq distinct-first), or NA. The
 * result is positional to the target ordering -- a valid target-platform beta.cg.
 *
 * The empirical liftOver.* quality mappings (mapping=) and impute=TRUE are not
 * ported; this is the mapping=NULL, impute=FALSE path.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026-present The Children's Hospital of Philadelphia
 */
#include "sesame.h"
#include "internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int is_modern(const char *p)
{
    return p && (!strcmp(p, "EPICv2") || !strcmp(p, "MSA"));
}
static int is_legacy(const char *p)
{
    return p && (!strcmp(p, "EPIC") || !strcmp(p, "HM450") || !strcmp(p, "HM27"));
}

/* Compare the prefixes of two IDs. A side with strip=1 ends its prefix at the
 * first '_' (or end); strip=0 uses the whole ID. Shorter-prefix-first ordering,
 * C locale -- matches R's strsplit("_")[[1]][1] keys. */
static int pcmp(const char *a, int as, const char *b, int bs)
{
    for (;;) {
        char ca = *a, cb = *b;
        int ta = (ca == '\0') || (as && ca == '_');
        int tb = (cb == '\0') || (bs && cb == '_');
        if (ta && tb) return 0;
        if (ta) return -1;
        if (tb) return 1;
        if (ca != cb) return (unsigned char)ca < (unsigned char)cb ? -1 : 1;
        a++; b++;
    }
}

typedef struct { const char *id; int32_t idx; } lo_ent;

/* Sort source entries by (prefix, idx) so the first matching entry for any
 * prefix carries the lowest source index -- R's distinct(.keep_all) first. */
static int src_strip_g;
static int lo_cmp(const void *x, const void *y)
{
    const lo_ent *a = (const lo_ent *)x, *b = (const lo_ent *)y;
    int c = pcmp(a->id, src_strip_g, b->id, src_strip_g);
    if (c) return c;
    return a->idx < b->idx ? -1 : a->idx > b->idx ? 1 : 0;
}

int sesame_liftover_betas(const char *src_platform, const sesame_index_t *src_ix,
                          const char *tgt_platform, const sesame_index_t *tgt_ix,
                          const double *mat_in, int32_t nsamp,
                          double **mat_out, sesame_err_t *err)
{
    int32_t nsrc, ntgt, i, j;
    int src_strip = 0, tgt_strip = 0;
    lo_ent *ent = NULL;
    int32_t *map = NULL;      /* target i -> source index, or -1 */
    double *out = NULL;
    int rc = SESAME_OK;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!src_ix || !tgt_ix || !mat_in || !mat_out)
        return sesame__fail(err, SESAME_ERR_IO, "null argument");
    nsrc = sesame_index_nprobes(src_ix);
    ntgt = sesame_index_nprobes(tgt_ix);
    if (nsrc <= 0 || ntgt <= 0) return sesame__fail(err, SESAME_ERR_FORMAT, "empty index");

    /* Direction: strip the suffix on the modern side of a modern<->legacy lift;
     * otherwise join on the full ID (same family, incl. identity). */
    if (is_legacy(tgt_platform) && is_modern(src_platform)) src_strip = 1;
    else if (is_modern(tgt_platform) && is_legacy(src_platform)) tgt_strip = 1;

    ent = (lo_ent *)malloc((size_t)nsrc * sizeof *ent);
    map = (int32_t *)malloc((size_t)ntgt * sizeof *map);
    if (!ent || !map) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto out; }

    for (i = 0; i < nsrc; i++) { ent[i].id = sesame_index_probe_id(src_ix, i); ent[i].idx = i; }
    src_strip_g = src_strip;
    qsort(ent, (size_t)nsrc, sizeof *ent, lo_cmp);

    /* For each target probe, lower-bound its prefix in the sorted source keys. */
    for (i = 0; i < ntgt; i++) {
        const char *tid = sesame_index_probe_id(tgt_ix, i);
        int32_t lo = 0, hi = nsrc;             /* first ent with prefix >= tid */
        while (lo < hi) {
            int32_t mid = lo + ((hi - lo) >> 1);
            if (pcmp(ent[mid].id, src_strip, tid, tgt_strip) < 0) lo = mid + 1;
            else hi = mid;
        }
        map[i] = (lo < nsrc && pcmp(ent[lo].id, src_strip, tid, tgt_strip) == 0)
                 ? ent[lo].idx : -1;
    }

    out = (double *)malloc((size_t)nsamp * (size_t)ntgt * sizeof(double));
    if (!out) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto out; }
    for (j = 0; j < nsamp; j++) {
        const double *si = mat_in + (size_t)j * (size_t)nsrc;
        double *ti = out + (size_t)j * (size_t)ntgt;
        for (i = 0; i < ntgt; i++) ti[i] = map[i] >= 0 ? si[map[i]] : NAN;
    }

    *mat_out = out; out = NULL;

out:
    free(out); free(ent); free(map);
    return rc;
}
