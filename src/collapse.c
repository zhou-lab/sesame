/* betasCollapseToPfx (R/sesame.R:157-170): average beta values that share a
 * probe-ID prefix (the part before the first '_'), collapsing EPICv2/MSA
 * replicate probes (cg00000029_TC21, ...) to one value per prefix.
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
#include "sesame.h"
#include "internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Compare two probe IDs by their prefix (everything up to the first '_', or the
 * whole string if none). '_' and '\0' both terminate a prefix; a shorter prefix
 * that is a proper head of the other sorts first. Matches R's split() ordering
 * of the sorted unique prefixes under the C locale. */
static int pfx_cmp(const char *a, const char *b)
{
    for (;;) {
        char ca = *a, cb = *b;
        int ta = (ca == '\0' || ca == '_'), tb = (cb == '\0' || cb == '_');
        if (ta && tb) return 0;
        if (ta) return -1;
        if (tb) return 1;
        if (ca != cb) return (unsigned char)ca < (unsigned char)cb ? -1 : 1;
        a++; b++;
    }
}

typedef struct { const char *id; int32_t orig; } pfx_ent;

static int pfx_ent_cmp(const void *x, const void *y)
{
    const pfx_ent *a = (const pfx_ent *)x, *b = (const pfx_ent *)y;
    int c = pfx_cmp(a->id, b->id);
    if (c) return c;
    /* stable within a group by original position (irrelevant to the mean, but
     * keeps the walk deterministic). */
    return a->orig < b->orig ? -1 : a->orig > b->orig ? 1 : 0;
}

/* Duplicate the prefix of `id` (up to the first '_' or end) into a fresh string. */
static char *pfx_dup(const char *id)
{
    size_t k = 0;
    char *out;
    while (id[k] && id[k] != '_') k++;
    out = (char *)malloc(k + 1);
    if (!out) return NULL;
    memcpy(out, id, k);
    out[k] = '\0';
    return out;
}

int sesame_betas_collapse_prefix(const sesame_index_t *ix,
                                 const double *mat_in, int32_t nsamp,
                                 double **mat_out, char ***pfx_out, int32_t *m_out,
                                 sesame_err_t *err)
{
    int32_t n, i, j, m = 0;
    pfx_ent *ent = NULL;
    int32_t *grp = NULL;       /* probe i -> group id (0..m-1) */
    char **pfx = NULL;
    double *out = NULL, *sum = NULL;
    int32_t *cnt = NULL;
    int rc = SESAME_OK;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!ix || !mat_in || !mat_out || !pfx_out || !m_out)
        return sesame__fail(err, SESAME_ERR_IO, "null argument");
    n = sesame_index_nprobes(ix);
    if (n <= 0) return sesame__fail(err, SESAME_ERR_FORMAT, "empty index");

    ent = (pfx_ent *)malloc((size_t)n * sizeof *ent);
    grp = (int32_t *)malloc((size_t)n * sizeof *grp);
    pfx = (char **)malloc((size_t)n * sizeof *pfx);   /* <= n groups */
    if (!ent || !grp || !pfx) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto out; }

    for (i = 0; i < n; i++) { ent[i].id = sesame_index_probe_id(ix, i); ent[i].orig = i; }
    qsort(ent, (size_t)n, sizeof *ent, pfx_ent_cmp);

    /* Assign each probe a group id in sorted-prefix order; collect unique pfxes. */
    for (i = 0; i < n; i++) {
        if (i == 0 || pfx_cmp(ent[i].id, ent[i-1].id) != 0) {
            if (!(pfx[m] = pfx_dup(ent[i].id))) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto out; }
            m++;
        }
        grp[ent[i].orig] = m - 1;
    }

    out = (double *)malloc((size_t)nsamp * (size_t)m * sizeof(double));
    sum = (double *)malloc((size_t)m * sizeof(double));
    cnt = (int32_t *)malloc((size_t)m * sizeof(int32_t));
    if (!out || !sum || !cnt) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto out; }

    for (j = 0; j < nsamp; j++) {
        const double *col = mat_in + (size_t)j * (size_t)n;
        int32_t g;
        for (g = 0; g < m; g++) { sum[g] = 0.0; cnt[g] = 0; }
        for (i = 0; i < n; i++) {
            double b = col[i];
            if (!isnan(b)) { sum[grp[i]] += b; cnt[grp[i]]++; }
        }
        for (g = 0; g < m; g++)
            out[(size_t)j * (size_t)m + (size_t)g] = cnt[g] ? sum[g] / (double)cnt[g] : NAN;
    }

    *mat_out = out; out = NULL;
    *pfx_out = pfx; pfx = NULL;
    *m_out = m;

out:
    if (pfx) { for (i = 0; i < m; i++) free(pfx[i]); free(pfx); }
    free(out); free(sum); free(cnt); free(grp); free(ent);
    return rc;
}
