/* pipeline_dump -- validation harness. Runs the library pipeline
 * (sesame_pipeline) on one sample and dumps a chosen quantity at full double
 * precision, so the differential tests keep their bit-identical / ULP gates,
 * which the float32 .cg product cannot carry. Not a product command.
 *
 *   pipeline_dump --index IDX [--platform P] --prep CODE
 *                 --what beta|total|pval|M|U|col [--f64] PREFIX
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int resolve(const char *pfx, const char *chan, char *out, size_t n)
{
    FILE *f;
    snprintf(out, n, "%s_%s.idat", pfx, chan);
    if ((f = fopen(out, "rb"))) { fclose(f); return 0; }
    snprintf(out, n, "%s_%s.idat.gz", pfx, chan);
    if ((f = fopen(out, "rb"))) { fclose(f); return 0; }
    return -1;
}

int main(int argc, char **argv)
{
    const char *idxpath = NULL, *platform = NULL, *prep = "", *what = "beta", *prefix = NULL;
    int f64 = 0, pneg = 0, i, rc = 1;
    char gp[4096], rp[4096];
    sesame_idat_t *g = NULL, *r = NULL;
    sesame_index_t *ix = NULL;
    sesame_sigdf_t *sig = NULL;
    uint8_t *qmask = NULL, *bgmask = NULL, *col = NULL;
    int32_t qn = 0, bgn = 0, n, k;
    double *out = NULL, *beta = NULL, *M = NULL, *U = NULL, *tot = NULL, *pval = NULL;
    sesame_err_t e;

    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--index") && i+1 < argc) idxpath = argv[++i];
        else if (!strcmp(argv[i], "--platform") && i+1 < argc) platform = argv[++i];
        else if (!strcmp(argv[i], "--prep") && i+1 < argc) prep = argv[++i];
        else if (!strcmp(argv[i], "--what") && i+1 < argc) what = argv[++i];
        else if (!strcmp(argv[i], "--f64")) f64 = 1;
        else if (!strcmp(argv[i], "--detection") && i+1 < argc) pneg = !strcmp(argv[++i], "pneg");
        else prefix = argv[i];
    }
    if (!prefix) {
        fprintf(stderr, "usage: pipeline_dump [--index IDX] [--platform P] --prep CODE "
                "--what beta|total|pval|M|U|col [--f64] PREFIX\n");
        return 2;
    }
    if (resolve(prefix, "Grn", gp, sizeof gp) || resolve(prefix, "Red", rp, sizeof rp)) {
        fprintf(stderr, "no IDAT for %s\n", prefix); return 1;
    }
    if (sesame_idat_read(gp, &g, &e) || sesame_idat_read(rp, &r, &e)) {
        fprintf(stderr, "%s\n", e.msg); goto out;
    }
    if (!platform) platform = sesame_platform_from_beads(g->n);
    if (!idxpath) {                              /* auto-locate in the store */
        static char resolved[4096];
        if (!platform || sesame_index_locate(platform, resolved, sizeof resolved)) {
            fprintf(stderr, "no index for %s; pass --index\n", platform ? platform : "?");
            goto out;
        }
        idxpath = resolved;
    }
    if (!(ix = sesame_index_open(idxpath, &e))) { fprintf(stderr, "%s\n", e.msg); goto out; }
    if (!(sig = sesame_sigdf_from_idats(g, r, ix, 0, &e))) { fprintf(stderr, "%s\n", e.msg); goto out; }

    if (strchr(prep, 'Q') && sesame_quality_mask(platform, NULL, &qmask, &qn, &e))
        { fprintf(stderr, "%s\n", e.msg); goto out; }
    if ((strchr(prep,'P') || strchr(prep,'B') || !strcmp(what,"pval")) &&
        sesame_background_mask(platform, NULL, &bgmask, &bgn, &e))
        { fprintf(stderr, "%s\n", e.msg); goto out; }

    n = sesame_index_nprobes(ix);
    out = (double *)malloc((size_t)n * sizeof(double));
    if (!out) { fprintf(stderr, "oom\n"); goto out; }
    if      (!strcmp(what, "beta"))  beta = out;
    else if (!strcmp(what, "total")) tot = out;
    else if (!strcmp(what, "pval"))  pval = out;
    else if (!strcmp(what, "M"))     M = out;
    else if (!strcmp(what, "U"))     U = out;
    else if (!strcmp(what, "col"))   col = (uint8_t *)malloc((size_t)n);
    else { fprintf(stderr, "unknown --what %s\n", what); goto out; }

    if (sesame_pipeline(sig, prep, qmask, qn, bgmask, bgn, 0, pneg,
                        beta, M, U, tot, pval, col, &e)) { fprintf(stderr, "%s\n", e.msg); goto out; }

    if (col) {
        static const char *nm[] = { "2", "G", "R" };
        for (k = 0; k < n; k++)
            printf("%s\t%s\n", sesame_index_probe_id(ix, k), nm[col[k]]);
    } else if (f64) {
        if (fwrite(out, sizeof(double), (size_t)n, stdout) != (size_t)n) goto out;
    } else {
        for (k = 0; k < n; k++) {
            const char *id = sesame_index_probe_id(ix, k);
            if (isnan(out[k])) printf("%s\tNA\n", id);
            else               printf("%s\t%.17g\n", id, out[k]);
        }
    }
    rc = 0;

out:
    free(out); free(col); free(qmask); free(bgmask);
    sesame_sigdf_free(sig); sesame_idat_free(g); sesame_idat_free(r);
    sesame_index_close(ix);
    return rc;
}
