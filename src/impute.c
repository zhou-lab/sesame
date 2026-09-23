/* imputeBetas* (R/impute.R): fill missing (NaN) beta values.
 *
 *  - mean (imputeBetasMatrixByMean): replace each NaN with the mean of its
 *    probe across samples (axis "probe", R axis=1) or of its sample across
 *    probes (axis "sample", R axis=2).
 *  - neighbors (imputeBetasByGenomicNeighbors): replace each missing probe's
 *    NaN with the mean of up to `max_neighbors` nearest non-missing probes
 *    within `max_dist` bp, searched in the 5'->3' direction of the probe's
 *    strand (R resizes each missing range start-anchored + strand-aware, then
 *    findOverlaps -- strand-specific -- against the non-missing probes). Ties at
 *    the neighbor cutoff are all kept (slice_min with_ties). Uses the store's
 *    coord table, which equals sesameData's manifest coords for every mapQ>=1
 *    probe, so mapped probes match R; the two disagree only on mapQ=0 (multi-
 *    mapping) probes. Unmapped probes are LEFT NA rather than copy R's artifact
 *    of placing them all at a shared 0 position and averaging. imputeBetas (the
 *    celltype-reference median path) needs a shipped reference panel; not ported.
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
#include <zlib.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int sesame_impute_mean(double *mat, int32_t nprobe, int32_t nsamp, int axis,
                       sesame_err_t *err)
{
    int32_t i, j;
    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!mat) return sesame__fail(err, SESAME_ERR_IO, "null matrix");

    size_t np = (size_t)nprobe;
    if (axis == 1) {                 /* per probe, mean across samples */
        for (i = 0; i < nprobe; i++) {
            double sum = 0.0; int32_t cnt = 0, k;
            for (j = 0; j < nsamp; j++) {
                double v = mat[(size_t)j*np + (size_t)i];
                if (!isnan(v)) { sum += v; cnt++; }
            }
            if (!cnt) continue;      /* all NaN -> mean NaN, leave as is */
            { double m = sum / cnt;
              for (k = 0; k < nsamp; k++)
                  if (isnan(mat[(size_t)k*np + (size_t)i])) mat[(size_t)k*np + (size_t)i] = m; }
        }
    } else if (axis == 2) {          /* per sample, mean across probes */
        for (j = 0; j < nsamp; j++) {
            double *col = mat + (size_t)j*np;
            double sum = 0.0; int32_t cnt = 0;
            for (i = 0; i < nprobe; i++) if (!isnan(col[i])) { sum += col[i]; cnt++; }
            if (!cnt) continue;
            { double m = sum / cnt;
              for (i = 0; i < nprobe; i++) if (isnan(col[i])) col[i] = m; }
        }
    } else {
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED, "axis must be 1 (probe) or 2 (sample)");
    }
    return SESAME_OK;
}

/* ------------------------------------------------------------- neighbors --- */

typedef struct { int32_t cid; long beg; int32_t idx; } coord_ent;

static int coord_cmp(const void *x, const void *y)
{
    const coord_ent *a = x, *b = y;
    if (a->cid != b->cid) return a->cid < b->cid ? -1 : 1;
    if (a->beg != b->beg) return a->beg < b->beg ? -1 : 1;
    return a->idx < b->idx ? -1 : a->idx > b->idx ? 1 : 0;
}

/* Load the probe coordinate table (CpG_chrm, CpG_beg, strand, mapQ -- the same
 * <plat>.<genome>.coord.tsv.gz `region` uses), positional in the ordering. cid[i]
 * is a small chromosome id (-1 if unmapped), beg[i]/end[i] the 1-based inclusive
 * CpG range (CpG_beg is 0-based, so beg = CpG_beg+1, end = CpG_beg+2 -- the
 * range GRanges reports), strnd[i] the strand code (0 '-', 1 '+', 2 other).
 *
 * This is the confidently-mapped position for every mapQ>=1 probe -- identical
 * to sesameData_getManifestGRanges there, so neighbour selection matches R. The
 * two annotations disagree only on mapQ=0 (multi-/non-mapping) probes, which
 * impute to noise either way. */
/* strand code: 0 = '-', 1 = '+', 2 = other/unstranded. */
static int load_coords_strand(const char *path, int32_t np, int32_t *cid,
                              long *beg, long *end, uint8_t *strnd, sesame_err_t *err)
{
    gzFile f = gzopen(path, "rb");
    char line[512];
    char (*names)[64] = NULL; int32_t ncap = 0, nchrom = 0, row = 0;

    if (!f) return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", path);
    if (!gzgets(f, line, sizeof line)) { gzclose(f); return sesame__fail(err, SESAME_ERR_FORMAT, "empty %s", path); }
    while (gzgets(f, line, sizeof line)) {
        char *tok[4] = {0}, *p = line; int nt = 0;
        if (row >= np) { row++; continue; }
        tok[nt++] = p;
        while (*p && nt < 4) { if (*p == '\t') { *p = '\0'; tok[nt++] = p+1; } p++; }
        if (nt < 3 || tok[0][0]=='\0' || !strcmp(tok[0],"NA") || !strcmp(tok[0],"*")) {
            cid[row]=-1; beg[row]=-1; end[row]=-1; strnd[row]=2;
        } else {
            int32_t k, found = -1;
            long cb = strtol(tok[1], NULL, 10);        /* CpG_beg, 0-based */
            for (k = 0; k < nchrom; k++) if (!strcmp(names[k], tok[0])) { found = k; break; }
            if (found < 0) {
                if (nchrom == ncap) {
                    int32_t nc = ncap ? ncap*2 : 32;
                    void *q = realloc(names, (size_t)nc*sizeof(*names));
                    if (!q) { free(names); gzclose(f); return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }
                    names = q; ncap = nc;
                }
                snprintf(names[nchrom], 64, "%s", tok[0]); found = nchrom++;
            }
            cid[row] = found;
            beg[row] = cb + 1;                          /* 1-based start */
            end[row] = cb + 2;                          /* CpG spans two bases */
            strnd[row] = (tok[2][0]=='-') ? 0 : (tok[2][0]=='+') ? 1 : 2;
        }
        row++;
    }
    free(names); gzclose(f);
    if (row != np)
        return sesame__fail(err, SESAME_ERR_FORMAT,
            "%s has %d rows, ordering has %d -- lineage mismatch", path, row, np);
    return SESAME_OK;
}

int sesame_impute_neighbors(double *mat, int32_t nprobe, int32_t nsamp,
                            const char *coords_path, int max_neighbors,
                            long max_dist, sesame_err_t *err)
{
    int32_t *cid = NULL; long *beg = NULL, *end = NULL; uint8_t *strnd = NULL;
    coord_ent *srt = NULL;
    double *dbuf = NULL, *impv = NULL; int32_t *impi = NULL;
    int rc = SESAME_OK;
    int32_t i, j;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!mat || !coords_path) return sesame__fail(err, SESAME_ERR_IO, "null argument");
    if (max_neighbors < 1) max_neighbors = 1;

    cid   = malloc((size_t)nprobe*sizeof *cid);
    beg   = malloc((size_t)nprobe*sizeof *beg);
    end   = malloc((size_t)nprobe*sizeof *end);
    strnd = malloc((size_t)nprobe*sizeof *strnd);
    srt   = malloc((size_t)nprobe*sizeof *srt);
    dbuf  = malloc((size_t)nprobe*sizeof *dbuf);     /* per-missing neighbor dists */
    impv  = malloc((size_t)nprobe*sizeof *impv);
    impi  = malloc((size_t)nprobe*sizeof *impi);
    if (!cid||!beg||!end||!strnd||!srt||!dbuf||!impv||!impi) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto out; }

    if ((rc = load_coords_strand(coords_path, nprobe, cid, beg, end, strnd, err)) != SESAME_OK) goto out;

    /* probes sorted by (chrom, beg); mapped probes only (cid>=0) come first. */
    for (i = 0; i < nprobe; i++) { srt[i].cid = cid[i]; srt[i].beg = beg[i]; srt[i].idx = i; }
    qsort(srt, (size_t)nprobe, sizeof *srt, coord_cmp);

    for (j = 0; j < nsamp; j++) {
        double *col = mat + (size_t)j*(size_t)nprobe;
        int32_t nimp = 0;
        for (i = 0; i < nprobe; i++) {
            long bi, ei, ws, we; int32_t lo, hi, s, e, k;
            double nb[512]; int32_t nn = 0;
            if (!isnan(col[i]) || cid[i] < 0) continue;   /* only missing, mapped */
            bi = beg[i]; ei = end[i];                       /* mapped 1-based range */
            if (strnd[i] != 0) { ws = bi; we = bi + max_dist - 1; }  /* + / *: start */
            else               { we = ei; ws = ei - max_dist + 1; }  /* -: end-anchor */

            /* binary-search this chrom's slice in srt for beg in [ws-1, we]. */
            lo = 0; hi = nprobe;                            /* first srt >= (cid[i], ws-1) */
            while (lo < hi) { int32_t m = (lo+hi)>>1;
                if (srt[m].cid < cid[i] || (srt[m].cid==cid[i] && srt[m].beg < ws-1)) lo = m+1; else hi = m; }
            s = lo;
            lo = s; hi = nprobe;                            /* first srt > (cid[i], we) */
            while (lo < hi) { int32_t m = (lo+hi)>>1;
                if (srt[m].cid < cid[i] || (srt[m].cid==cid[i] && srt[m].beg <= we)) lo = m+1; else hi = m; }
            e = lo;
            for (k = s; k < e; k++) {
                int32_t n = srt[k].idx; long bn, en; double d1, d2, d;
                if (n == i || isnan(col[n])) continue;      /* non-missing only */
                /* strand-aware overlap: same strand, or either side unstranded (*) */
                if (strnd[n] != 2 && strnd[i] != 2 && strnd[n] != strnd[i]) continue;
                bn = beg[n]; en = end[n];
                if (!(bn <= we && en >= ws)) continue;      /* overlap the window */
                d1 = (double)(bi - en - 1); d2 = (double)(bn - ei - 1);
                d = d1 > d2 ? d1 : d2;
                if (nn < (int32_t)(sizeof nb/sizeof nb[0])) { dbuf[nn] = d; nb[nn] = col[n]; nn++; }
            }
            if (nn == 0) continue;
            /* slice_min(n=max_neighbors, with_ties): keep dist <= the k-th smallest. */
            {
                double thr, tmp[512]; int32_t t;
                double sum = 0.0; int32_t cnt = 0;
                for (t = 0; t < nn; t++) tmp[t] = dbuf[t];
                sesame__sort(tmp, nn);
                thr = tmp[nn <= max_neighbors ? nn-1 : max_neighbors-1];
                for (t = 0; t < nn; t++) if (dbuf[t] <= thr) { sum += nb[t]; cnt++; }
                impi[nimp] = i; impv[nimp] = sum / cnt; nimp++;
            }
        }
        for (i = 0; i < nimp; i++) col[impi[i]] = impv[i];   /* apply after the pass */
    }

out:
    free(cid); free(beg); free(end); free(strnd); free(srt); free(dbuf); free(impv); free(impi);
    return rc;
}
