/* cgwrite.c -- write per-probe values as a YAME format-4 .cg file (one float32
 * per probe, NA encoded as -1.0, RLE over NA runs), plus a companion .idx of
 * "<sample>\t<bgzf-offset>" lines. Links YAME directly (both AGPL), same as the
 * mask reader.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present The Children's Hospital of Philadelphia
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "sesame.h"
#include "internal.h"

#include "cfile.h"    /* cdata_write1, open_cfile, read_cdata1, BGZF */
#include "cdata.h"    /* cdata_t, cdata_compress, decompress, bgzf_open2/tell/close */
#include "index.h"    /* loadSampleNamesFromIndex, cleanSampleNames2, snames_t */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Write the companion "<path>.idx" of "<name>\t<bgzf-offset>" lines. */
static void write_idx(const char *path, char *const *names,
                      const int64_t *offs, int32_t nsamp)
{
    size_t z = strlen(path) + 5;
    char *idxpath = (char *)malloc(z);
    FILE *ix;
    int32_t j;
    if (!idxpath) return;
    snprintf(idxpath, z, "%s.idx", path);
    if ((ix = fopen(idxpath, "w"))) {
        for (j = 0; j < nsamp; j++)
            fprintf(ix, "%s\t%lld\n", names[j], (long long)offs[j]);
        fclose(ix);
    }
    free(idxpath);
}

/* Read a format-4 .cg into a sample-major matrix (sample j, probe i at
 * mat[j*nprobe+i]) and the sample names from <path>.idx. NA (negative float) ->
 * NaN. Caller frees mat, names[i], and names. */
int sesame_read_cg(const char *path, double **mat_out, int32_t *nprobe_out,
                   int32_t *nsamp_out, char ***names_out, sesame_err_t *err)
{
    cfile_t cf;
    snames_t sn;
    double *mat = NULL;
    char **names = NULL, pathbuf[4096];
    int32_t np = -1, ns = 0, cap = 0, i;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    snprintf(pathbuf, sizeof pathbuf, "%s", path);
    cf = open_cfile(pathbuf);
    if (!cf.fh) return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", path);
    sn = loadSampleNamesFromIndex(pathbuf);

    for (;;) {
        cdata_t c = read_cdata1(&cf), d;
        float *s;
        if (c.n == 0) break;                     /* EOF */
        d = decompress(c);
        if (d.fmt != '4') {
            free_cdata(&c); free_cdata(&d); bgzf_close(cf.fh); cleanSampleNames2(sn);
            free(mat); free(names);
            return sesame__fail(err, SESAME_ERR_FORMAT, "%s is format %c, expected 4", path, d.fmt);
        }
        if (np < 0) np = (int32_t)d.n;
        else if ((int32_t)d.n != np) {
            free_cdata(&c); free_cdata(&d); bgzf_close(cf.fh); cleanSampleNames2(sn);
            free(mat); free(names);
            return sesame__fail(err, SESAME_ERR_FORMAT, "inconsistent probe count in %s", path);
        }
        if (ns >= cap) {
            cap = cap ? cap * 2 : 8;
            mat = (double *)realloc(mat, (size_t)cap * (size_t)np * sizeof(double));
            names = (char **)realloc(names, (size_t)cap * sizeof(char *));
        }
        s = (float *)d.s;
        for (i = 0; i < np; i++) {
            float v = s[i];
            mat[(size_t)ns * (size_t)np + (size_t)i] = v < 0.0f ? NAN : (double)v;
        }
        names[ns] = strdup(ns < sn.n ? sn.s[ns] : "");
        ns++;
        free_cdata(&c); free_cdata(&d);
    }
    bgzf_close(cf.fh);
    cleanSampleNames2(sn);
    *mat_out = mat; *nprobe_out = np; *nsamp_out = ns; *names_out = names;
    return SESAME_OK;
}

/* Read a .cg's per-probe TOTAL intensity into a sample-major matrix, regardless
 * of format: format 4 gives the stored float (a total_intensity.cg), format 3
 * gives M+U (a cnvnormals.cg). NA (negative float, or M==U==0) -> NaN. This is
 * what CNV consumes for both the target and the normal reference. Caller frees
 * mat, names[i], and names. */
int sesame_read_cg_total(const char *path, double **mat_out, int32_t *nprobe_out,
                         int32_t *nsamp_out, char ***names_out, sesame_err_t *err)
{
    cfile_t cf;
    snames_t sn;
    double *mat = NULL;
    char **names = NULL, pathbuf[4096];
    int32_t np = -1, ns = 0, cap = 0, i;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    snprintf(pathbuf, sizeof pathbuf, "%s", path);
    cf = open_cfile(pathbuf);
    if (!cf.fh) return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", path);
    sn = loadSampleNamesFromIndex(pathbuf);

    for (;;) {
        cdata_t c = read_cdata1(&cf), d;
        if (c.n == 0) break;                         /* EOF */
        d = decompress(c);
        if (d.fmt != '3' && d.fmt != '4') {
            free_cdata(&c); free_cdata(&d); bgzf_close(cf.fh); cleanSampleNames2(sn);
            free(mat); free(names);
            return sesame__fail(err, SESAME_ERR_FORMAT,
                "%s is format %c, expected 3 or 4", path, d.fmt);
        }
        if (np < 0) np = (int32_t)d.n;
        else if ((int32_t)d.n != np) {
            free_cdata(&c); free_cdata(&d); bgzf_close(cf.fh); cleanSampleNames2(sn);
            free(mat); free(names);
            return sesame__fail(err, SESAME_ERR_FORMAT, "inconsistent probe count in %s", path);
        }
        if (ns >= cap) {
            cap = cap ? cap * 2 : 8;
            mat = (double *)realloc(mat, (size_t)cap * (size_t)np * sizeof(double));
            names = (char **)realloc(names, (size_t)cap * sizeof(char *));
        }
        for (i = 0; i < np; i++) {
            double v;
            if (d.fmt == '4') {
                float f = ((float *)d.s)[i];
                v = f < 0.0f ? NAN : (double)f;
            } else {
                uint64_t mu = f3_get_mu(&d, (uint64_t)i);
                uint64_t M = mu >> 32, U = (mu << 32) >> 32;
                v = (M == 0 && U == 0) ? NAN : (double)(M + U);
            }
            mat[(size_t)ns * (size_t)np + (size_t)i] = v;
        }
        names[ns] = strdup(ns < sn.n ? sn.s[ns] : "");
        ns++;
        free_cdata(&c); free_cdata(&d);
    }
    bgzf_close(cf.fh);
    cleanSampleNames2(sn);
    *mat_out = mat; *nprobe_out = np; *nsamp_out = ns; *names_out = names;
    return SESAME_OK;
}

/* Format 3 (M/U counts), the default. Rounds M,U to integers -- exact for raw
 * IDAT signal (uint16 means), sub-integer rounding for preprocessed floats
 * (noise). NA (either allele NaN) is stored as (0,0), which yame reads as
 * beta=NA / cov=0. yame derives beta (MU2beta) and total intensity (MU2cov).
 * matM/matU are sample-major: sample j, probe i at mat[j*nprobe+i]. */
int sesame_write_cg_mu(const char *path, const double *matM, const double *matU,
                       int32_t nprobe, int32_t nsamp, char *const *names,
                       sesame_err_t *err)
{
    BGZF *fp;
    int64_t *offs;
    int32_t j, i;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!(fp = bgzf_open2(path, "w")))
        return sesame__fail(err, SESAME_ERR_IO, "cannot open %s for writing", path);
    if (!(offs = (int64_t *)malloc((size_t)nsamp * sizeof(int64_t)))) {
        bgzf_close(fp); return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
    }
    for (j = 0; j < nsamp; j++) {
        cdata_t c;
        memset(&c, 0, sizeof c);
        c.fmt = '3'; c.compressed = 0; c.n = (uint64_t)nprobe; c.unit = 8;
        c.s = (uint8_t *)calloc((size_t)nprobe, 8);
        if (!c.s) { free(offs); bgzf_close(fp);
                    return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }
        for (i = 0; i < nprobe; i++) {
            double m = matM[(size_t)j*(size_t)nprobe + (size_t)i];
            double u = matU[(size_t)j*(size_t)nprobe + (size_t)i];
            if (isnan(m) || isnan(u)) f3_set_mu(&c, (uint64_t)i, 0, 0);   /* NA */
            else f3_set_mu(&c, (uint64_t)i, (uint64_t)llround(m), (uint64_t)llround(u));
        }
        offs[j] = bgzf_tell(fp);
        cdata_compress(&c);
        cdata_write1(fp, &c);
        free(c.s);
    }
    bgzf_close(fp);
    write_idx(path, names, offs, nsamp);
    free(offs);
    return SESAME_OK;
}

/* Format 4 (one float32 per probe, NA=-1.0), on request -- exact values (e.g.
 * betas or total intensity) with no rounding. */
int sesame_write_cg(const char *path, const double *mat, int32_t nprobe,
                    int32_t nsamp, char *const *names, sesame_err_t *err)
{
    BGZF *fp;
    int64_t *offs;
    int32_t j, i;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!(fp = bgzf_open2(path, "w")))
        return sesame__fail(err, SESAME_ERR_IO, "cannot open %s for writing", path);
    if (!(offs = (int64_t *)malloc((size_t)nsamp * sizeof(int64_t)))) {
        bgzf_close(fp); return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
    }
    for (j = 0; j < nsamp; j++) {
        cdata_t c;
        float *s = (float *)malloc((size_t)nprobe * sizeof(float));
        if (!s) { free(offs); bgzf_close(fp);
                  return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }
        for (i = 0; i < nprobe; i++) {
            double v = mat[(size_t)j * (size_t)nprobe + (size_t)i];
            s[i] = isnan(v) ? -1.0f : (float)v;     /* NA -> negative (bit31 set) */
        }
        memset(&c, 0, sizeof c);
        c.fmt = '4'; c.compressed = 0; c.n = (uint64_t)nprobe;
        c.unit = sizeof(float); c.s = (uint8_t *)s;
        offs[j] = bgzf_tell(fp);
        cdata_compress(&c);                         /* frees s, sets c.s compressed */
        cdata_write1(fp, &c);
        free(c.s);
    }
    bgzf_close(fp);
    write_idx(path, names, offs, nsamp);
    free(offs);
    return SESAME_OK;
}
