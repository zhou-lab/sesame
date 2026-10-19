/* Titration calibration: beta -> modification fraction.
 *
 * Each probe's response on the methylation standards (0..100 % mixtures) was
 * fitted with the two-bead curve
 *     y(f) = floor + (top - floor) * rho f / (rho f + 1 - f)
 * and shipped per platform as <PLAT>.titration.tsv.gz (InfiniumAnnotation;
 * columns floor, top, rho; positional to the ordering; NA where a standard was
 * missing). A sample's response uses the nonlinearity compressed by one global
 * exponent, rho_eff = rho^gamma, and inverts in closed form:
 *     b' = (beta - floor) / amp,   amp = max(top - floor, 0.05)
 *     f  = b' / (b' + rho^gamma (1 - b')),   clamped to [0, 1]
 * gamma = 0.15 was fitted on EPICv2 bACE arrays against ACE-seq
 * (labjournal 2024/20240318_titration_gamma.R); the amp floor of 0.05 is the
 * reference fit's, so flat probes invert the same way here as there.
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
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double parse_num(const char *s, char **end)
{
    if (strncmp(s, "NA", 2) == 0) { *end = (char *)s + 2; return NAN; }
    return strtod(s, end);
}

int sesame_titration_load(const char *path, int32_t nprobe,
                          sesame_titration_t *out, sesame_err_t *err)
{
    gzFile f;
    char line[512];
    int32_t row = 0;
    sesame_titration_t t = { NULL, NULL, NULL, 0 };

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!(f = gzopen(path, "rb")))
        return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", path);
    if (!gzgets(f, line, sizeof line) || strncmp(line, "floor\ttop\trho", 13) != 0) {
        gzclose(f);
        return sesame__fail(err, SESAME_ERR_FORMAT,
            "%s: header must be floor<TAB>top<TAB>rho", path);
    }
    t.floor = (double *)malloc((size_t)nprobe * sizeof(double));
    t.amp   = (double *)malloc((size_t)nprobe * sizeof(double));
    t.rho   = (double *)malloc((size_t)nprobe * sizeof(double));
    if (!t.floor || !t.amp || !t.rho) {
        gzclose(f); sesame_titration_free(&t);
        return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
    }
    while (gzgets(f, line, sizeof line)) {
        char *p = line, *e;
        double fl, tp, rh;
        fl = parse_num(p, &e); if (*e != '\t') goto bad; p = e + 1;
        tp = parse_num(p, &e); if (*e != '\t') goto bad; p = e + 1;
        rh = parse_num(p, &e);
        if (row < nprobe) {
            if (isnan(fl) || isnan(tp) || isnan(rh)) {
                t.floor[row] = t.amp[row] = t.rho[row] = NAN;
            } else {
                t.floor[row] = fl;
                t.amp[row] = (tp - fl > 0.05) ? tp - fl : 0.05;
                t.rho[row] = rh;
            }
        }
        row++;
        continue;
    bad:
        gzclose(f); sesame_titration_free(&t);
        return sesame__fail(err, SESAME_ERR_FORMAT,
            "%s: row %d is not three tab-separated numbers", path, row + 2);
    }
    gzclose(f);
    if (row != nprobe) {
        sesame_titration_free(&t);
        return sesame__fail(err, SESAME_ERR_FORMAT,
            "%s has %d rows, the beta matrix has %d probes -- lineage mismatch",
            path, row, nprobe);
    }
    t.n = nprobe;
    *out = t;
    return SESAME_OK;
}

void sesame_titration_free(sesame_titration_t *t)
{
    if (!t) return;
    free(t->floor); free(t->amp); free(t->rho);
    t->floor = t->amp = t->rho = NULL; t->n = 0;
}

int32_t sesame_titration_apply(double *mat, int32_t nprobe, int32_t nsamp,
                               const sesame_titration_t *t, double gamma,
                               int32_t *ncal, int32_t *nzero)
{
    int32_t i, j, nval = 0;
    double *rg = (double *)malloc((size_t)nprobe * sizeof(double));

    for (i = 0; i < nprobe; i++) {
        rg[i] = isnan(t->rho[i]) ? NAN : pow(t->rho[i], gamma);
        if (!isnan(rg[i])) nval++;
    }
    for (j = 0; j < nsamp; j++) {
        double *b = mat + (size_t)j * (size_t)nprobe;
        int32_t nc = 0, nz = 0;
        for (i = 0; i < nprobe; i++) {
            double bp;
            if (isnan(b[i]) || isnan(rg[i])) { b[i] = NAN; continue; }
            nc++;
            bp = (b[i] - t->floor[i]) / t->amp[i];
            if (bp <= 0)      { b[i] = 0; nz++; }
            else if (bp >= 1) b[i] = 1;
            else              b[i] = bp / (bp + rg[i] * (1 - bp));
        }
        if (ncal)  ncal[j] = nc;
        if (nzero) nzero[j] = nz;
    }
    free(rg);
    return nval;
}
