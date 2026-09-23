/* dml.c -- per-probe differential methylation (sesame's DML).
 *
 * For each probe, an ordinary least-squares fit of its betas across samples on a
 * design matrix, giving per-coefficient estimates / t-tests and a holdout F-test
 * per categorical variable. The fit is a Householder QR least squares; p-values
 * come from the t and F tail CDFs (numerics.c). Unlike the preprocessing steps
 * this has no data-lineage caveat: the same betas matrix in gives the same stats
 * as R's lm, to ~1e-9.
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

#include <math.h>
#include <stdlib.h>
#include <string.h>

struct sesame_dml_work_t {
    int32_t m, p;
    double *Xs, *ys;          /* compacted (non-NA) design + response */
    double *qr, *qty, *beta;  /* QR working copy, Q'y, coefficients */
    double *R, *Rinv;         /* p*p upper-triangular factor and its inverse */
    double *v;                /* Householder reflector, length m */
    double *Xred, *qtyr, *betar, *Rr;  /* reduced-model scratch (F-test) */
};

sesame_dml_work_t *sesame_dml_work_new(int32_t m, int32_t p)
{
    sesame_dml_work_t *w = (sesame_dml_work_t *)calloc(1, sizeof *w);
    size_t ms = (size_t)m, ps = (size_t)p, mp = ms * ps, pp = ps * ps;
    if (!w) return NULL;
    w->m = m; w->p = p;
    w->Xs   = (double *)malloc(mp * sizeof(double));
    w->ys   = (double *)malloc(ms * sizeof(double));
    w->qr   = (double *)malloc(mp * sizeof(double));
    w->qty  = (double *)malloc(ms * sizeof(double));
    w->beta = (double *)malloc(ps * sizeof(double));
    w->R    = (double *)malloc(pp * sizeof(double));
    w->Rinv = (double *)malloc(pp * sizeof(double));
    w->v    = (double *)malloc(ms * sizeof(double));
    w->Xred = (double *)malloc(mp * sizeof(double));
    w->qtyr = (double *)malloc(ms * sizeof(double));
    w->betar= (double *)malloc(ps * sizeof(double));
    w->Rr   = (double *)malloc(pp * sizeof(double));
    if (!w->Xs || !w->ys || !w->qr || !w->qty || !w->beta || !w->R || !w->Rinv ||
        !w->v || !w->Xred || !w->qtyr || !w->betar || !w->Rr) {
        sesame_dml_work_free(w); return NULL;
    }
    return w;
}

void sesame_dml_work_free(sesame_dml_work_t *w)
{
    if (!w) return;
    free(w->Xs); free(w->ys); free(w->qr); free(w->qty); free(w->beta);
    free(w->R); free(w->Rinv); free(w->v); free(w->Xred); free(w->qtyr);
    free(w->betar); free(w->Rr);
    free(w);
}

/* Householder QR least squares. A is n*k row-major (destroyed), b is n
 * (destroyed, becomes Q'b). On success fills beta[k], *rss, and the k*k
 * upper-triangular R (row-major). v is scratch of length >= n. Returns 0, or -1
 * if rank-deficient (a zero pivot). Shared with cnv.c (target ~ normals fit). */
int sesame__ols(double *A, int32_t n, int32_t k, double *b,
                double *beta, double *rss, double *R, double *v)
{
    int32_t i, j, c;
    for (j = 0; j < k; j++) {
        double norm = 0.0, ajj, alpha, vn2, r;
        for (i = j; i < n; i++) { double a = A[i*k+j]; norm += a*a; }
        norm = sqrt(norm);
        if (norm == 0.0) return -1;
        ajj = A[j*k+j];
        alpha = (ajj >= 0.0) ? -norm : norm;
        for (i = j; i < n; i++) v[i] = A[i*k+j];
        v[j] -= alpha;
        vn2 = 0.0; for (i = j; i < n; i++) vn2 += v[i]*v[i];
        if (vn2 > 0.0) {
            for (c = j; c < k; c++) {
                double dot = 0.0, f;
                for (i = j; i < n; i++) dot += v[i]*A[i*k+c];
                f = 2.0*dot/vn2;
                for (i = j; i < n; i++) A[i*k+c] -= f*v[i];
            }
            { double dot = 0.0, f;
              for (i = j; i < n; i++) dot += v[i]*b[i];
              f = 2.0*dot/vn2;
              for (i = j; i < n; i++) b[i] -= f*v[i]; }
        }
        for (c = 0; c < j; c++) R[j*k+c] = 0.0;
        for (c = j; c < k; c++) R[j*k+c] = A[j*k+c];
        r = R[j*k+j];
        if (r == 0.0 || !isfinite(r)) return -1;
    }
    for (j = k-1; j >= 0; j--) {                 /* back-substitute R beta = b */
        double s = b[j];
        for (c = j+1; c < k; c++) s -= R[j*k+c]*beta[c];
        beta[j] = s / R[j*k+j];
    }
    { double s = 0.0; for (i = k; i < n; i++) s += b[i]*b[i]; *rss = s; }
    return 0;
}

/* Invert the k*k upper-triangular R into Rinv (row-major). */
static void invert_upper(const double *R, double *Rinv, int32_t k)
{
    int32_t i, j, l;
    for (i = 0; i < k; i++) for (j = 0; j < k; j++) Rinv[i*k+j] = 0.0;
    for (j = 0; j < k; j++) Rinv[j*k+j] = 1.0 / R[j*k+j];
    for (i = k-1; i >= 0; i--)
        for (j = i+1; j < k; j++) {
            double s = 0.0;
            for (l = i+1; l <= j; l++) s += R[i*k+l]*Rinv[l*k+j];
            Rinv[i*k+j] = -s / R[i*k+i];
        }
}

int32_t sesame_dml_fit(const sesame_dml_design_t *d, sesame_dml_work_t *w,
                       const double *y, double *est, double *pval,
                       double *fpval, double *eff)
{
    int32_t m = d->m, p = d->p, nobs = 0, i, j, v, df;
    double rss, sigma2;

    for (j = 0; j < p; j++) { est[j] = NAN; pval[j] = NAN; }
    for (v = 0; v < d->nvar; v++) { fpval[v] = NAN; eff[v] = NAN; }

    for (i = 0; i < m; i++) {                     /* drop NA samples */
        if (isnan(y[i])) continue;
        memcpy(w->Xs + (size_t)nobs*(size_t)p, d->X + (size_t)i*(size_t)p,
               (size_t)p*sizeof(double));
        w->ys[nobs] = y[i];
        nobs++;
    }
    df = nobs - p;
    if (df <= 0) return 0;

    memcpy(w->qr, w->Xs, (size_t)nobs*(size_t)p*sizeof(double));
    memcpy(w->qty, w->ys, (size_t)nobs*sizeof(double));
    if (sesame__ols(w->qr, nobs, p, w->qty, w->beta, &rss, w->R, w->v) != 0) return 0;

    sigma2 = rss / (double)df;
    invert_upper(w->R, w->Rinv, p);
    for (j = 0; j < p; j++) {
        double d2 = 0.0, se, t;
        int32_t l;
        for (l = j; l < p; l++) d2 += w->Rinv[j*p+l]*w->Rinv[j*p+l];  /* (X'X)^-1_jj */
        se = sqrt(sigma2 * d2);
        est[j] = w->beta[j];
        t = w->beta[j] / se;
        pval[j] = sesame__pt_2sided(t, (double)df);
    }

    for (v = 0; v < d->nvar; v++) {               /* holdout F-test per variable */
        int32_t k0 = d->var_off[v], k1 = d->var_off[v+1], kv = k1 - k0;
        int32_t pr = p - kv, cc = 0, l;
        double rssr, F, mx = 0.0, mn = 0.0;
        for (j = 0; j < p; j++) {                 /* copy columns not in var v */
            int in = 0;
            for (l = k0; l < k1; l++) if (d->var_col[l] == j) { in = 1; break; }
            if (in) continue;
            for (i = 0; i < nobs; i++)
                w->Xred[(size_t)i*(size_t)pr + (size_t)cc] = w->Xs[(size_t)i*(size_t)p + (size_t)j];
            cc++;
        }
        memcpy(w->qtyr, w->ys, (size_t)nobs*sizeof(double));
        if (sesame__ols(w->Xred, nobs, pr, w->qtyr, w->betar, &rssr, w->Rr, w->v) != 0)
            continue;
        F = ((rssr - rss) / (double)kv) / (rss / (double)df);
        if (F < 0.0) F = 0.0;                     /* guard tiny negative from roundoff */
        fpval[v] = sesame__pf_upper(F, (double)kv, (double)df);
        for (l = k0; l < k1; l++) {               /* effect size: max(est,0)-min(est,0) */
            double e = est[d->var_col[l]];
            if (e > mx) mx = e;
            if (e < mn) mn = e;
        }
        eff[v] = mx - mn;
    }
    return nobs;
}
