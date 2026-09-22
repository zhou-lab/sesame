/* internal.h -- shared across libsesame TUs. Not installed, not public.
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present The Children's Hospital of Philadelphia
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#ifndef SESAME_INTERNAL_H
#define SESAME_INTERNAL_H

#include "sesame.h"

/* Fills err (if non-NULL) and returns code, so callers can `return
 * sesame__fail(...)`. Never longjmps, never writes to stderr. */
int sesame__fail(sesame_err_t *err, int code, const char *fmt, ...);

/* --- numerics: each replicates a specific R semantic exactly --- */

/* Copy non-NaN values of x into out (>= n); returns the count. Mirrors R's
 * na.rm=TRUE / sort(na.last=NA): NAs drop *before* n is taken. */
int32_t sesame__drop_na(const double *x, int32_t n, double *out);

/* R's quantile(x, p, type=7). _sorted takes ascending NaN-free data; the other
 * sorts x in place first. */
double sesame__quantile7_sorted(const double *x, int32_t n, double p);
double sesame__quantile7(double *x, int32_t n, double p);

/* R's median() -- identical to quantile type 7 at p=0.5. */
double sesame__median_sorted(const double *x, int32_t n);

void sesame__sort(double *x, int32_t n);

/* preprocessCore::normalize.quantiles.use.target, CLEAN-ROOM -- characterized by
 * black-box probing rather than vendored, to keep sesame-cli single-licence and
 * free of LGPL's relink obligation on static linking. See NUMERICS.md (D8).
 * x_sorted must be ascending and NaN-free; tgt is sorted in place; out holds n
 * doubles. */
void sesame__qnorm_use_target(const double *x_sorted, int32_t n,
                              double *tgt, int32_t m, double *out);

/* R's pmax(a, b) with NO na.rm -- NA propagates. Note C's fmax() does the
 * opposite and would silently disagree. */
double sesame__pmax2(double a, double b);

/* R's pmax(a,b,na.rm=TRUE): NA dropped; both NA -> NA. */
double sesame__pmax2_narm(double a, double b);

/* # elements <= x in an ascending NaN-free array (k in ecdf F(x)=k/n). */
int32_t sesame__count_le(const double *sorted, int32_t n, double x);

/* Inverse Mills ratio phi(t)/Phi(t) (D5): erfc for t>-5, Laplace continued
 * fraction below. Better conditioned than R's difference-of-logs. */
double sesame__inv_mills(double t);

/* R's normExpSignal via the inverse Mills ratio; floors negatives at 1e-6. Does
 * NOT add the noob offset. NaN x -> NaN. */
double sesame__norm_exp_signal(double mu, double sigma, double alpha, double x);

/* MASS::huber fixed-scale location. y sorted in place & must be NaN-free;
 * scratch >= n. D6: mad==0 falls back to max(IQR/1.349, 1e-8), *mad0 set. */
void sesame__huber(double *y, int32_t n, double *scratch,
                   double k, double tol, double *mu_out, double *s_out,
                   int *mad0);

/* Case-insensitive: does Probe_ID contain "negative" (R's neg-control grep)? */
int sesame__is_neg_control(const char *id);

/* Regularized incomplete beta I_x(a,b) and the t/F tail CDFs built on it (DML). */
double sesame__betai(double a, double b, double x);
double sesame__pt_2sided(double t, double df);       /* Pr(>|t|) */
double sesame__pf_upper(double f, double d1, double d2);

/* Householder QR least squares (dml.c). A is n*k row-major (destroyed), b is n
 * (destroyed). Fills beta[k], *rss, upper-triangular R[k*k]; v is scratch >= n.
 * Returns 0, or -1 if rank-deficient. Shared by dml and cnv. */
int sesame__ols(double *A, int32_t n, int32_t k, double *b,
                double *beta, double *rss, double *R, double *v);

/* Circular binary segmentation (cbs.c), a deterministic port of DNAcopy's core.
 * Segments x[0..n-1]; fills seg_end[] (cumulative 1-based segment ends, ascending)
 * and returns the segment count. min_width/alpha/kmax/nmin/ngrid/tol mirror
 * DNAcopy::segment; trimmed_sd>0 and undo_sd>0 enable the sdundo merge. */
int32_t sesame__cbs(const double *x, int32_t n, int32_t min_width, double alpha,
                    int32_t kmax, int32_t nmin, int ngrid, double tol,
                    double trimmed_sd, double undo_sd,
                    int32_t *seg_end, int32_t max_seg);
double sesame__trimmed_sd(const double *x, int32_t n, double trim);

/* Registry lookup (registry.h must be included by the TU that uses the type). */
struct sesame_reg_t;

/* Index accessors -- keep the on-disk layout private to index.c. */
const uint32_t *sesame__index_M(const sesame_index_t *ix);
const uint32_t *sesame__index_U(const sesame_index_t *ix);
const uint8_t  *sesame__index_col(const sesame_index_t *ix);
const uint8_t  *sesame__index_mask(const sesame_index_t *ix);

#endif /* SESAME_INTERNAL_H */
