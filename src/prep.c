/* prep.c -- preprocessing steps (the prepSesame codes).
 *
 * Implemented: Q, C, D, E, P, B (the full QCDPB pipeline, plus E = dyeBiasL).
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present The Children's Hospital of Philadelphia
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "sesame.h"
#include "internal.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------- Q ---
 *
 * qualityMask (R/mask.R:286-297): OR a precomputed recommended-mask vector into
 * the SigDF mask. The vector is loaded once per platform via sesame_quality_mask
 * (which shells out to yame); here we just apply it. It is aligned to the
 * ordering, so mask[i] corresponds to probe i. */
int sesame_prep_quality_mask(sesame_sigdf_t *s, const uint8_t *qmask,
                             int32_t qn, sesame_err_t *err)
{
    int32_t i;
    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!s || !qmask) return sesame__fail(err, SESAME_ERR_IO, "null argument");
    if (qn != s->n)
        return sesame__fail(err, SESAME_ERR_FORMAT,
            "mask length %d != probe count %d (ordering/mask mismatch)", qn, s->n);
    for (i = 0; i < s->n; i++) if (qmask[i]) s->mask[i] = 1;
    return SESAME_OK;
}

/* ---------------------------------------------------------------- P ---
 *
 * pOOBAH (R/detection.R:141-170) -- detection-p masking from the empirical
 * distribution of out-of-band + negative-control signal.
 *
 * Background per channel:
 *   bgG = Grn signal (MG,UG) of col==R Inf-I probes NOT in the background mask
 *         + Grn (UG) of negative-control probes;
 *   bgR = Red signal (MR,UR) of col==G Inf-I probes not in the background mask
 *         + Red (UR) of negative controls.
 * If a channel has <= 100 non-NA background values, R substitutes the prior
 * 1:1000 (R/detection.R:155-156).
 *
 * Per probe: redmax = pmax(MR,UR, na.rm), grnmax = pmax(MG,UG, na.rm); the
 * detection p is the more significant of 1-F_R(redmax), 1-F_G(grnmax), where F
 * is the ecdf of the channel background. Mask if p > pval_threshold.
 *
 * D2: R computes pmin() with NO na.rm, so a channel whose alleles are both NA
 * makes p NA, and `NA > 0.05` does not mask -- the probe silently survives. The
 * sibling ELBAR sets such p to 1.0. We do the same: pmin ignores an NA channel,
 * and if BOTH are NA the probe is masked (p = 1).
 */
int sesame__is_neg_control(const char *id)
{
    /* case-insensitive substring "negative" (R greps tolower(Probe_ID)) */
    const char *n = "negative";
    size_t li = strlen(id), ln = strlen(n), i, j;
    if (li < ln) return 0;
    for (i = 0; i + ln <= li; i++) {
        for (j = 0; j < ln; j++) {
            char c = id[i + j];
            if (c >= 'A' && c <= 'Z') c = (char)(c - 'A' + 'a');
            if (c != n[j]) break;
        }
        if (j == ln) return 1;
    }
    return 0;
}

/* 1 - F(x): detection p from a sorted, NaN-free background. R evaluates
 * `1 - k/n` (not (n-k)/n) -- match that arithmetic bit-for-bit. */
static double detect_p(const double *sorted, int32_t n, double x)
{
    int32_t k;
    if (isnan(x)) return NAN;
    k = sesame__count_le(sorted, n, x);
    return 1.0 - (double)k / (double)n;
}

/* Detection p-value per probe (pOOBAH) into pout[n]. Same background pools as the
 * P step; D2: a probe with no signal in either channel gets p=1 (so it would
 * mask). Shared by the P step and the `pval` output; does not mutate the SigDF. */
int sesame_poobah_pvals(const sesame_sigdf_t *s, const uint8_t *bgmask, int32_t bn,
                        int combine_neg, double *pout, sesame_err_t *err)
{
    int32_t n, i, nG = 0, nR = 0;
    double *bgG = NULL, *bgR = NULL, *tmp = NULL;
    int rc = SESAME_OK;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!s || !bgmask || !pout) return sesame__fail(err, SESAME_ERR_IO, "null argument");
    n = s->n;
    if (bn != n)
        return sesame__fail(err, SESAME_ERR_FORMAT,
            "background mask length %d != probe count %d", bn, n);

    /* bg pools: oob (2 per Inf-I probe) + neg controls (1 per probe), so 3n is a
     * safe upper bound for each channel. */
    bgG = (double *)malloc((size_t)n * 3 * sizeof(double));
    bgR = (double *)malloc((size_t)n * 3 * sizeof(double));
    tmp = (double *)malloc((size_t)n * 3 * sizeof(double));
    if (!bgG || !bgR || !tmp) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto out; }

    for (i = 0; i < n; i++) {
        if (bgmask[i]) continue;                       /* excluded from bg */
        if (s->col[i] == SESAME_COL_R) {               /* oobG */
            bgG[nG++] = s->MG[i]; bgG[nG++] = s->UG[i];
        } else if (s->col[i] == SESAME_COL_G) {        /* oobR */
            bgR[nR++] = s->MR[i]; bgR[nR++] = s->UR[i];
        }
    }
    if (combine_neg) {
        for (i = 0; i < n; i++) {
            if (sesame__is_neg_control(sesame_index_probe_id(s->ix, i))) {
                bgG[nG++] = s->UG[i];
                bgR[nR++] = s->UR[i];
            }
        }
    }

    /* drop NA (na.last=NA), then the <=100 empirical-prior fallback, then sort */
    nG = sesame__drop_na(bgG, nG, tmp); memcpy(bgG, tmp, (size_t)nG * sizeof(double));
    nR = sesame__drop_na(bgR, nR, tmp); memcpy(bgR, tmp, (size_t)nR * sizeof(double));
    if (nG <= 100) { for (i = 0; i < 1000; i++) bgG[i] = i + 1; nG = 1000; }
    if (nR <= 100) { for (i = 0; i < 1000; i++) bgR[i] = i + 1; nR = 1000; }
    sesame__sort(bgG, nG);
    sesame__sort(bgR, nR);

    for (i = 0; i < n; i++) {
        double pr = detect_p(bgR, nR, sesame__pmax2_narm(s->MR[i], s->UR[i]));
        double pg = detect_p(bgG, nG, sesame__pmax2_narm(s->MG[i], s->UG[i]));
        if (isnan(pr) && isnan(pg)) pout[i] = 1.0;     /* D2: both NA -> mask */
        else if (isnan(pr)) pout[i] = pg;              /* D2: use the other  */
        else if (isnan(pg)) pout[i] = pr;
        else pout[i] = pr < pg ? pr : pg;              /* pmin */
    }

out:
    free(bgG); free(bgR); free(tmp);
    return rc;
}

int sesame_prep_poobah(sesame_sigdf_t *s, const uint8_t *bgmask, int32_t bn,
                       double pval_threshold, int combine_neg, sesame_err_t *err)
{
    double *pv;
    int32_t i, n;
    int rc;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!s) return sesame__fail(err, SESAME_ERR_IO, "null sigdf");
    n = s->n;
    if (!(pv = (double *)malloc((size_t)n * sizeof(double))))
        return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
    rc = sesame_poobah_pvals(s, bgmask, bn, combine_neg, pv, err);
    if (rc == SESAME_OK)
        for (i = 0; i < n; i++) if (pv[i] > pval_threshold) s->mask[i] = 1;
    free(pv);
    return rc;
}

/* detectionPnegEcdf (R/detection.R:100-118): detection p from the empirical CDF
 * of the negative-control probes -- funcG over their Grn (UG), funcR over their
 * Red (UR). Per probe p = pmin(1-funcR(pmax(MR,UR)), 1-funcG(pmax(MG,UG))) with
 * na.rm (an NA channel keeps the other's p), and no signal in either channel
 * gives p = 1. Into pout[n]; does not mutate the SigDF. An alternative to pOOBAH
 * that needs only the negative controls (no out-of-band pool). */
int sesame_detection_pneg_pvals(const sesame_sigdf_t *s, double *pout,
                                sesame_err_t *err)
{
    int32_t n, i, nG = 0, nR = 0;
    double *bgG = NULL, *bgR = NULL, *tmp = NULL;
    int rc = SESAME_OK;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!s || !pout) return sesame__fail(err, SESAME_ERR_IO, "null argument");
    n = s->n;

    /* neg-control pools: G = UG, R = UR of every "negative" probe. */
    bgG = (double *)malloc((size_t)n * sizeof(double));
    bgR = (double *)malloc((size_t)n * sizeof(double));
    tmp = (double *)malloc((size_t)n * sizeof(double));
    if (!bgG || !bgR || !tmp) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto out; }

    for (i = 0; i < n; i++)
        if (sesame__is_neg_control(sesame_index_probe_id(s->ix, i))) {
            bgG[nG++] = s->UG[i];
            bgR[nR++] = s->UR[i];
        }
    if (nG == 0 || nR == 0) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT,
                          "no negative-control probes for pneg detection");
        goto out;
    }
    /* R's ecdf sorts and drops NA (na.last=NA). */
    nG = sesame__drop_na(bgG, nG, tmp); memcpy(bgG, tmp, (size_t)nG * sizeof(double));
    nR = sesame__drop_na(bgR, nR, tmp); memcpy(bgR, tmp, (size_t)nR * sizeof(double));
    sesame__sort(bgG, nG);
    sesame__sort(bgR, nR);

    for (i = 0; i < n; i++) {
        double pr = detect_p(bgR, nR, sesame__pmax2_narm(s->MR[i], s->UR[i]));
        double pg = detect_p(bgG, nG, sesame__pmax2_narm(s->MG[i], s->UG[i]));
        if (isnan(pr) && isnan(pg)) pout[i] = 1.0;     /* both NA -> mask */
        else if (isnan(pr)) pout[i] = pg;
        else if (isnan(pg)) pout[i] = pr;
        else pout[i] = pr < pg ? pr : pg;              /* pmin */
    }

out:
    free(bgG); free(bgR); free(tmp);
    return rc;
}

/* Per-sample pipeline: apply prep to a copy of raw, filling whichever of
 * beta/M/U/total/pval are non-NULL (each nprobes long). Signal (M/U/total) is
 * from the preprocessed copy unless raw_signal, then from raw. pval is the
 * detection p -- pOOBAH by default, or (pneg_detection) detectionPnegEcdf --
 * taken at the P-step lineage point, or on the final signal if prep has no P.
 * Masking always uses pOOBAH. Shared by `preprocess` and the f64 harness. */
int sesame_pipeline(const sesame_sigdf_t *raw, const char *prep,
                    const uint8_t *qmask, int32_t qn,
                    const uint8_t *bgmask, int32_t bgn, int raw_signal,
                    int pneg_detection,
                    double *beta, double *M, double *U, double *total,
                    double *pval, uint8_t *col, sesame_err_t *err)
{
    int32_t n, k;
    sesame_sigdf_t *sig = NULL;
    double *pv = NULL;
    int rc = SESAME_OK, pval_done = 0;
    const char *p;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!raw || !prep) return sesame__fail(err, SESAME_ERR_IO, "null argument");
    n = raw->n;

    if (raw_signal && (M || U || total)) {
        if (M || U) rc = sesame_signal_mu(raw, M, U, err);
        if (rc == SESAME_OK && total) {
            if (M && U) for (k = 0; k < n; k++) total[k] = M[k] + U[k];
            else rc = sesame_total_intensities(raw, total, err);
        }
        if (rc != SESAME_OK) return rc;
    }

    if (!(sig = sesame_sigdf_dup(raw)))
        return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
    if (strchr(prep, 'P')) {   /* pv holds pOOBAH values for the P-step mask */
        if (!(pv = (double *)malloc((size_t)n * sizeof(double)))) {
            sesame_sigdf_free(sig);
            return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
        }
    }

    for (p = prep; *p && rc == SESAME_OK; p++) {
        switch (*p) {
        case 'Q': rc = sesame_prep_quality_mask(sig, qmask, qn, err); break;
        case 'C': rc = sesame_prep_infer_channel(sig, 0, 0, err); break;
        case 'D': rc = sesame_prep_dye_bias_nl(sig, err); break;
        case 'E': rc = sesame_prep_dye_bias_l(sig, err); break;
        case 'P':
            /* Masking is always pOOBAH (that is what the P code means). The pval
             * OUTPUT is taken at this lineage point too -- pneg if requested. */
            rc = sesame_poobah_pvals(sig, bgmask, bgn, 1, pv, err);
            if (rc == SESAME_OK) {
                for (k = 0; k < n; k++) if (pv[k] > 0.05) sig->mask[k] = 1;
                if (pval) {
                    if (pneg_detection)
                        rc = sesame_detection_pneg_pvals(sig, pval, err);
                    else
                        memcpy(pval, pv, (size_t)n * sizeof(double));
                    pval_done = 1;
                }
            }
            break;
        case 'B': rc = sesame_prep_noob(sig, bgmask, bgn, 1, 15.0, err); break;
        default:  rc = sesame__fail(err, SESAME_ERR_UNSUPPORTED, "prep code '%c'", *p);
        }
    }
    if (rc == SESAME_OK && pval && !pval_done)   /* no P: detect on final signal */
        rc = pneg_detection ? sesame_detection_pneg_pvals(sig, pval, err)
                            : sesame_poobah_pvals(sig, bgmask, bgn, 1, pval, err);
    if (rc == SESAME_OK && col)  memcpy(col, sig->col, (size_t)n);
    if (rc == SESAME_OK && beta) rc = sesame_get_betas(sig, 1, beta, err);
    if (rc == SESAME_OK && !raw_signal && (M || U || total)) {
        if (M || U) rc = sesame_signal_mu(sig, M, U, err);
        if (rc == SESAME_OK && total) {
            if (M && U) for (k = 0; k < n; k++) total[k] = M[k] + U[k];
            else rc = sesame_total_intensities(sig, total, err);
        }
    }

    free(pv);
    sesame_sigdf_free(sig);
    return rc;
}

/* ---------------------------------------------------------------- C ---
 *
 * inferInfiniumIChannel (R/channel_inference.R:20-55).
 *
 * For every Infinium-I probe, take the larger of the two alleles in each
 * channel and assign the probe to whichever channel is brighter. The
 * "background" is then the pooled *out-of-band* signal under that new
 * assignment, and any probe whose brighter channel does not clear the 95th
 * percentile of it is deemed failed.
 *
 * Four R semantics that a naive port gets wrong, all verified against the
 * source:
 *
 *   1. `pmax(MR, UR)` has NO na.rm (:26-27) -- an NA allele makes the max NA.
 *      C's fmax() would instead ignore the NaN and return the other operand.
 *      sesame__pmax2 replicates R.
 *   2. `ifelse(red_max > grn_max, "R", "G")` (:28) yields **NA**, not "G", when
 *      the comparison is NA. Those probes are neither R nor G at this stage.
 *   3. `sdf1[new_col == "R",]` (:30-31) indexes with a logical containing NAs,
 *      which yields NA-*filled rows*; `quantile(..., na.rm=TRUE)` then drops
 *      them. Net effect: NA probes contribute nothing to the background pool.
 *      Not an accident to be tidied -- it is the behaviour.
 *   4. `red_max > grn_max` is strict, so exact ties go to **G**.
 *
 * With the defaults (switch_failed=0, mask_failed=0) failed probes revert to
 * the manifest channel and are not masked.
 */
int sesame_prep_infer_channel(sesame_sigdf_t *s, int switch_failed,
                              int mask_failed, sesame_err_t *err)
{
    int32_t n, i, n_inf1 = 0, n_pool = 0, n_ok;
    double *rmax = NULL, *gmax = NULL, *pool = NULL, *buf = NULL;
    uint8_t *newc = NULL;
    double bg_max;
    int rc = SESAME_OK;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!s) return sesame__fail(err, SESAME_ERR_IO, "null sigdf");
    n = s->n;

    rmax = (double  *)malloc((size_t)n * sizeof(double));
    gmax = (double  *)malloc((size_t)n * sizeof(double));
    newc = (uint8_t *)malloc((size_t)n);
    /* pool holds at most 2 values per Inf-I probe */
    pool = (double  *)malloc((size_t)n * 2 * sizeof(double));
    buf  = (double  *)malloc((size_t)n * 2 * sizeof(double));
    if (!rmax || !gmax || !newc || !pool || !buf) {
        rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom in C step");
        goto out;
    }

    /* Pass 1: per-probe maxima and the tentative assignment.
     * newc is SESAME_COL_II for "NA" here -- Inf-I probes can never legitimately
     * be II, so it is an unambiguous sentinel for R's NA. */
    for (i = 0; i < n; i++) {
        if (s->col[i] == SESAME_COL_II) { newc[i] = SESAME_COL_II; continue; }
        n_inf1++;
        rmax[i] = sesame__pmax2(s->MR[i], s->UR[i]);
        gmax[i] = sesame__pmax2(s->MG[i], s->UG[i]);
        if (isnan(rmax[i]) || isnan(gmax[i]))
            newc[i] = SESAME_COL_II;               /* R's NA (semantic 2) */
        else
            newc[i] = (rmax[i] > gmax[i]) ? SESAME_COL_R  /* strict: ties -> G */
                                          : SESAME_COL_G;
    }
    if (n_inf1 == 0) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT, "no Infinium-I probes");
        goto out;
    }

    /* Pass 2: background = out-of-band signal under the NEW assignment.
     * Probes whose newc is the NA sentinel contribute nothing (semantic 3). */
    for (i = 0; i < n; i++) {
        if (newc[i] == SESAME_COL_R) {          /* inferred red -> green is oob */
            pool[n_pool++] = s->MG[i];
            pool[n_pool++] = s->UG[i];
        } else if (newc[i] == SESAME_COL_G && s->col[i] != SESAME_COL_II) {
            pool[n_pool++] = s->MR[i];          /* inferred green -> red is oob */
            pool[n_pool++] = s->UR[i];
        }
    }

    n_ok = sesame__drop_na(pool, n_pool, buf);
    if (n_ok == 0) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT,
                          "no out-of-band signal for channel inference");
        goto out;
    }
    bg_max = sesame__quantile7(buf, n_ok, 0.95);

    /* Pass 3: failed probes revert to the manifest channel (default). */
    for (i = 0; i < n; i++) {
        int failed;
        if (s->col[i] == SESAME_COL_II) continue;

        failed = isnan(rmax[i]) || isnan(gmax[i]) ||
                 (sesame__pmax2(rmax[i], gmax[i]) < bg_max);

        if (newc[i] == SESAME_COL_II) {
            /* R's NA. With switch_failed=TRUE R would leave col=NA -- a
             * degenerate state our encoding cannot represent (and which would
             * silently turn an Inf-I probe into Inf-II). Revert instead. Moot
             * at the default switch_failed=FALSE, where `failed` is already
             * TRUE for these and reverts them anyway. */
            newc[i] = s->col[i];
        } else if (failed && !switch_failed) {
            newc[i] = s->col[i];
        }
        if (failed && mask_failed) s->mask[i] = 1;
        s->col[i] = newc[i];
    }

out:
    free(rmax); free(gmax); free(newc); free(pool); free(buf);
    return rc;
}

/* ---------------------------------------------------------------- D ---
 *
 * dyeBiasNL (R/dye_bias.R:118-167).
 *
 * Pull the Red and Green Infinium-I distributions to their common midpoint:
 * quantile-normalize Red onto Green, average that with Red-as-is, and use the
 * resulting IR1 -> IRmid map as a lookup applied to every Red value. Green is
 * handled symmetrically. Both maps are built from the ORIGINAL signals before
 * either is applied (R builds both closures first, R/dye_bias.R:164-165).
 */

/* A monotone lookup table: knots (x, y), linearly interpolated. */
typedef struct { double *x, *y; int32_t n; } knots_t;

/* Build knots from a sorted src and its aligned midpoint values, collapsing
 * tied src values. This is R's approx(ties=mean) (R/dye_bias.R:143): approx
 * needs unique x, so tied x collapse to one knot carrying the mean of their y.
 * Here the mean is a formality -- qnorm_use_target gives every element of a
 * tied run the same value (average rank), so y is constant across the run.
 * Verified on real data: 0 of 11259 tied runs on HM450 had non-constant y. */
static void knots_build(knots_t *k, const double *src, const double *mid,
                        int32_t n)
{
    int32_t a = 0;
    k->n = 0;
    while (a < n) {
        int32_t b = a, i;
        double sum = 0.0;
        while (b + 1 < n && src[b + 1] == src[a]) b++;
        for (i = a; i <= b; i++) sum += mid[i];
        k->x[k->n] = src[a];
        k->y[k->n] = sum / (double)(b - a + 1);
        k->n++;
        a = b + 1;
    }
}

/* Linear interpolation on the knots. Callers guarantee k->x[0] <= v <=
 * k->x[n-1], which is R's approx(rule=1) guarded by the explicit in-support
 * test, so no extrapolation is possible here. */
static double knots_interp(const knots_t *k, double v)
{
    int32_t lo = 0, hi = k->n - 1;

    if (v <= k->x[0])          return k->y[0];
    if (v >= k->x[k->n - 1])   return k->y[k->n - 1];

    while (hi - lo > 1) {                      /* largest lo with x[lo] <= v */
        int32_t mid = lo + ((hi - lo) >> 1);
        if (k->x[mid] <= v) lo = mid; else hi = mid;
    }
    if (k->x[lo] == v) return k->y[lo];
    return k->y[lo] + (v - k->x[lo]) / (k->x[lo + 1] - k->x[lo]) *
                      (k->y[lo + 1] - k->y[lo]);
}

/* R's fitfun (R/dye_bias.R:139-147): interpolate inside the support, shift
 * additively above it, scale multiplicatively below it. NA passes through. */
static double fitfun(const knots_t *k, double v,
                     double lo_in, double hi_in, double lo_mid, double hi_mid)
{
    if (isnan(v))   return v;
    if (v > hi_in)  return v - hi_in + hi_mid;          /* oversupp  */
    if (v < lo_in)  return lo_mid / lo_in * v;          /* undersupp */
    return knots_interp(k, v);                          /* insupp    */
}

/* Total intensity per probe, as totalIntensities(sdf, mask=FALSE) via signalMU:
 * col G -> MG+UG, col R -> MR+UR, col II -> UG+UR (R/sesame.R:111-115,
 * R/SigDFMethods.R:69-81). */
static double total_intensity(const sesame_sigdf_t *s, int32_t i)
{
    switch (s->col[i]) {
    case SESAME_COL_G: return s->MG[i] + s->UG[i];
    case SESAME_COL_R: return s->MR[i] + s->UR[i];
    default:           return s->UG[i] + s->UR[i];
    }
}

/* RGdistort = (topR/topG) / (medR/medG), where med* is the median total
 * intensity of that channel's Inf-I probes and top* the median of its largest
 * 20 (R/QC.R:359-366). NA or > 10 means the green channel has failed, and
 * dyeBiasNL gives up and masks all Inf-I green probes. */
static double rg_distort(const sesame_sigdf_t *s, double *buf)
{
    double medR, medG, topR, topG;
    int32_t i, nR = 0, nG = 0;

    for (i = 0; i < s->n; i++)
        if (s->col[i] == SESAME_COL_R) {
            double v = total_intensity(s, i);
            if (!isnan(v)) buf[nR++] = v;
        }
    sesame__sort(buf, nR);
    if (nR == 0) return NAN;
    medR = sesame__median_sorted(buf, nR);
    topR = sesame__median_sorted(buf + (nR > 20 ? nR - 20 : 0),
                                 nR > 20 ? 20 : nR);

    for (i = 0; i < s->n; i++)
        if (s->col[i] == SESAME_COL_G) {
            double v = total_intensity(s, i);
            if (!isnan(v)) buf[nG++] = v;
        }
    sesame__sort(buf, nG);
    if (nG == 0) return NAN;
    medG = sesame__median_sorted(buf, nG);
    topG = sesame__median_sorted(buf + (nG > 20 ? nG - 20 : 0),
                                 nG > 20 ? 20 : nG);

    return (topR / topG) / (medR / medG);
}

int sesame_prep_dye_bias_nl(sesame_sigdf_t *s, sesame_err_t *err)
{
    int32_t n, i, nIG = 0, nIR = 0, nIG1, nIR1;
    double *IG0 = NULL, *IR0 = NULL, *IG1 = NULL, *IR1 = NULL;
    double *IG2 = NULL, *IR2 = NULL, *IGmid = NULL, *IRmid = NULL, *tbuf = NULL;
    knots_t kR = {0}, kG = {0};
    double maxIG, minIG, maxIR, minIR, maxIGmid, minIGmid, maxIRmid, minIRmid;
    double dist;
    int rc = SESAME_OK;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!s) return sesame__fail(err, SESAME_ERR_IO, "null sigdf");
    n = s->n;

    tbuf = (double *)malloc((size_t)n * sizeof(double));
    if (!tbuf) return sesame__fail(err, SESAME_ERR_NOMEM, "oom in D step");

    /* Guard: green channel failed completely -> mask Inf-I green and stop
     * (R/dye_bias.R:122-123, maskIG at :94-97). R does this silently; we set a
     * status bit so a pipeline can count it. */
    dist = rg_distort(s, tbuf);
    if (isnan(dist) || dist > 10.0) {
        for (i = 0; i < n; i++)
            if (s->col[i] == SESAME_COL_G) s->mask[i] = 1;
        s->status |= SESAME_STAT_DYEBIAS_FAILED;
        free(tbuf);
        return SESAME_OK;
    }

    /* Pools. mask=TRUE is the default and means "include masked probes"
     * (R/dye_bias.R:126) -- so this is every Inf-I probe, masked or not. */
    IG0 = (double *)malloc((size_t)n * 2 * sizeof(double));
    IR0 = (double *)malloc((size_t)n * 2 * sizeof(double));
    if (!IG0 || !IR0) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto out; }
    for (i = 0; i < n; i++) {
        if (s->col[i] == SESAME_COL_G) {
            IG0[nIG++] = s->MG[i]; IG0[nIG++] = s->UG[i];
        } else if (s->col[i] == SESAME_COL_R) {
            IR0[nIR++] = s->MR[i]; IR0[nIR++] = s->UR[i];
        }
    }
    if (nIG == 0 || nIR == 0) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT,
                          "no Infinium-I %s probes for dye bias",
                          nIG == 0 ? "green" : "red");
        goto out;
    }

    /* min/max over the raw pools, na.rm=TRUE (R/dye_bias.R:130-131). */
    maxIG = -INFINITY; minIG = INFINITY; maxIR = -INFINITY; minIR = INFINITY;
    for (i = 0; i < nIG; i++) if (!isnan(IG0[i])) {
        if (IG0[i] > maxIG) maxIG = IG0[i];
        if (IG0[i] < minIG) minIG = IG0[i];
    }
    for (i = 0; i < nIR; i++) if (!isnan(IR0[i])) {
        if (IR0[i] > maxIR) maxIR = IR0[i];
        if (IR0[i] < minIR) minIR = IR0[i];
    }
    if (maxIG <= 0.0 || maxIR <= 0.0) {   /* R/dye_bias.R:133 */
        rc = SESAME_OK;
        goto out;
    }

    IG1 = (double *)malloc((size_t)nIG * sizeof(double));
    IR1 = (double *)malloc((size_t)nIR * sizeof(double));
    IG2 = (double *)malloc((size_t)nIG * sizeof(double));
    IR2 = (double *)malloc((size_t)nIR * sizeof(double));
    IGmid = (double *)malloc((size_t)nIG * sizeof(double));
    IRmid = (double *)malloc((size_t)nIR * sizeof(double));
    kR.x = (double *)malloc((size_t)nIR * sizeof(double));
    kR.y = (double *)malloc((size_t)nIR * sizeof(double));
    kG.x = (double *)malloc((size_t)nIG * sizeof(double));
    kG.y = (double *)malloc((size_t)nIG * sizeof(double));
    if (!IG1 || !IR1 || !IG2 || !IR2 || !IGmid || !IRmid ||
        !kR.x || !kR.y || !kG.x || !kG.y) {
        rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto out;
    }

    /* sort() drops NA (na.last=NA), so IR1 can be shorter than IR0. */
    nIR1 = sesame__drop_na(IR0, nIR, IR1); sesame__sort(IR1, nIR1);
    nIG1 = sesame__drop_na(IG0, nIG, IG1); sesame__sort(IG1, nIG1);
    if (nIR1 == 0 || nIG1 == 0) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT, "no non-NA Inf-I signal");
        goto out;
    }

    /* Red onto Green, and Green onto Red. R re-sorts the result; the output is
     * already non-decreasing, so that sort is a no-op. tgt is sorted in place,
     * so pass a scratch copy of the raw pool. */
    memcpy(tbuf, IG0, (size_t)nIG * sizeof(double));
    {
        int32_t mt = sesame__drop_na(tbuf, nIG, tbuf);   /* NA target dropped */
        sesame__qnorm_use_target(IR1, nIR1, tbuf, mt, IR2);
    }
    for (i = 0; i < nIR1; i++) IRmid[i] = (IR1[i] + IR2[i]) / 2.0;

    memcpy(tbuf, IR0, (size_t)nIR * sizeof(double));
    {
        int32_t mt = sesame__drop_na(tbuf, nIR, tbuf);
        sesame__qnorm_use_target(IG1, nIG1, tbuf, mt, IG2);
    }
    for (i = 0; i < nIG1; i++) IGmid[i] = (IG1[i] + IG2[i]) / 2.0;

    minIRmid = INFINITY; maxIRmid = -INFINITY;
    for (i = 0; i < nIR1; i++) {
        if (IRmid[i] < minIRmid) minIRmid = IRmid[i];
        if (IRmid[i] > maxIRmid) maxIRmid = IRmid[i];
    }
    minIGmid = INFINITY; maxIGmid = -INFINITY;
    for (i = 0; i < nIG1; i++) {
        if (IGmid[i] < minIGmid) minIGmid = IGmid[i];
        if (IGmid[i] > maxIGmid) maxIGmid = IGmid[i];
    }

    knots_build(&kR, IR1, IRmid, nIR1);
    knots_build(&kG, IG1, IGmid, nIG1);

    /* Both maps are built from the original signals before either is applied. */
    for (i = 0; i < n; i++) {
        s->MR[i] = fitfun(&kR, s->MR[i], minIR, maxIR, minIRmid, maxIRmid);
        s->UR[i] = fitfun(&kR, s->UR[i], minIR, maxIR, minIRmid, maxIRmid);
        s->MG[i] = fitfun(&kG, s->MG[i], minIG, maxIG, minIGmid, maxIGmid);
        s->UG[i] = fitfun(&kG, s->UG[i], minIG, maxIG, minIGmid, maxIGmid);
    }

out:
    free(IG0); free(IR0); free(IG1); free(IR1); free(IG2); free(IR2);
    free(IGmid); free(IRmid); free(tbuf);
    free(kR.x); free(kR.y); free(kG.x); free(kG.y);
    return rc;
}

/* ---------------------------------------------------------------- E ---
 *
 * dyeBiasL (R/dye_bias.R:191-208) -- linear dye-bias correction.
 *
 * ref = meanIntensity(sdf)                 (mean of unmasked in-band M,U)
 * fR  = ref / median(unmasked col==R {MR,UR})
 * fG  = ref / median(unmasked col==G {MG,UG})
 * then MG,UG *= fG and MR,UR *= fR over every probe (R scales whole columns).
 * signalMU/meanIntensity use mask=TRUE, i.e. masked probes are excluded from the
 * reference and the per-channel medians (but scaling still applies to all rows).
 */
int sesame_prep_dye_bias_l(sesame_sigdf_t *s, sesame_err_t *err)
{
    int32_t i, n, nR = 0, nG = 0;
    double *buf = NULL, ref, sum = 0.0, medR, medG, fR, fG;
    long cnt = 0;
    int rc = SESAME_OK;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!s) return sesame__fail(err, SESAME_ERR_IO, "null sigdf");
    n = s->n;

    /* ref = mean of unmasked in-band M,U over all probes (meanIntensity). */
    for (i = 0; i < n; i++) {
        double m, u;
        if (s->mask[i]) continue;
        if (s->col[i] == SESAME_COL_G)      { m = s->MG[i]; u = s->UG[i]; }
        else if (s->col[i] == SESAME_COL_R) { m = s->MR[i]; u = s->UR[i]; }
        else                                { m = s->UG[i]; u = s->UR[i]; }
        if (!isnan(m)) { sum += m; cnt++; }
        if (!isnan(u)) { sum += u; cnt++; }
    }
    if (cnt == 0) return SESAME_OK;   /* nothing to normalize; R -> NaN ref */
    ref = sum / (double)cnt;

    /* per-channel medians of unmasked in-band Inf-I {M,U}. */
    buf = (double *)malloc((size_t)n * 2 * sizeof(double));
    if (!buf) return sesame__fail(err, SESAME_ERR_NOMEM, "oom in E step");

    for (i = 0; i < n; i++)
        if (!s->mask[i] && s->col[i] == SESAME_COL_R) {
            if (!isnan(s->MR[i])) buf[nR++] = s->MR[i];
            if (!isnan(s->UR[i])) buf[nR++] = s->UR[i];
        }
    sesame__sort(buf, nR);
    medR = nR ? sesame__median_sorted(buf, nR) : NAN;

    for (i = 0; i < n; i++)
        if (!s->mask[i] && s->col[i] == SESAME_COL_G) {
            if (!isnan(s->MG[i])) buf[nG++] = s->MG[i];
            if (!isnan(s->UG[i])) buf[nG++] = s->UG[i];
        }
    sesame__sort(buf, nG);
    medG = nG ? sesame__median_sorted(buf, nG) : NAN;
    free(buf);

    fR = ref / medR;
    fG = ref / medG;

    /* R multiplies the whole MG/UG/MR/UR columns (NaN passes through). If a
     * channel has no Inf-I signal its factor is NaN; guard so the other channel
     * (and the Inf-II readout it drives) is not clobbered -- R would propagate
     * NaN, but a missing channel is a degenerate array, not a normal result. */
    if (isfinite(fG)) for (i = 0; i < n; i++) { s->MG[i] *= fG; s->UG[i] *= fG; }
    if (isfinite(fR)) for (i = 0; i < n; i++) { s->MR[i] *= fR; s->UR[i] *= fR; }

    return rc;
}

/* ---------------------------------------------------------------- B ---
 *
 * noob (R/background.R:85-123) -- normal-exponential background subtraction.
 *
 * Background (Normal) and foreground (Exponential) are parameterized per channel
 * from probes NOT in the background mask, then every signal is deconvolved via
 * the inverse Mills ratio (D5) and shifted by offset.
 *
 * Background pool (same construction as pOOBAH):
 *   bgG = oob Grn (MG,UG of col==R) + neg-control UG;
 *   bgR = oob Red (MR,UR of col==G) + neg-control UR.
 * Foreground pool (in-band signal):
 *   ibG = MG,UG of col==G + UG of col==II;
 *   ibR = MR,UR of col==R + UR of col==II.
 * R's SigDF carries control rows as col=="2" (chipAddressToSignal), so controls
 * legitimately enter the col==II foreground in both R and here.
 *
 * Fit: mu,sigma = huber(bg_capped); alpha = max(huber(ib).mu - mu, 10). Applied
 * to the FULL columns (all probes, masked or not), NaN passing through.
 */

/* Keep the non-NA background values strictly below median + 10*IQR (both na.rm),
 * R's cap against multi-mapping probes (R/background.R:103-104). R leaves the NA
 * elements in place but huber drops them, so we just drop them here. Survivors
 * are written back to v (sorted); returns the count. */
static int32_t bg_cap(double *v, int32_t n, double *scr)
{
    int32_t m = sesame__drop_na(v, n, scr), i, k = 0;
    double med, iqr, thr;
    if (m == 0) return 0;
    sesame__sort(scr, m);
    med = sesame__median_sorted(scr, m);
    iqr = sesame__quantile7_sorted(scr, m, 0.75) -
          sesame__quantile7_sorted(scr, m, 0.25);
    thr = med + 10.0 * iqr;
    for (i = 0; i < m; i++) if (scr[i] < thr) v[k++] = scr[i];
    return k;
}

int sesame_prep_noob(sesame_sigdf_t *s, const uint8_t *bgmask, int32_t bn,
                     int combine_neg, double offset, sesame_err_t *err)
{
    int32_t n, i, nbG = 0, nbR = 0, niG = 0, niR = 0, pG = 0, pR = 0, m;
    double *bgG = NULL, *bgR = NULL, *ibG = NULL, *ibR = NULL;
    double *scr = NULL, *scr2 = NULL;
    double muG, sgG, alG, muR, sgR, alR, mIbG, mIbR, dummy;
    int mad0, rc = SESAME_OK;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!s || !bgmask) return sesame__fail(err, SESAME_ERR_IO, "null argument");
    n = s->n;
    if (bn != n)
        return sesame__fail(err, SESAME_ERR_FORMAT,
            "background mask length %d != probe count %d", bn, n);

    bgG  = (double *)malloc((size_t)n * 3 * sizeof(double));
    bgR  = (double *)malloc((size_t)n * 3 * sizeof(double));
    ibG  = (double *)malloc((size_t)n * 2 * sizeof(double));
    ibR  = (double *)malloc((size_t)n * 2 * sizeof(double));
    scr  = (double *)malloc((size_t)n * 3 * sizeof(double));
    scr2 = (double *)malloc((size_t)n * 2 * sizeof(double));
    if (!bgG || !bgR || !ibG || !ibR || !scr || !scr2) {
        rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto out;
    }

    /* background: out-of-band of non-bgmask Inf-I probes */
    for (i = 0; i < n; i++) {
        if (bgmask[i]) continue;
        if (s->col[i] == SESAME_COL_R) {
            bgG[nbG++] = s->MG[i]; bgG[nbG++] = s->UG[i];
        } else if (s->col[i] == SESAME_COL_G) {
            bgR[nbR++] = s->MR[i]; bgR[nbR++] = s->UR[i];
        }
    }
    if (combine_neg) {
        for (i = 0; i < n; i++)
            if (sesame__is_neg_control(sesame_index_probe_id(s->ix, i))) {
                bgG[nbG++] = s->UG[i]; bgR[nbR++] = s->UR[i];
            }
    }

    /* foreground: in-band signal of non-bgmask probes (II incl. controls) */
    for (i = 0; i < n; i++) {
        if (bgmask[i]) continue;
        if (s->col[i] == SESAME_COL_G) {
            ibG[niG++] = s->MG[i]; ibG[niG++] = s->UG[i];
        } else if (s->col[i] == SESAME_COL_R) {
            ibR[niR++] = s->MR[i]; ibR[niR++] = s->UR[i];
        } else {                                   /* Inf-II */
            ibG[niG++] = s->UG[i]; ibR[niR++] = s->UR[i];
        }
    }

    /* R/background.R:98 -- need enough positive background (strict >0, NA out) */
    for (i = 0; i < nbG; i++) if (!isnan(bgG[i]) && bgG[i] > 0.0) pG++;
    for (i = 0; i < nbR; i++) if (!isnan(bgR[i]) && bgR[i] > 0.0) pR++;
    if (pG < 100 || pR < 100) {                    /* D7: flag, don't fail */
        s->status |= SESAME_STAT_NOOB_SKIPPED;
        goto out;
    }

    /* zeros -> 1 (R/background.R:101), then cap background at median+10*IQR */
    for (i = 0; i < nbG; i++) if (bgG[i] == 0.0) bgG[i] = 1.0;
    for (i = 0; i < nbR; i++) if (bgR[i] == 0.0) bgR[i] = 1.0;
    nbG = bg_cap(bgG, nbG, scr);
    nbR = bg_cap(bgR, nbR, scr);

    /* foreground zeros -> 1 (R/background.R:109-110) */
    for (i = 0; i < niG; i++) if (ibG[i] == 0.0) ibG[i] = 1.0;
    for (i = 0; i < niR; i++) if (ibR[i] == 0.0) ibR[i] = 1.0;

    /* fits: mu,sigma from background; alpha from foreground location */
    sesame__huber(bgG, nbG, scr, 1.5, 1e-6, &muG, &sgG, &mad0);
    if (mad0) s->status |= SESAME_STAT_NOOB_MAD0;
    m = sesame__drop_na(ibG, niG, scr);
    sesame__huber(scr, m, scr2, 1.5, 1e-6, &mIbG, &dummy, &mad0);
    if (mad0) s->status |= SESAME_STAT_NOOB_MAD0;
    alG = mIbG - muG; if (alG < 10.0) alG = 10.0;

    sesame__huber(bgR, nbR, scr, 1.5, 1e-6, &muR, &sgR, &mad0);
    if (mad0) s->status |= SESAME_STAT_NOOB_MAD0;
    m = sesame__drop_na(ibR, niR, scr);
    sesame__huber(scr, m, scr2, 1.5, 1e-6, &mIbR, &dummy, &mad0);
    if (mad0) s->status |= SESAME_STAT_NOOB_MAD0;
    alR = mIbR - muR; if (alR < 10.0) alR = 10.0;

    if (getenv("SESAME_DEBUG_NOOB"))
        fprintf(stderr, "NOOB fitG mu=%.10g s=%.10g a=%.10g | fitR mu=%.10g s=%.10g a=%.10g"
                " | nbG=%d nbR=%d niG=%d niR=%d\n",
                muG, sgG, alG, muR, sgR, alR, nbG, nbR, niG, niR);

    /* apply to every probe; NaN (Inf-II M, missing addresses) passes through */
    for (i = 0; i < n; i++) {
        s->MG[i] = sesame__norm_exp_signal(muG, sgG, alG, s->MG[i]) + offset;
        s->UG[i] = sesame__norm_exp_signal(muG, sgG, alG, s->UG[i]) + offset;
        s->MR[i] = sesame__norm_exp_signal(muR, sgR, alR, s->MR[i]) + offset;
        s->UR[i] = sesame__norm_exp_signal(muR, sgR, alR, s->UR[i]) + offset;
    }

out:
    free(bgG); free(bgR); free(ibG); free(ibR); free(scr); free(scr2);
    return rc;
}
