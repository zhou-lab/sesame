/* sigdf.c -- IDAT pair + index -> signals -> betas.
 *
 * Ports chipAddressToSignal (R/sesame.R:463-508) and getBetas
 * (R/sesame.R:192-225).
 *
 * The address->signal mapping, which is easy to get wrong:
 *
 *   Infinium-I (col == G or R): TWO address lookups.
 *       M address -> MG (green) and MR (red)
 *       U address -> UG (green) and UR (red)
 *     Both channels are recorded at both addresses; `col` only says which is
 *     in-band. For col==G, MG/UG are in-band and MR/UR are the out-of-band
 *     (oobR) signal; for col==R it is the reverse. The Grn/Red split is NOT
 *     applied at read time.
 *
 *   Infinium-II (col == 2): ONE address lookup, on the U address only.
 *       MG = MR = NA
 *       UG = methylated allele, UR = unmethylated allele
 *     The manifest's M address is NA and is ignored.
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

/* IDAT addresses are sorted ascending and unique (verified across every test
 * array on EPIC/EPIC+/EPICv2/HM27/HM450/Mammal40/MM285/MSA), so a binary
 * search beats building a hash. Returns -1 when absent. */
static int32_t addr_find(const uint32_t *a, int32_t n, uint32_t key)
{
    int32_t lo = 0, hi = n - 1;
    while (lo <= hi) {
        int32_t mid = lo + ((hi - lo) >> 1);
        if (a[mid] == key) return mid;
        if (a[mid] < key) lo = mid + 1;
        else              hi = mid - 1;
    }
    return -1;
}

sesame_sigdf_t *sesame_sigdf_from_idats(const sesame_idat_t *grn,
                                        const sesame_idat_t *red,
                                        const sesame_index_t *ix,
                                        int min_beads,
                                        sesame_err_t *err)
{
    sesame_sigdf_t *s = NULL;
    const uint32_t *M, *U;
    const uint8_t *col, *dmask;
    int32_t n, i;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!grn || !red || !ix) {
        sesame__fail(err, SESAME_ERR_IO, "null argument");
        return NULL;
    }

    /* D1: R cbinds Grn and Red positionally, taking rownames from Grn only
     * (R/sesame.R:293-298) -- no join on IlluminaID. The assumption holds for
     * well-formed files, but if a caller pairs mismatched IDATs (wrong sample
     * or chip) R silently emits garbage betas array-wide. Verify instead. */
    if (grn->n != red->n ||
        memcmp(grn->addr, red->addr, (size_t)grn->n * sizeof(uint32_t)) != 0) {
        sesame__fail(err, SESAME_ERR_FORMAT,
            "Grn/Red address vectors differ (n=%d vs %d) -- mismatched IDAT pair?",
            grn->n, red->n);
        return NULL;
    }

    n     = sesame_index_nprobes(ix);
    M     = sesame__index_M(ix);
    U     = sesame__index_U(ix);
    col   = sesame__index_col(ix);
    dmask = sesame__index_mask(ix);
    (void)dmask;  /* deliberately unused -- see below */

    s = (sesame_sigdf_t *)calloc(1, sizeof(*s));
    if (!s) { sesame__fail(err, SESAME_ERR_NOMEM, "oom"); return NULL; }
    s->ix = ix;
    s->n  = n;
    s->MG = (double *)malloc((size_t)n * sizeof(double));
    s->MR = (double *)malloc((size_t)n * sizeof(double));
    s->UG = (double *)malloc((size_t)n * sizeof(double));
    s->UR = (double *)malloc((size_t)n * sizeof(double));
    s->col  = (uint8_t *)malloc((size_t)n);
    s->mask = (uint8_t *)malloc((size_t)n);
    if (!s->MG || !s->MR || !s->UG || !s->UR || !s->col || !s->mask) {
        sesame_sigdf_free(s);
        sesame__fail(err, SESAME_ERR_NOMEM, "oom sizing %d probes", n);
        return NULL;
    }

    for (i = 0; i < n; i++) {
        int32_t im, iu;
        int miss = 0, low = 0;

        s->col[i]  = col[i];
        /* Start unmasked. The ordering table carries a design mask column
         * (create_default_mask()$ref_issue, written by sesameAnno.R:170), but
         * chipAddressToSignal explicitly sets mask=FALSE (R/sesame.R:475,498)
         * and never reads it -- the column is written but never consumed
         * anywhere in the pipeline. Seeding from it here would mask ~257
         * (HM450) / ~2263 (EPICv2) / ~3120 (MSA) probes that R keeps. Design
         * masking is qualityMask's job (the "Q" step). */
        s->mask[i] = 0;

        if (col[i] == SESAME_COL_II) {
            /* One lookup, U address only. */
            iu = (U[i] == 0) ? -1 : addr_find(grn->addr, grn->n, U[i]);
            s->MG[i] = s->MR[i] = NAN;
            if (iu < 0) {
                s->UG[i] = s->UR[i] = NAN;
                miss = 1;
            } else {
                s->UG[i] = (double)grn->mean[iu];
                s->UR[i] = (double)red->mean[iu];
                if (min_beads > 0 &&
                    (grn->nbeads[iu] < min_beads || red->nbeads[iu] < min_beads))
                    low = 1;
            }
        } else {
            /* Two lookups: M address -> MG/MR, U address -> UG/UR. */
            im = (M[i] == 0) ? -1 : addr_find(grn->addr, grn->n, M[i]);
            iu = (U[i] == 0) ? -1 : addr_find(grn->addr, grn->n, U[i]);
            if (im < 0) { s->MG[i] = s->MR[i] = NAN; miss = 1; }
            else {
                s->MG[i] = (double)grn->mean[im];
                s->MR[i] = (double)red->mean[im];
                if (min_beads > 0 &&
                    (grn->nbeads[im] < min_beads || red->nbeads[im] < min_beads))
                    low = 1;
            }
            if (iu < 0) { s->UG[i] = s->UR[i] = NAN; miss = 1; }
            else {
                s->UG[i] = (double)grn->mean[iu];
                s->UR[i] = (double)red->mean[iu];
                if (min_beads > 0 &&
                    (grn->nbeads[iu] < min_beads || red->nbeads[iu] < min_beads))
                    low = 1;
            }
        }

        if (miss) {
            s->n_addr_missing++;
            s->status |= SESAME_STAT_ADDR_MISSING;
            /* R yields all-NA rows here and only masks them when min_beads is
             * set (the NA bead count trips the < min_beads test). Match that:
             * mask on missing only when bead filtering is on. */
            if (min_beads > 0) s->mask[i] = 1;
        }
        if (low) s->mask[i] = 1;
    }

    return s;
}

void sesame_sigdf_free(sesame_sigdf_t *s)
{
    if (!s) return;
    free(s->MG); free(s->MR); free(s->UG); free(s->UR);
    free(s->col); free(s->mask);
    free(s);
}

/* Deep copy (shares the read-only index). Used to keep a raw SigDF alongside a
 * preprocessed one. Returns NULL on OOM. */
sesame_sigdf_t *sesame_sigdf_dup(const sesame_sigdf_t *s)
{
    sesame_sigdf_t *d;
    size_t z;
    if (!s) return NULL;
    d = (sesame_sigdf_t *)calloc(1, sizeof *d);
    if (!d) return NULL;
    z = (size_t)s->n;
    d->ix = s->ix; d->n = s->n; d->status = s->status;
    d->n_addr_missing = s->n_addr_missing;
    d->MG = (double *)malloc(z * sizeof(double));
    d->MR = (double *)malloc(z * sizeof(double));
    d->UG = (double *)malloc(z * sizeof(double));
    d->UR = (double *)malloc(z * sizeof(double));
    d->col  = (uint8_t *)malloc(z);
    d->mask = (uint8_t *)malloc(z);
    if (!d->MG || !d->MR || !d->UG || !d->UR || !d->col || !d->mask) {
        sesame_sigdf_free(d); return NULL;
    }
    memcpy(d->MG, s->MG, z * sizeof(double));
    memcpy(d->MR, s->MR, z * sizeof(double));
    memcpy(d->UG, s->UG, z * sizeof(double));
    memcpy(d->UR, s->UR, z * sizeof(double));
    memcpy(d->col, s->col, z);
    memcpy(d->mask, s->mask, z);
    return d;
}

/* getBetas (R/sesame.R:192-225).
 *
 *   Inf-I  G : max(MG,1) / max(MG+UG, 2)
 *   Inf-I  R : max(MR,1) / max(MR+UR, 2)
 *   Inf-II   : max(UG,1) / max(UG+UR, 2)
 *
 * Note sesame uses NO Illumina +100 offset. R's pmax() here has no na.rm, so
 * an NA in either allele propagates to an NA beta -- replicated exactly. */
int sesame_get_betas(const sesame_sigdf_t *s, int apply_mask,
                     double *out, sesame_err_t *err)
{
    int32_t i;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!s || !out) return sesame__fail(err, SESAME_ERR_IO, "null argument");

    for (i = 0; i < s->n; i++) {
        double m, u;

        switch (s->col[i]) {
        case SESAME_COL_G:  m = s->MG[i]; u = s->UG[i]; break;
        case SESAME_COL_R:  m = s->MR[i]; u = s->UR[i]; break;
        default:            m = s->UG[i]; u = s->UR[i]; break; /* Inf-II */
        }

        if (isnan(m) || isnan(u)) {
            out[i] = NAN;
        } else {
            double num = m > 1.0 ? m : 1.0;
            double den = (m + u) > 2.0 ? (m + u) : 2.0;
            out[i] = num / den;
        }

        if (apply_mask && s->mask[i]) out[i] = NAN;
    }
    return SESAME_OK;
}

/* Per-probe methylated/unmethylated signal (R's signalMU): col G -> (MG,UG),
 * col R -> (MR,UR), Inf-II -> (UG,UR). NaN preserved. Either out may be NULL. */
int sesame_signal_mu(const sesame_sigdf_t *s, double *M, double *U,
                     sesame_err_t *err)
{
    int32_t i;
    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!s || (!M && !U)) return sesame__fail(err, SESAME_ERR_IO, "null argument");
    for (i = 0; i < s->n; i++) {
        double m, u;
        switch (s->col[i]) {
        case SESAME_COL_G: m = s->MG[i]; u = s->UG[i]; break;
        case SESAME_COL_R: m = s->MR[i]; u = s->UR[i]; break;
        default:           m = s->UG[i]; u = s->UR[i]; break;
        }
        if (M) M[i] = m;
        if (U) U[i] = u;
    }
    return SESAME_OK;
}

/* Total intensity M+U per probe (R's totalIntensities, mask=FALSE): col G ->
 * MG+UG, col R -> MR+UR, Inf-II -> UG+UR. NaN if either allele is NA. This is
 * the CNV signal input and the quantity summarized by the QC intensity group. */
int sesame_total_intensities(const sesame_sigdf_t *s, double *out, sesame_err_t *err)
{
    int32_t i;
    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!s || !out) return sesame__fail(err, SESAME_ERR_IO, "null argument");
    for (i = 0; i < s->n; i++) {
        double m, u;
        switch (s->col[i]) {
        case SESAME_COL_G: m = s->MG[i]; u = s->UG[i]; break;
        case SESAME_COL_R: m = s->MR[i]; u = s->UR[i]; break;
        default:           m = s->UG[i]; u = s->UR[i]; break;
        }
        out[i] = (isnan(m) || isnan(u)) ? NAN : m + u;
    }
    return SESAME_OK;
}
