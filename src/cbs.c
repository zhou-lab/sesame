/* cbs.c -- circular binary segmentation, a clean-room deterministic port of the
 * DNAcopy core (`fndcpt`/`changepoints`) sesame's cnSegmentation runs on the bin
 * log2 ratios. GPL(>=2) DNAcopy was read only to reproduce the algorithm; the
 * formulas here come from Olshen-Venkatraman (2004) and Siegmund (1988).
 *
 * DNAcopy assesses each candidate split by PERMUTATION against R's RNG, which
 * cnSegmentation never seeds -- so its borderline splits are not even reproducible
 * run-to-run. We keep the parts that ARE deterministic and well-defined:
 *   * the change-point location is the arc (i,j] maximising the CBS statistic
 *     BSS(i,j) = n (S_j-S_i)^2 / (L (n-L)), L = j-i -- exact, no RNG (a naive
 *     O(m^2) scan finds the same argmax as DNAcopy's block-decomposed tmaxo, as
 *     real-valued data has no ties);
 *   * significance is the Siegmund (1988) analytic tail probability `tailp` that
 *     DNAcopy itself uses for large segments (p.method="hybrid"), applied over the
 *     full valid arc range (delta = min.width/m) so it also covers the small arcs
 *     DNAcopy would permute -- deterministic, and a strict improvement in
 *     reproducibility;
 *   * the edge (1-vs-2 change-point) test uses the same two-sample statistic as
 *     `tpermp`, scored against the t distribution instead of a permutation;
 *   * the sdundo merge and trimmed-SD are already deterministic in DNAcopy.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present The Children's Hospital of Philadelphia
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "sesame.h"
#include "internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

static int cmp_i32(const void *a, const void *b)
{ int32_t x = *(const int32_t*)a, y = *(const int32_t*)b; return x<y?-1:(x>y?1:0); }

/* standard normal CDF (DNAcopy's fpnorm = pnorm) */
static double cbs_pnorm(double x) { return 0.5 * erfc(-x * M_SQRT1_2); }

/* Siegmund's nu(x): the correction function in the CBS tail probability. */
static double cbs_nu(double x, double tol)
{
    if (x > 0.01) {
        double lnu1 = log(2.0) - 2.0 * log(x), lnu0 = lnu1, dk = 0.0;
        int i, k = 2;
        for (i = 1; i <= k; i++) { dk += 1.0; lnu1 -= 2.0 * cbs_pnorm(-x * sqrt(dk) / 2.0) / dk; }
        while (fabs((lnu1 - lnu0) / lnu1) > tol) {
            lnu0 = lnu1;
            for (i = 1; i <= k; i++) { dk += 1.0; lnu1 -= 2.0 * cbs_pnorm(-x * sqrt(dk) / 2.0) / dk; }
            k *= 2;
        }
        return exp(lnu1);
    }
    return exp(-0.583 * x);
}

/* integral of 1/(t(1-t))^2 from x to x+a (closed form). */
static double cbs_it1tsq(double x, double a)
{
    double y = x + a - 0.5, r;
    r = (8.0 * y) / (1.0 - 4.0 * y * y) + 2.0 * log((1.0 + 2.0 * y) / (1.0 - 2.0 * y));
    y = x - 0.5;
    r -= (8.0 * y) / (1.0 - 4.0 * y * y) + 2.0 * log((1.0 + 2.0 * y) / (1.0 - 2.0 * y));
    return r;
}

/* Two-sided analytic tail probability P(max CBS statistic > b) for a segment of
 * length m, integrating the arc fraction t over [delta, 1-delta]. */
static double cbs_tailp(double b, double delta, int32_t m, int ngrid, double tol)
{
    double dincr = (0.5 - delta) / ngrid, bsqrtm = b / sqrt((double)m);
    double tl = 0.5 - dincr, t = 0.5 - 0.5 * dincr, tp = 0.0, x, nux;
    int i;
    if (dincr <= 0.0) return 1.0;                 /* degenerate: no valid arc */
    for (i = 0; i < ngrid; i++) {
        tl += dincr; t += dincr;
        x = bsqrtm / sqrt(t * (1.0 - t));
        nux = cbs_nu(x, tol);
        tp += nux * nux * cbs_it1tsq(tl, dincr);
    }
    tp = 9.973557e-2 * (b * b * b) * exp(-b * b / 2.0) * tp;
    return 2.0 * tp;
}

/* Maximise BSS(i,j) = n (S_j-S_i)^2 / (L(n-L)) over arcs with al0 <= L <= n-al0,
 * then return the F-standardised statistic BSS / ((tss-BSS)/(n-2)) -- exactly
 * DNAcopy tmaxo's output, so sqrt(ostat) is a t-scale value the 0.1/7.0 cut-offs
 * and tailp expect. x is centred (mean ~0); *pi,*pj get the arc (0..n). */
static double cbs_tmax(const double *x, int32_t n, int32_t al0, int32_t *pi, int32_t *pj)
{
    double *S = (double *)malloc((size_t)(n + 1) * sizeof(double));
    double rn = (double)n, best = 0.0, tss = 0.0;
    int32_t i, j, bi = 0, bj = n;
    if (!S) { *pi = 0; *pj = n; return 0.0; }
    S[0] = 0.0;
    for (i = 1; i <= n; i++) { S[i] = S[i-1] + x[i-1]; tss += x[i-1] * x[i-1]; }
    for (i = 0; i <= n - al0; i++) {
        int32_t jmax = i + (n - al0);
        if (jmax > n) jmax = n;
        for (j = i + al0; j <= jmax; j++) {
            int32_t L = j - i;
            double d = S[j] - S[i];
            double bss = rn * d * d / ((double)L * (double)(n - L));
            if (bss > best) { best = bss; bi = i; bj = j; }
        }
    }
    free(S);
    *pi = bi; *pj = bj;
    if (tss <= best + 1e-4) tss = best + 1.0;         /* DNAcopy guard */
    return best / ((tss - best) / (rn - 2.0));        /* F-standardise */
}

/* DNAcopy's tpermp statistic, scored analytically: is there a change at index n1
 * within x[0..n-1] (n = n1 + n2)? Returns a two-sided p-value. */
static double cbs_edge_p(const double *x, int32_t n1, int32_t n2)
{
    int32_t n = n1 + n2, i, m1;
    double xsum1 = 0, xsum2 = 0, tss = 0, xbar, ostat, tstat, rn = n, rn1 = n1, rn2 = n2, rm1;
    if (n1 <= 1 || n2 <= 1) return 1.0;
    for (i = 0; i < n1; i++) { xsum1 += x[i]; tss += x[i]*x[i]; }
    for (i = n1; i < n; i++) { xsum2 += x[i]; tss += x[i]*x[i]; }
    xbar = (xsum1 + xsum2) / rn;
    tss -= rn * xbar * xbar;
    if (n1 <= n2) { m1 = n1; rm1 = rn1; ostat = fabs(xsum1/rn1 - xbar); tstat = ostat*ostat*rn1*rn/rn2; }
    else          { m1 = n2; rm1 = rn2; ostat = fabs(xsum2/rn2 - xbar); tstat = ostat*ostat*rn2*rn/rn1; }
    (void)rm1;
    if (tss <= tstat) return 0.0;                 /* perfect split */
    tstat = tstat / ((tss - tstat) / (rn - 2.0)); /* two-sample t^2 (F_{1,n-2}) */
    if (m1 < 2) return 1.0;
    return sesame__pf_upper(tstat, 1.0, rn - 2.0);
}

/* trimmed-variance SD estimate (DNAcopy trimmed.variance), the sdundo scale.
 * inflfact is DNAcopy's variance-inflation correction; it depends only on trim,
 * and cnSegmentation always uses trim=0.025, so it is the compiled-in constant
 * (= inflfact(0.025), verified against DNAcopy). */
double sesame__trimmed_sd(const double *x, int32_t n, double trim)
{
    const double inflfact = 1.31779804080133;     /* DNAcopy inflfact(0.025) */
    int32_t nkeep, i;
    double *d, s = 0.0;
    if (n < 3) return 0.0;
    nkeep = (int32_t)lround((1.0 - 2.0*trim) * (n - 1));
    if (nkeep < 1) nkeep = 1;
    d = (double *)malloc((size_t)(n - 1) * sizeof(double));
    if (!d) return 0.0;
    for (i = 0; i < n - 1; i++) d[i] = fabs(x[i+1] - x[i]);
    sesame__sort(d, n - 1);
    for (i = 0; i < nkeep; i++) s += d[i]*d[i];
    free(d);
    return sqrt(inflfact * s / (2.0 * nkeep));
}

/* median of x[a..b] (inclusive), non-destructive. */
static double seg_median(const double *x, int32_t a, int32_t b)
{
    int32_t k = b - a + 1, i;
    double *t = (double *)malloc((size_t)k * sizeof(double)), m;
    if (!t) return NAN;
    for (i = 0; i < k; i++) t[i] = x[a + i];
    sesame__sort(t, k);
    m = sesame__median_sorted(t, k);
    free(t);
    return m;
}

/* sdundo: while the smallest |median diff| of adjacent segments is below
 * undo_sd * trimmed_sd, drop that boundary. seg_end holds cumulative ends. */
static int32_t cbs_sdundo(const double *x, int32_t *seg_end, int32_t nseg,
                          double trimmed_sd, double undo_sd)
{
    double thr = trimmed_sd * undo_sd;
    while (nseg > 1) {
        int32_t i, imin = -1, lo, hi;
        double dmin = 0.0, prev_med, med;
        lo = 0; hi = seg_end[0] - 1;
        prev_med = seg_median(x, lo, hi);
        for (i = 1; i < nseg; i++) {
            lo = seg_end[i-1]; hi = seg_end[i] - 1;
            med = seg_median(x, lo, hi);
            { double ad = fabs(med - prev_med);
              if (imin < 0 || ad < dmin) { dmin = ad; imin = i; } }
            prev_med = med;
        }
        if (imin < 0 || dmin >= thr) break;
        memmove(&seg_end[imin-1], &seg_end[imin], (size_t)(nseg - imin) * sizeof(int32_t));
        nseg--;
    }
    return nseg;
}

/* Segment x[0..n-1]; fills seg_end[] (cumulative 1-based ends) and returns nseg.
 * Deterministic CBS: recursive max-BSS split gated by the Siegmund tail. */
int32_t sesame__cbs(const double *x, int32_t n, int32_t min_width, double alpha,
                    int32_t kmax, int32_t nmin, int ngrid, double tol,
                    double trimmed_sd, double undo_sd,
                    int32_t *seg_end, int32_t max_seg)
{
    int32_t *stk = (int32_t *)malloc((size_t)(n + 2) * sizeof(int32_t));
    int32_t nstk = 2, ncl = 0, nseg;
    double *sub = (double *)malloc((size_t)n * sizeof(double));
    (void)kmax;
    if (!stk || !sub) { free(stk); free(sub); if (max_seg > 0) { seg_end[0] = n; } return n > 0; }
    stk[0] = 0; stk[1] = n;                        /* stack of segment ends: {0, n} */

    while (nstk > 1) {
        int32_t lo = stk[nstk-2], hi = stk[nstk-1], m = hi - lo, i, j, ncpt = 0, ic1 = 0, ic2 = 0;
        double mean = 0.0, mn, mx, ostat, ostat1, delta;

        if (m >= 2 * min_width) {
            for (i = 0; i < m; i++) mean += x[lo+i];
            mean /= m;
            mn = mx = x[lo];
            for (i = 0; i < m; i++) { double v = x[lo+i]; sub[i] = v - mean; if (v < mn) mn = v; if (v > mx) mx = v; }
            if (mx > mn) {
                ostat = cbs_tmax(sub, m, min_width, &i, &j);
                ostat1 = sqrt(ostat);
                if (ostat1 > 0.1) {
                    int32_t arclen = j - i, l = arclen < m - arclen ? arclen : m - arclen;
                    int split;
                    /* Siegmund tail over the full valid arc range [min_width,
                     * m-min_width]; covers the small arcs DNAcopy would permute. */
                    delta = (double)min_width / (double)m;
                    if (ostat1 >= 7.0 && l >= 10) split = 1;
                    else split = (cbs_tailp(ostat1, delta, m, ngrid, tol) <= alpha);
                    if (split) {
                        if (j == m) { ncpt = 1; ic1 = i; }        /* arc touches right end */
                        else if (i == 0) { ncpt = 1; ic1 = j; }   /* arc touches left end  */
                        else {                                    /* interior bump: test both edges */
                            if (cbs_edge_p(sub, i, j - i) <= alpha) { ic1 = i; ncpt = 1; }
                            if (cbs_edge_p(sub + i, j - i, m - j) <= alpha) { if (ncpt == 0) ic1 = j; else ic2 = j; ncpt++; }
                        }
                    }
                }
            }
        }
        (void)nmin;
        if (ncpt == 0) {
            if (ncl < max_seg) seg_end[ncl] = hi;   /* finalize this segment's end */
            ncl++;
            nstk--;                                  /* pop */
        } else if (ncpt == 1) {
            stk[nstk-1] = lo + ic1;                  /* [lo, lo+ic1] then [lo+ic1, hi] */
            stk[nstk] = hi;
            nstk++;
        } else {                                     /* two change points */
            stk[nstk-1] = lo + ic1;
            stk[nstk] = lo + ic2;
            stk[nstk+1] = hi;
            nstk += 2;
        }
    }
    free(stk); free(sub);

    /* ends were finalized right-to-left; sort ascending into segment order */
    nseg = ncl < max_seg ? ncl : max_seg;
    qsort(seg_end, (size_t)nseg, sizeof(int32_t), cmp_i32);
    if (trimmed_sd > 0.0 && undo_sd > 0.0 && nseg > 1)
        nseg = cbs_sdundo(x, seg_end, nseg, trimmed_sd, undo_sd);
    return nseg;
}
