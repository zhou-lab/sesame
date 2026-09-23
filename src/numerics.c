/* numerics.c -- the handful of R numerical primitives the prep steps rely on.
 *
 * Each one replicates a specific R semantic exactly; the comments say which,
 * because "close enough" here shows up as a different beta.
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
#include "internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    if (x < y) return -1;
    if (x > y) return 1;
    return 0;
}

/* Copy the non-NaN values of x into out (caller-allocated, >= n). Returns the
 * count. Mirrors R's na.rm=TRUE / sort()'s na.last=NA: NAs are dropped *before*
 * n is taken, so they never affect the quantile position. */
int32_t sesame__drop_na(const double *x, int32_t n, double *out)
{
    int32_t m = 0, i;
    for (i = 0; i < n; i++)
        if (!isnan(x[i])) out[m++] = x[i];
    return m;
}

/* R's quantile(x, p, type=7) on data already sorted ascending and NaN-free.
 *
 * Transcribed from stats:::quantile.default, which for type 7 does:
 *
 *     index <- 1 + max(n - 1, 0) * probs
 *     lo <- floor(index);  hi <- ceiling(index)
 *     qs <- x[lo]
 *     i  <- which(index > lo & x[hi] != qs)
 *     h  <- (index - lo)[i]
 *     qs[i] <- (1 - h) * qs[i] + h * x[hi[i]]
 *
 * The interpolation MUST be written `(1-h)*a + h*b`, not the algebraically
 * equivalent `a + h*(b-a)`: the two disagree in the last bit for ~24% of random
 * inputs, which showed up as ~1% of betas differing by 1-2 ULP. Likewise the
 * `x[hi] != qs` guard is load-bearing -- when the two order statistics are equal
 * R returns qs untouched rather than recombining them.
 *
 * Returns NaN for n == 0 (R errors; callers guard). */
double sesame__quantile7_sorted(const double *x, int32_t n, double p)
{
    double index, h, qs;
    int32_t lo, hi;

    if (n <= 0) return NAN;
    if (n == 1) return x[0];

    index = 1.0 + (double)(n - 1) * p;   /* R is 1-based */
    lo    = (int32_t)floor(index);
    hi    = (int32_t)ceil(index);
    if (lo < 1) lo = 1;
    if (hi > n) hi = n;

    qs = x[lo - 1];                       /* to 0-based */
    if ((double)lo < index && x[hi - 1] != qs) {
        h  = index - (double)lo;
        qs = (1.0 - h) * qs + h * x[hi - 1];
    }
    return qs;
}

/* Same, but sorts x in place first. */
double sesame__quantile7(double *x, int32_t n, double p)
{
    if (n <= 0) return NAN;
    qsort(x, (size_t)n, sizeof(double), cmp_double);
    return sesame__quantile7_sorted(x, n, p);
}

void sesame__sort(double *x, int32_t n)
{
    if (n > 1) qsort(x, (size_t)n, sizeof(double), cmp_double);
}

/* R's median() on sorted NaN-free data. Identical to quantile type 7 at p=0.5
 * (even n -> mean of the two middle order statistics), so it is not a separate
 * implementation. */
double sesame__median_sorted(const double *x, int32_t n)
{
    return sesame__quantile7_sorted(x, n, 0.5);
}

/* preprocessCore::normalize.quantiles.use.target(matrix(x), target).
 *
 * CLEAN-ROOM, characterized purely by black-box probing.
 *
 * Vendoring preprocessCore's qnorm.c would buy bit-exactness but is not open
 * to us: it is LGPL, and sesame-cli is BSD-2 with commercial licensing on
 * inquiry, so pulling it in would impose LGPL obligations on a work we
 * relicense. The clean-room stands on its own engineering merits anyway: the
 * clean-room already agrees to ~2 ULP, it is ~40 lines, and vendoring would pull
 * in an LGPL dependency and preprocessCore's LinkingTo: plumbing. See NUMERICS.md
 * (D8). Observations:
 *
 *   n==m, no ties  -> sorted(target)[rank(x)]
 *   n!=m           -> quantile(target, (rank-1)/(n-1), type=7)
 *   ties in x      -> quantile at the AVERAGE rank  (verified against the
 *                     competing "mean of the quantiles at each tied rank" rule
 *                     with a non-linear target: actual 100, rule-B 106.67)
 *   NA in target   -> dropped before ranking
 *   NA in x        -> preserved; ranking uses non-NA only
 *
 * This implementation requires x already sorted ascending and NaN-free, which
 * is how dyeBiasNL calls it (R/dye_bias.R:135). Because ties take the average
 * rank, every element of a tied run maps to the identical output -- confirmed
 * on real data (0 of 11259 tied runs on HM450 had non-constant output). That is
 * also why the subsequent approx(ties=mean) never actually averages anything.
 *
 * tgt is sorted in place. out must hold n doubles. Output is non-decreasing,
 * which is why R's subsequent sort() of the result is a no-op. */
void sesame__qnorm_use_target(const double *x_sorted, int32_t n,
                              double *tgt, int32_t m, double *out)
{
    int32_t a = 0;

    sesame__sort(tgt, m);
    while (a < n) {
        int32_t b = a, i;
        double r, p, q;
        while (b + 1 < n && x_sorted[b + 1] == x_sorted[a]) b++;
        /* average of the 1-based ranks spanned by the tied run [a, b] */
        r = ((double)(a + 1) + (double)(b + 1)) / 2.0;
        p = (n == 1) ? 0.0 : (r - 1.0) / (double)(n - 1);
        q = sesame__quantile7_sorted(tgt, m, p);
        for (i = a; i <= b; i++) out[i] = q;
        a = b + 1;
    }
}

/* R's pmax(a, b) with NO na.rm: NA in either operand propagates.
 *
 * C's fmax() does the opposite -- fmax(NaN, x) returns x -- which would
 * silently disagree with R at exactly the probes that matter (missing
 * addresses). Hence the explicit test. */
double sesame__pmax2(double a, double b)
{
    if (isnan(a) || isnan(b)) return NAN;
    return a > b ? a : b;
}

/* R's pmax(a, b, na.rm=TRUE): NA is dropped; both NA -> NA (verified against R,
 * which returns NA not -Inf). Used for the per-channel allele max in pOOBAH. */
double sesame__pmax2_narm(double a, double b)
{
    if (isnan(a)) return b;      /* b may also be NaN -> NaN, matching R */
    if (isnan(b)) return a;
    return a > b ? a : b;
}

/* Number of elements <= x in an ascending, NaN-free array (upper_bound). This
 * is k in R's ecdf: F(x) = k/n, so the pOOBAH detection p-value is 1 - k/n. */
int32_t sesame__count_le(const double *sorted, int32_t n, double x)
{
    int32_t lo = 0, hi = n;      /* first index with sorted[idx] > x */
    while (lo < hi) {
        int32_t mid = lo + ((hi - lo) >> 1);
        if (sorted[mid] <= x) lo = mid + 1;
        else                  hi = mid;
    }
    return lo;
}

/* Inverse Mills ratio lambda(t) = phi(t)/Phi(t) -- the conditional-expectation
 * factor in normal-exponential deconvolution (D5).
 *
 * R's normExpSignal forms this as exp(dnorm(...,log) - pnorm(...,log)), a
 * difference of logs that loses precision in the left tail (where the 1e-6 floor
 * then fires). We evaluate lambda directly:
 *
 *   - t > -5:  Phi(t) = 0.5*erfc(-t/sqrt2) is far above underflow, so
 *              phi(t)/Phi(t) is computed as written.
 *   - t <= -5: phi and Phi both underflow but lambda(t) -> -t. Use the Laplace
 *              continued fraction for the upper-tail Mills ratio
 *              m(s) = (1-Phi(s))/phi(s) = 1/(s + 1/(s + 2/(s + ...))), s = -t;
 *              then lambda = 1/m(s) = s + f, where f is the fraction after the
 *              leading s. Converges in a handful of terms for s >= 5.
 *
 * NaN in -> NaN out (the CF branch is unreachable for NaN since NaN <= -5 is
 * false, so NaN flows through the erfc branch). */
double sesame__inv_mills(double t)
{
    static const double INV_SQRT_2PI = 0.3989422804014326779399460599343819;
    static const double INV_SQRT_2   = 0.7071067811865475244008443621048490;

    if (t > -5.0) {
        double phi = INV_SQRT_2PI * exp(-0.5 * t * t);
        double Phi = 0.5 * erfc(-t * INV_SQRT_2);
        return phi / Phi;
    } else {
        double s = -t, f = 0.0;
        int j;
        for (j = 60; j >= 1; j--) f = (double)j / (s + f);
        return s + f;
    }
}

/* R's normExpSignal (R/background.R:143-167) written via the inverse Mills ratio:
 *   sigma2 = sigma^2;  mu.sf = x - mu - sigma2/alpha
 *   signal = mu.sf + sigma * lambda(mu.sf / sigma)
 * A non-NA signal that comes out negative is floored at 1e-6, exactly as R does.
 * The offset is NOT added here -- noob's caller adds it. NaN x -> NaN. */
double sesame__norm_exp_signal(double mu, double sigma, double alpha, double x)
{
    double sigma2 = sigma * sigma;
    double mu_sf  = x - mu - sigma2 / alpha;
    double signal = mu_sf + sigma * sesame__inv_mills(mu_sf / sigma);
    if (!isnan(signal) && signal < 0.0) signal = 1e-6;
    return signal;
}

/* Regularized incomplete beta I_x(a,b) -- the CDF workhorse for the t and F
 * distributions in DML. Lentz continued fraction (Numerical Recipes betacf) with
 * a lgamma prefactor; agrees with R's pbeta to ~1e-12 across the usual range.
 * Returns NaN for x outside [0,1]. */
static double betacf(double a, double b, double x)
{
    const int MAXIT = 300;
    const double EPS = 3e-16, FPMIN = 1e-300;
    double qab = a + b, qap = a + 1.0, qam = a - 1.0;
    double c = 1.0, d = 1.0 - qab * x / qap, aa, del;
    int m;
    if (fabs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    del = d;
    for (m = 1; m <= MAXIT; m++) {
        double m2 = 2.0 * m;
        aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d; if (fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c; if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d; del *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d; if (fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c; if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d; del *= (d * c);
        if (fabs(d * c - 1.0) <= EPS) break;
    }
    return del;
}

double sesame__betai(double a, double b, double x)
{
    double bt;
    if (isnan(x)) return NAN;
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    bt = exp(lgamma(a + b) - lgamma(a) - lgamma(b)
             + a * log(x) + b * log1p(-x));
    if (x < (a + 1.0) / (a + b + 2.0))
        return bt * betacf(a, b, x) / a;
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

/* Two-sided Student-t tail Pr(>|t|) = I_{df/(df+t^2)}(df/2, 1/2). */
double sesame__pt_2sided(double t, double df)
{
    if (isnan(t) || df <= 0.0) return NAN;
    return sesame__betai(0.5 * df, 0.5, df / (df + t * t));
}

/* Upper-tail F: P(F_{d1,d2} > f) = I_{d2/(d2+d1 f)}(d2/2, d1/2). */
double sesame__pf_upper(double f, double d1, double d2)
{
    if (isnan(f) || f < 0.0 || d1 <= 0.0 || d2 <= 0.0) return NAN;
    return sesame__betai(0.5 * d2, 0.5 * d1, d2 / (d2 + d1 * f));
}

/* MASS::huber(y, k, tol): the fixed-scale Huber M-estimator of location.
 *
 * mu starts at median(y); the scale s = mad(y) = 1.4826*median(|y-median(y)|) is
 * computed once and held fixed; the winsorized mean is then iterated until it
 * moves by less than tol*s. Note MASS returns the mu whose *successor* first fell
 * within tolerance -- i.e. the pre-update value -- so the break precedes the
 * assignment here too. R's sum() accumulates in long double, matched below.
 *
 * D6: MASS errors when mad == 0; we instead fall back to
 * s = max(IQR/1.349, 1e-8) and set *mad0, so a single degenerate array cannot
 * abort a whole-cohort run.
 *
 * y is sorted in place and must be NaN-free (caller drops NA); scratch holds >= n
 * doubles. */
void sesame__huber(double *y, int32_t n, double *scratch,
                   double k, double tol, double *mu_out, double *s_out,
                   int *mad0)
{
    double mu, s;
    int32_t i, iter;

    if (mad0) *mad0 = 0;
    sesame__sort(y, n);
    mu = sesame__median_sorted(y, n);
    for (i = 0; i < n; i++) scratch[i] = fabs(y[i] - mu);
    sesame__sort(scratch, n);
    s = 1.4826 * sesame__median_sorted(scratch, n);
    if (s == 0.0) {                                 /* D6 */
        double iqr = sesame__quantile7_sorted(y, n, 0.75)
                   - sesame__quantile7_sorted(y, n, 0.25);
        double alt = iqr / 1.349;
        s = alt > 1e-8 ? alt : 1e-8;
        if (mad0) *mad0 = 1;
    }
    for (iter = 0; iter < 1000; iter++) {           /* cap is a safety backstop */
        long double sum = 0.0L;
        double lo = mu - k * s, hi = mu + k * s, mu1;
        for (i = 0; i < n; i++) {
            double v = y[i];
            if (v < lo) v = lo; else if (v > hi) v = hi;
            sum += (long double)v;
        }
        mu1 = (double)(sum / (long double)n);
        if (fabs(mu - mu1) < tol * s) break;        /* return current mu, not mu1 */
        mu = mu1;
    }
    *mu_out = mu;
    *s_out  = s;
}
