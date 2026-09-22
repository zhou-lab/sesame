/* cbs_test -- read a bin-signal vector (one double per line) from stdin, run the
 * deterministic CBS with sesame's cnSegmentation parameters, and print the
 * segments as "start end nmark mean" (1-based, inclusive). Used to validate
 * sesame__cbs against DNAcopy::segment on identical input.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026-present The Children's Hospital of Philadelphia
 */
#include "sesame.h"
#include "internal.h"
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    double *x = NULL, v;
    int32_t n = 0, cap = 0, *seg_end, i, nseg, lo;
    double tsd;

    while (scanf("%lf", &v) == 1) {
        if (n >= cap) { cap = cap ? cap*2 : 1024; x = realloc(x, (size_t)cap*sizeof(double)); }
        x[n++] = v;
    }
    if (n == 0) return 0;
    seg_end = malloc((size_t)n * sizeof(int32_t));

    tsd = sesame__trimmed_sd(x, n, 0.025);
    /* sesame's segmentBins params: min.width=5, alpha=0.001, undo.SD=2.2 */
    nseg = sesame__cbs(x, n, 5, 0.001, 25, 200, 100, 1e-6, tsd, 2.2, seg_end, n);

    lo = 0;
    for (i = 0; i < nseg; i++) {
        int32_t hi = seg_end[i], k;
        double s = 0.0;
        for (k = lo; k < hi; k++) s += x[k];
        printf("%d\t%d\t%d\t%.6g\n", lo + 1, hi, hi - lo, s / (hi - lo));
        lo = hi;
    }
    free(x); free(seg_end);
    return 0;
}
