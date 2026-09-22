/* normexp_test -- unit harness exposing the two B-step numeric primitives so the
 * R oracle can compare them on identical inputs, isolating the D5/D6 arithmetic
 * from the mask-lineage divergence that the integration test carries.
 *
 *   normexp_test normexp   < "mu sigma alpha x" per line -> signal per line
 *   normexp_test huber      < one value per line (a vector) -> "mu<TAB>s<TAB>mad0"
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026-present The Children's Hospital of Philadelphia
 */
#include "internal.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: normexp_test normexp|huber\n"); return 2; }

    if (strcmp(argv[1], "normexp") == 0) {
        double mu, sigma, alpha, x;
        while (scanf("%lf %lf %lf %lf", &mu, &sigma, &alpha, &x) == 4)
            printf("%.17g\n", sesame__norm_exp_signal(mu, sigma, alpha, x));
        return 0;
    }

    if (strcmp(argv[1], "huber") == 0) {
        int32_t cap = 1024, n = 0, mad0;
        double *y = (double *)malloc((size_t)cap * sizeof(double)), *scr, v, mu, s;
        if (!y) return 1;
        while (scanf("%lf", &v) == 1) {
            if (n == cap) { cap *= 2; y = (double *)realloc(y, (size_t)cap * sizeof(double)); if (!y) return 1; }
            y[n++] = v;
        }
        scr = (double *)malloc((size_t)(n ? n : 1) * sizeof(double));
        if (!scr) return 1;
        sesame__huber(y, n, scr, 1.5, 1e-6, &mu, &s, &mad0);
        printf("%.17g\t%.17g\t%d\n", mu, s, mad0);
        free(y); free(scr);
        return 0;
    }

    fprintf(stderr, "unknown mode %s\n", argv[1]);
    return 2;
}
