/* mu2cg -- convert an M/U TSV (Probe_ID + <name>_M / <name>_U column pairs, in
 * ordering order) to a YAME format-3 .cg (+ .idx). Used to pack the CNV normal
 * reference exported from sesameData. Links libsesame.
 *
 *   mu2cg <in.mu.tsv[.gz]> <out.cg>
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026-present The Children's Hospital of Philadelphia
 */
#include "sesame.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>     /* gzopen/gzgets read the .tsv.gz (and plain .tsv) directly */

static long read_line(gzFile f, char **buf, size_t *cap)
{
    size_t len = 0;
    if (!*buf) { *cap = 1 << 16; *buf = (char *)malloc(*cap); if (!*buf) return -1; }
    for (;;) {
        if (len + 1 >= *cap) { char *n = realloc(*buf, *cap*2); if (!n) return -1; *buf = n; *cap *= 2; }
        if (!gzgets(f, *buf + len, (int)(*cap - len))) break;
        len += strlen(*buf + len);
        if (len && (*buf)[len-1] == '\n') { (*buf)[--len] = 0; return (long)len; }
        if (gzeof(f)) break;
    }
    return len ? (long)len : -1;
}

static int split_tabs(char *s, char **tok, int max)
{
    int n = 0;
    for (;;) { char *t = strchr(s, '\t'); if (n < max) tok[n] = s; n++; if (!t) break; *t = 0; s = t+1; }
    return n;
}

static double pn(const char *s)
{ char *e; double v; if (!s || !*s || !strcmp(s,"NA") || !strcmp(s,"NaN")) return NAN; v = strtod(s,&e); return e==s?NAN:v; }

int main(int argc, char **argv)
{
    gzFile f;
    char *line = NULL, **tok = NULL, **names = NULL;
    size_t cap = 0, caprow = 0;
    int maxtok, ncol, nsamp, s, rc = 1;
    int32_t np = 0, p;
    double *pmM = NULL, *pmU = NULL, *smM = NULL, *smU = NULL; /* probe-major, then sample-major */
    sesame_err_t e;

    if (argc < 3) { fprintf(stderr, "usage: mu2cg <in.mu.tsv[.gz]> <out.cg>\n"); return 2; }
    if (!(f = gzopen(argv[1], "rb"))) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    if (read_line(f, &line, &cap) < 0) { fprintf(stderr, "empty %s\n", argv[1]); return 1; }
    maxtok = (int)cap;
    tok = (char **)malloc((size_t)maxtok * sizeof(char *));
    ncol = split_tabs(line, tok, maxtok);
    if (ncol < 3 || (ncol - 1) % 2 != 0) { fprintf(stderr, "expected Probe_ID + M/U column pairs\n"); return 1; }
    nsamp = (ncol - 1) / 2;
    names = (char **)malloc((size_t)nsamp * sizeof(char *));
    for (s = 0; s < nsamp; s++) {              /* name = "<name>_M" minus "_M" */
        char *h = tok[1 + 2*s]; size_t l = strlen(h);
        if (l > 2 && !strcmp(h + l - 2, "_M")) h[l-2] = 0;
        names[s] = strdup(h);
    }

    while (read_line(f, &line, &cap) >= 0) {
        int nf = split_tabs(line, tok, maxtok);
        if (nf < 1 || !tok[0][0]) continue;
        if ((size_t)np >= caprow) {
            size_t nc = caprow ? caprow*2 : 8192;
            pmM = realloc(pmM, nc*(size_t)nsamp*sizeof(double));
            pmU = realloc(pmU, nc*(size_t)nsamp*sizeof(double));
            caprow = nc;
        }
        for (s = 0; s < nsamp; s++) {
            pmM[(size_t)np*(size_t)nsamp + (size_t)s] = (1+2*s   < nf) ? pn(tok[1+2*s])   : NAN;
            pmU[(size_t)np*(size_t)nsamp + (size_t)s] = (2+2*s   < nf) ? pn(tok[2+2*s])   : NAN;
        }
        np++;
    }
    gzclose(f);

    /* probe-major -> sample-major for sesame_write_cg_mu */
    smM = (double *)malloc((size_t)np*(size_t)nsamp*sizeof(double));
    smU = (double *)malloc((size_t)np*(size_t)nsamp*sizeof(double));
    for (p = 0; p < np; p++) for (s = 0; s < nsamp; s++) {
        smM[(size_t)s*(size_t)np + (size_t)p] = pmM[(size_t)p*(size_t)nsamp + (size_t)s];
        smU[(size_t)s*(size_t)np + (size_t)p] = pmU[(size_t)p*(size_t)nsamp + (size_t)s];
    }
    if (sesame_write_cg_mu(argv[2], smM, smU, np, nsamp, names, &e)) { fprintf(stderr, "%s\n", e.msg); goto out; }
    fprintf(stderr, "mu2cg: wrote %s (%d probes x %d samples, format 3) + %s.idx\n", argv[2], np, nsamp, argv[2]);
    rc = 0;
out:
    free(pmM); free(pmU); free(smM); free(smU); free(line); free(tok);
    if (names) { for (s = 0; s < nsamp; s++) free(names[s]); free(names); }
    return rc;
}
