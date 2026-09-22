/* cnv.c -- copy-number signal, mirroring sesame's cnSegmentation up to the bin
 * log2 ratios (CBS segmentation is left to a downstream tool).
 *
 * For one target sample: regress its per-probe total intensity on a panel of
 * normal samples' total intensities (OLS with intercept, the shared Householder
 * QR in dml.c), then log2(target / max(fitted, 1)) per probe -- the copy-number
 * signal. Probes are binned along the genome exactly as R does: tile each
 * chromosome at `tilewidth`, subtract assembly gaps, left/right-merge bins until
 * each holds >= min_probes, and take the median probe signal per bin.
 *
 * Target, normals, and the per-probe coordinate table are all positional in the
 * ordering, so probe i is the same row across all three -- no id join needed.
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
#include <zlib.h>

/* ---------------------------------------------------------- text load ---- */

static int read_gzline(gzFile f, char **buf, size_t *cap)
{
    size_t len;
    if (!gzgets(f, *buf, (int)*cap)) return 0;
    len = strlen(*buf);
    while (len + 1 == *cap && (*buf)[len-1] != '\n') {
        char *n = (char *)realloc(*buf, *cap * 2);
        if (!n) return -1;
        *buf = n; *cap *= 2;
        if (!gzgets(f, *buf + len, (int)(*cap - len))) break;
        len += strlen(*buf + len);
    }
    (*buf)[strcspn(*buf, "\r\n")] = '\0';
    return 1;
}

/* Per-probe coordinates, positional in the ordering. chrom[i] is a strdup'd
 * chromosome ("" if unmapped), pos[i] the start (or -1). */
static int load_coords(const char *path, int32_t np, char ***chrom_out,
                       int32_t **pos_out, sesame_err_t *err)
{
    gzFile f = gzopen(path, "rb");
    char *buf, *tab, *tab2;
    size_t cap = 1 << 16;
    char **chrom = NULL;
    int32_t *pos = NULL, row = 0, r;

    if (!f) return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", path);
    if (!(buf = (char *)malloc(cap))) { gzclose(f);
        return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }
    chrom = (char **)malloc((size_t)np * sizeof(char *));
    pos = (int32_t *)malloc((size_t)np * sizeof(int32_t));
    if (!chrom || !pos) { free(buf); free(chrom); free(pos); gzclose(f);
        return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }

    r = read_gzline(f, &buf, &cap);              /* header */
    while ((r = read_gzline(f, &buf, &cap)) == 1) {
        if (row >= np) { row++; continue; }      /* count overflow, report below */
        tab = strchr(buf, '\t');
        if (tab) *tab = '\0';
        if (buf[0] == '\0' || !strcmp(buf, "*") || !strcmp(buf, "NA")) {
            chrom[row] = strdup(""); pos[row] = -1;
        } else {
            chrom[row] = strdup(buf);
            tab2 = tab ? strchr(tab + 1, '\t') : NULL;
            if (tab2) *tab2 = '\0';
            pos[row] = tab ? (int32_t)strtol(tab + 1, NULL, 10) : -1;
        }
        row++;
    }
    free(buf); gzclose(f);
    if (r < 0 || row != np) {
        for (int32_t i = 0; i < row && i < np; i++) free(chrom[i]);
        free(chrom); free(pos);
        if (r < 0) return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
        return sesame__fail(err, SESAME_ERR_FORMAT,
            "%s has %d data rows, ordering has %d -- lineage mismatch", path, row, np);
    }
    *chrom_out = chrom; *pos_out = pos;
    return SESAME_OK;
}

/* seqinfo.tsv.gz (chrom<TAB>length) and gaps.tsv.gz (chrom<TAB>start<TAB>end).
 * Both keep chromosome order as written (seqinfo order drives the output). */
typedef struct { int32_t n; char **chrom; long *len; } seqinfo_t;
typedef struct { int32_t n; char **chrom; long *start; long *end; } gaps_t;

static int load_seqinfo(const char *path, seqinfo_t *s, sesame_err_t *err)
{
    gzFile f = gzopen(path, "rb");
    char *buf, *tab; size_t cap = 1 << 16;
    int32_t cap_n = 0, r;
    memset(s, 0, sizeof *s);
    if (!f) return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", path);
    if (!(buf = (char *)malloc(cap))) { gzclose(f); return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }
    read_gzline(f, &buf, &cap);                  /* header */
    while ((r = read_gzline(f, &buf, &cap)) == 1) {
        if (!buf[0]) continue;
        tab = strchr(buf, '\t'); if (!tab) continue; *tab = '\0';
        if (s->n >= cap_n) { cap_n = cap_n ? cap_n*2 : 64;
            s->chrom = realloc(s->chrom, (size_t)cap_n*sizeof(char*));
            s->len = realloc(s->len, (size_t)cap_n*sizeof(long)); }
        s->chrom[s->n] = strdup(buf);
        s->len[s->n] = strtol(tab+1, NULL, 10);
        s->n++;
    }
    free(buf); gzclose(f);
    if (r < 0) return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
    return SESAME_OK;
}

static int load_gaps(const char *path, gaps_t *g, sesame_err_t *err)
{
    gzFile f = gzopen(path, "rb");
    char *buf, *t1, *t2; size_t cap = 1 << 16;
    int32_t cap_n = 0, r;
    memset(g, 0, sizeof *g);
    if (!f) return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", path);
    if (!(buf = (char *)malloc(cap))) { gzclose(f); return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }
    read_gzline(f, &buf, &cap);                  /* header */
    while ((r = read_gzline(f, &buf, &cap)) == 1) {
        if (!buf[0]) continue;
        t1 = strchr(buf, '\t'); if (!t1) continue; *t1 = '\0';
        t2 = strchr(t1+1, '\t'); if (!t2) continue; *t2 = '\0';
        if (g->n >= cap_n) { cap_n = cap_n ? cap_n*2 : 256;
            g->chrom = realloc(g->chrom, (size_t)cap_n*sizeof(char*));
            g->start = realloc(g->start, (size_t)cap_n*sizeof(long));
            g->end = realloc(g->end, (size_t)cap_n*sizeof(long)); }
        g->chrom[g->n] = strdup(buf);
        g->start[g->n] = strtol(t1+1, NULL, 10);
        g->end[g->n] = strtol(t2+1, NULL, 10);
        g->n++;
    }
    free(buf); gzclose(f);
    if (r < 0) return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
    return SESAME_OK;
}

static void free_seqinfo(seqinfo_t *s) { for (int32_t i=0;i<s->n;i++) free(s->chrom[i]); free(s->chrom); free(s->len); }
static void free_gaps(gaps_t *g) { for (int32_t i=0;i<g->n;i++) free(g->chrom[i]); free(g->chrom); free(g->start); free(g->end); }

/* -------------------------------------------------------------- binning -- */

typedef struct { long start, end; double *sig; int32_t n, cap; } bin_t;

static void bin_push(bin_t *b, double v)
{
    if (b->n >= b->cap) { b->cap = b->cap ? b->cap*2 : 8;
        b->sig = realloc(b->sig, (size_t)b->cap*sizeof(double)); }
    b->sig[b->n++] = v;
}

/* One probe on this chromosome (sorted by pos). */
typedef struct { long pos; double sig; } cprobe_t;
static int cprobe_cmp(const void *a, const void *b)
{ long x = ((const cprobe_t*)a)->pos, y = ((const cprobe_t*)b)->pos; return x<y?-1:(x>y?1:0); }

/* Append this chromosome's merged bins to the growing output arrays. gs/gn are
 * this chrom's gaps (sorted by start). Mirrors getBinCoordinates +
 * leftRightMerge1: tile at tilewidth, subtract gaps, merge to >= min_probes. */
static void bin_chrom(const char *chrom, long chrlen, cprobe_t *pr, int32_t npr,
                      const long *gstart, const long *gend, int32_t gn,
                      int32_t tilewidth, int32_t min_probes,
                      char ***bc, int32_t **bs, int32_t **be, int32_t **bnp,
                      double **bsig, int32_t *nbin, int32_t *capbin)
{
    bin_t *bins = NULL; int32_t nb = 0, cb = 0;
    int32_t gi = 0, pi = 0;
    long t;

    qsort(pr, (size_t)npr, sizeof *pr, cprobe_cmp);

    /* tiles [t+1, min(t+tw, chrlen)] (1-based inclusive), each minus gaps */
    for (t = 0; t < chrlen; t += tilewidth) {
        long ts = t + 1, te = t + tilewidth; if (te > chrlen) te = chrlen;
        long cur = ts;
        /* advance gap cursor to the first gap that could overlap this tile */
        while (gi > 0 && gend[gi-1] >= ts) gi--;      /* rewind if needed */
        while (gi < gn && gend[gi] < ts) gi++;
        int32_t gj = gi;
        for (;;) {
            long fs = cur, fe = te;                    /* candidate fragment */
            /* clip against the next overlapping gap */
            while (gj < gn && gend[gj] < cur) gj++;
            if (gj < gn && gstart[gj] <= te) {         /* a gap intrudes */
                if (gstart[gj] > cur) fe = gstart[gj] - 1;  /* fragment before gap */
                else {                                  /* gap covers cur: skip past */
                    cur = gend[gj] + 1; gj++;
                    if (cur > te) break;
                    continue;
                }
            }
            if (fe >= fs) {
                bin_t b; memset(&b, 0, sizeof b); b.start = fs; b.end = fe;
                /* probes in [fs, fe] (pr sorted; pi is a forward cursor) */
                while (pi < npr && pr[pi].pos < fs) pi++;
                while (pi < npr && pr[pi].pos <= fe) { bin_push(&b, pr[pi].sig); pi++; }
                if (nb >= cb) { cb = cb ? cb*2 : 64; bins = realloc(bins, (size_t)cb*sizeof(bin_t)); }
                bins[nb++] = b;
            }
            if (fe >= te || (gj < gn && gstart[gj] > te)) break;
            cur = fe + 1;
            if (cur > te) break;
        }
    }

    /* leftRightMerge1: repeatedly absorb the smallest under-min bin into a
     * coordinate-adjacent neighbour (left preferred, matching R). */
    for (;;) {
        int32_t k, mn = -1;
        int32_t least = min_probes;
        for (k = 0; k < nb; k++) if (bins[k].n < least) { least = bins[k].n; mn = k; }
        if (mn < 0) break;
        int ml = (mn > 0 && bins[mn].start - 1 == bins[mn-1].end);
        int mr = (mn < nb-1 && bins[mn].end + 1 == bins[mn+1].start);
        if (ml) {
            bins[mn-1].end = bins[mn].end;
            for (k = 0; k < bins[mn].n; k++) bin_push(&bins[mn-1], bins[mn].sig[k]);
        } else if (mr) {
            bins[mn+1].start = bins[mn].start;
            for (k = 0; k < bins[mn].n; k++) bin_push(&bins[mn+1], bins[mn].sig[k]);
        }
        free(bins[mn].sig);
        memmove(&bins[mn], &bins[mn+1], (size_t)(nb-mn-1)*sizeof(bin_t));
        nb--;
    }

    for (int32_t k = 0; k < nb; k++) {
        if (*nbin >= *capbin) {
            *capbin = *capbin ? *capbin*2 : 1024;
            *bc = realloc(*bc, (size_t)*capbin*sizeof(char*));
            *bs = realloc(*bs, (size_t)*capbin*sizeof(int32_t));
            *be = realloc(*be, (size_t)*capbin*sizeof(int32_t));
            *bnp = realloc(*bnp, (size_t)*capbin*sizeof(int32_t));
            *bsig = realloc(*bsig, (size_t)*capbin*sizeof(double));
        }
        sesame__sort(bins[k].sig, bins[k].n);
        (*bc)[*nbin] = strdup(chrom);
        (*bs)[*nbin] = (int32_t)bins[k].start;
        (*be)[*nbin] = (int32_t)bins[k].end;
        (*bnp)[*nbin] = bins[k].n;
        (*bsig)[*nbin] = bins[k].n ? sesame__median_sorted(bins[k].sig, bins[k].n) : NAN;
        (*nbin)++;
        free(bins[k].sig);
    }
    free(bins);
    (void)gstart;
}

/* ------------------------------------------------------------- per sample - */

/* Compute one target sample's CNV into `r`. target/normals are this sample's and
 * the panel's totals (probe-major within a sample: target[i], normals[j*np+i]).
 * chrom/pos are per-probe coords; seq/gap drive binning. */
static int cnv_one(const char *sample, const double *target, const double *normals,
                   int32_t np, int32_t nn, const sesame_index_t *ix,
                   char *const *chrom, const int32_t *pos,
                   const seqinfo_t *seq, const gaps_t *gap,
                   int32_t tilewidth, int32_t min_probes,
                   sesame_cnv_t *r, sesame_err_t *err)
{
    int32_t *pb = NULL, npb = 0, i, j, k = nn + 1;
    double *A = NULL, *y = NULL, *beta = NULL, *R = NULL, *v = NULL, rss;
    int rc = SESAME_ERR_NOMEM;

    memset(r, 0, sizeof *r);
    r->sample = strdup(sample ? sample : "");
    r->n_normal = nn;

    /* pb = probes with target, all normals, and a coordinate present */
    pb = (int32_t *)malloc((size_t)np * sizeof(int32_t));
    if (!pb) { sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto done; }
    for (i = 0; i < np; i++) {
        int ok = !isnan(target[i]) && chrom[i][0] && pos[i] >= 0;
        for (j = 0; ok && j < nn; j++) if (isnan(normals[(size_t)j*(size_t)np + (size_t)i])) ok = 0;
        if (ok) pb[npb++] = i;
    }
    if (npb <= k) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT,
            "sample %s: only %d usable probes for %d normals", sample, npb, nn);
        goto done;
    }

    /* OLS design: [1, normal_1..normal_nn], response = target */
    A = (double *)malloc((size_t)npb * (size_t)k * sizeof(double));
    y = (double *)malloc((size_t)npb * sizeof(double));
    beta = (double *)malloc((size_t)k * sizeof(double));
    R = (double *)malloc((size_t)k * (size_t)k * sizeof(double));
    v = (double *)malloc((size_t)npb * sizeof(double));
    if (!A || !y || !beta || !R || !v) { sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto done; }
    for (i = 0; i < npb; i++) {
        int32_t p = pb[i];
        A[(size_t)i*(size_t)k + 0] = 1.0;
        for (j = 0; j < nn; j++) A[(size_t)i*(size_t)k + 1 + (size_t)j] = normals[(size_t)j*(size_t)np + (size_t)p];
        y[i] = target[p];
    }
    if (sesame__ols(A, npb, k, y, beta, &rss, R, v) != 0) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT, "sample %s: normal panel is rank-deficient", sample);
        goto done;
    }

    /* per-probe log2 ratio = log2( target / max(fitted, 1) ) */
    r->n_probe = npb; r->n_used = npb;
    r->probe_id = (char **)malloc((size_t)npb*sizeof(char*));
    r->chrom = (char **)malloc((size_t)npb*sizeof(char*));
    r->pos = (int32_t *)malloc((size_t)npb*sizeof(int32_t));
    r->log2ratio = (double *)malloc((size_t)npb*sizeof(double));
    if (!r->probe_id || !r->chrom || !r->pos || !r->log2ratio) { sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto done; }
    for (i = 0; i < npb; i++) {
        int32_t p = pb[i];
        double fitted = beta[0];
        for (j = 0; j < nn; j++) fitted += beta[1+j] * normals[(size_t)j*(size_t)np + (size_t)p];
        if (fitted < 1.0) fitted = 1.0;
        r->probe_id[i] = strdup(sesame_index_probe_id(ix, p));
        r->chrom[i] = strdup(chrom[p]);
        r->pos[i] = pos[p];
        r->log2ratio[i] = log2(target[p] / fitted);
    }

    /* binning, per chromosome in seqinfo order */
    {
        int32_t capbin = 0;
        cprobe_t *cp = (cprobe_t *)malloc((size_t)npb * sizeof(cprobe_t));
        long *gs = (long *)malloc((size_t)(gap->n ? gap->n : 1) * sizeof(long));
        long *ge = (long *)malloc((size_t)(gap->n ? gap->n : 1) * sizeof(long));
        if (!cp || !gs || !ge) { free(cp); free(gs); free(ge); sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto done; }
        for (int32_t c = 0; c < seq->n; c++) {
            int32_t m = 0, ng = 0;
            for (i = 0; i < npb; i++)
                if (!strcmp(r->chrom[i], seq->chrom[c])) { cp[m].pos = r->pos[i]; cp[m].sig = r->log2ratio[i]; m++; }
            if (m == 0) continue;
            for (i = 0; i < gap->n; i++)
                if (!strcmp(gap->chrom[i], seq->chrom[c])) { gs[ng] = gap->start[i]; ge[ng] = gap->end[i]; ng++; }
            bin_chrom(seq->chrom[c], seq->len[c], cp, m, gs, ge, ng, tilewidth, min_probes,
                      &r->bin_chrom, &r->bin_start, &r->bin_end, &r->bin_nprobe, &r->bin_log2ratio,
                      &r->n_bin, &capbin);
        }
        free(cp); free(gs); free(ge);
    }

    /* CBS segmentation over the bins, per chromosome. trimmed.SD (the sdundo
     * scale) is genome-wide over the finite bin signals, as in DNAcopy. */
    if (r->n_bin > 0) {
        int32_t nb = r->n_bin, na = 0, cs = 0, capseg = 0, bk;
        double *allsig = (double *)malloc((size_t)nb * sizeof(double));
        double tsd;
        if (!allsig) { sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto done; }
        for (bk = 0; bk < nb; bk++) if (!isnan(r->bin_log2ratio[bk])) allsig[na++] = r->bin_log2ratio[bk];
        tsd = sesame__trimmed_sd(allsig, na, 0.025);
        free(allsig);
        while (cs < nb) {
            int32_t ce = cs, cm, nf = 0, *bidx, *se, nseg, ns, lo;
            double *sv;
            while (ce < nb && !strcmp(r->bin_chrom[ce], r->bin_chrom[cs])) ce++;
            cm = ce - cs;
            sv = (double *)malloc((size_t)cm * sizeof(double));
            bidx = (int32_t *)malloc((size_t)cm * sizeof(int32_t));
            se = (int32_t *)malloc((size_t)cm * sizeof(int32_t));
            if (!sv || !bidx || !se) { free(sv); free(bidx); free(se);
                sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto done; }
            for (bk = cs; bk < ce; bk++)
                if (!isnan(r->bin_log2ratio[bk])) { sv[nf] = r->bin_log2ratio[bk]; bidx[nf] = bk; nf++; }
            nseg = nf > 0 ? sesame__cbs(sv, nf, 5, 0.001, 25, 200, 100, 1e-6, tsd, 2.2, se, nf) : 0;
            lo = 0;
            for (ns = 0; ns < nseg; ns++) {
                int32_t hi = se[ns], first = bidx[lo], last = bidx[hi-1];
                double sm = 0.0;
                for (bk = lo; bk < hi; bk++) sm += sv[bk];
                if (r->n_seg >= capseg) {
                    capseg = capseg ? capseg * 2 : 256;
                    r->seg_chrom = realloc(r->seg_chrom, (size_t)capseg * sizeof(char*));
                    r->seg_start = realloc(r->seg_start, (size_t)capseg * sizeof(int32_t));
                    r->seg_end = realloc(r->seg_end, (size_t)capseg * sizeof(int32_t));
                    r->seg_nbin = realloc(r->seg_nbin, (size_t)capseg * sizeof(int32_t));
                    r->seg_mean = realloc(r->seg_mean, (size_t)capseg * sizeof(double));
                }
                r->seg_chrom[r->n_seg] = strdup(r->bin_chrom[cs]);
                r->seg_start[r->n_seg] = r->bin_start[first];
                r->seg_end[r->n_seg] = r->bin_end[last];
                r->seg_nbin[r->n_seg] = hi - lo;
                r->seg_mean[r->n_seg] = sm / (hi - lo);
                r->n_seg++;
                lo = hi;
            }
            free(sv); free(bidx); free(se);
            cs = ce;
        }
    }
    rc = SESAME_OK;

done:
    free(pb); free(A); free(y); free(beta); free(R); free(v);
    return rc;
}

/* --------------------------------------------------------------- driver --- */

int sesame_cnv_run(const char *target_cg, const char *normals_cg,
                   const sesame_index_t *ix, const char *coords_path,
                   const char *seqinfo_path, const char *gaps_path,
                   int32_t tilewidth, int32_t min_probes,
                   sesame_cnv_t **out, int32_t *n_out, sesame_err_t *err)
{
    double *tmat = NULL, *nmat = NULL;
    char **tnames = NULL, **nnames = NULL, **chrom = NULL;
    int32_t tnp = 0, tns = 0, nnp = 0, nns = 0, np, *pos = NULL, s;
    seqinfo_t seq; gaps_t gap;
    sesame_cnv_t *res = NULL;
    int rc = SESAME_ERR_IO;

    memset(&seq, 0, sizeof seq); memset(&gap, 0, sizeof gap);
    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    np = sesame_index_nprobes(ix);

    if (sesame_read_cg_total(target_cg, &tmat, &tnp, &tns, &tnames, err) != SESAME_OK) goto done;
    if (sesame_read_cg_total(normals_cg, &nmat, &nnp, &nns, &nnames, err) != SESAME_OK) goto done;
    if (tnp != np || nnp != np) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT,
            "probe count mismatch: ordering %d, target %d, normals %d", np, tnp, nnp);
        goto done;
    }
    if (nns < 1) { rc = sesame__fail(err, SESAME_ERR_FORMAT, "no normals in %s", normals_cg); goto done; }
    if (load_coords(coords_path, np, &chrom, &pos, err) != SESAME_OK) goto done;
    if (load_seqinfo(seqinfo_path, &seq, err) != SESAME_OK) goto done;
    if (load_gaps(gaps_path, &gap, err) != SESAME_OK) goto done;

    res = (sesame_cnv_t *)calloc((size_t)tns, sizeof(sesame_cnv_t));
    if (!res) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto done; }
    for (s = 0; s < tns; s++) {
        /* target sample s (sample-major: tmat[s*np + i]); normals nmat[j*np + i] */
        if (cnv_one(tnames[s], tmat + (size_t)s*(size_t)np, nmat, np, nns, ix,
                    chrom, pos, &seq, &gap, tilewidth, min_probes, &res[s], err) != SESAME_OK) {
            sesame_cnv_free_array(res, s + 1);
            res = NULL; rc = err ? err->code : SESAME_ERR_IO; goto done;
        }
    }
    *out = res; *n_out = tns; res = NULL;
    rc = SESAME_OK;

done:
    free(tmat); free(nmat);
    if (tnames) { for (s = 0; s < tns; s++) free(tnames[s]); free(tnames); }
    if (nnames) { for (s = 0; s < nns; s++) free(nnames[s]); free(nnames); }
    if (chrom) { for (int32_t i = 0; i < np; i++) free(chrom[i]); free(chrom); }
    free(pos);
    free_seqinfo(&seq); free_gaps(&gap);
    if (res) sesame_cnv_free_array(res, *n_out);
    return rc;
}

void sesame_cnv_free_array(sesame_cnv_t *a, int32_t n)
{
    if (!a) return;
    for (int32_t s = 0; s < n; s++) {
        sesame_cnv_t *r = &a[s];
        int32_t i;
        free(r->sample);
        for (i = 0; i < r->n_probe; i++) { free(r->probe_id[i]); free(r->chrom[i]); }
        free(r->probe_id); free(r->chrom); free(r->pos); free(r->log2ratio);
        for (i = 0; i < r->n_bin; i++) free(r->bin_chrom[i]);
        free(r->bin_chrom); free(r->bin_start); free(r->bin_end);
        free(r->bin_nprobe); free(r->bin_log2ratio);
        for (i = 0; i < r->n_seg; i++) free(r->seg_chrom[i]);
        free(r->seg_chrom); free(r->seg_start); free(r->seg_end);
        free(r->seg_nbin); free(r->seg_mean);
    }
    free(a);
}
