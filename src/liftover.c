/* mLiftOver: move a .cg (or a .cm mask) from one row space to another.
 *
 * Every lift here is the same object -- an index-to-index map between a
 * SOURCE row space and a TARGET row space -- applied by one routine. What
 * differs between lifts is only how the map is built:
 *
 *   prefix join  array -> array   (R/mLiftOver.R, mapping=NULL/impute=FALSE)
 *   coordinates  array -> genome  (<plat>.<genome>.coord.tsv.gz vs the genome's
 *                                  cpg_nocontig.cr, YAME's CpG universe)
 *   inversion    genome -> array  (the coordinate map with its columns swapped)
 *
 * The map is a list of (src, tgt) PAIRS sorted by target, because that is the
 * one representation every case fits: one-to-one, many-to-one (array
 * replicates on one CpG), one-to-many (one CpG under several probes), and a
 * miss (a row in no pair). Applying it is a single pass over the pairs: each
 * run of equal targets is reduced onto that target -- numeric formats by the
 * mean of the non-NA sources (real M/U counts are pooled instead), bit masks
 * by OR, 2-bit codes by max, anything else first-wins -- and a target in no
 * pair keeps the format's own "absent" (NA, 0/0, bit 0, code 0). Output is
 * positional to the TARGET space, whatever it is.
 *
 * The coordinate map is built in memory each run from the two store files,
 * via YAME's own row index over the .cr (init_finder / row_finder_search:
 * coarse 2^17-bp bins per chromosome, then a short scan) -- ~2 s for 284k
 * probes against 29.4M CpGs, ~4 MB retained, nothing persisted. No rowmap
 * asset to publish or go stale. row_finder_search takes a 1-BASED position;
 * CpG_beg is 0-based, hence +1. Get that wrong and the map comes back nearly
 * empty rather than wrong, which is why the hit counts are always reported.
 *
 * A probe with no CpG row (rs, ch, unmapped) is in no pair: dropped and
 * counted. Genome -> array is lossy by construction: CpGs under no probe are
 * gone. The empirical liftOver.* quality mappings and impute=TRUE of the R
 * function are not ported.
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <zlib.h>
#include "cfile.h"    /* open_cfile, read_cdata1, cdata_write1, BGZF */
#include "cdata.h"    /* cdata_t, decompress, f3_get_mu/f3_set_mu, row_finder */
#include "index.h"    /* loadSampleNamesFromIndex, cleanSampleNames2 */
#include "summary.h"  /* prepare_mask: a format-0/1 block -> fmt0 bitset */

/* ------------------------------------------------------- prefix join --- */

static int is_modern(const char *p)
{
    return p && (!strcmp(p, "EPICv2") || !strcmp(p, "MSA"));
}
static int is_legacy(const char *p)
{
    return p && (!strcmp(p, "EPIC") || !strcmp(p, "HM450") || !strcmp(p, "HM27"));
}

/* Compare the prefixes of two IDs. A side with strip=1 ends its prefix at the
 * first '_' (or end); strip=0 uses the whole ID. Shorter-prefix-first ordering,
 * C locale -- matches R's strsplit("_")[[1]][1] keys. */
static int pcmp(const char *a, int as, const char *b, int bs)
{
    for (;;) {
        char ca = *a, cb = *b;
        int ta = (ca == '\0') || (as && ca == '_');
        int tb = (cb == '\0') || (bs && cb == '_');
        if (ta && tb) return 0;
        if (ta) return -1;
        if (tb) return 1;
        if (ca != cb) return (unsigned char)ca < (unsigned char)cb ? -1 : 1;
        a++; b++;
    }
}

typedef struct { const char *id; int32_t idx; } lo_ent;

/* Sort source entries by (prefix, idx) so the first matching entry for any
 * prefix carries the lowest source index -- R's distinct(.keep_all) first. */
static int src_strip_g;
static int lo_cmp(const void *x, const void *y)
{
    const lo_ent *a = (const lo_ent *)x, *b = (const lo_ent *)y;
    int c = pcmp(a->id, src_strip_g, b->id, src_strip_g);
    if (c) return c;
    return a->idx < b->idx ? -1 : a->idx > b->idx ? 1 : 0;
}

/* --------------------------------------------------------------- map --- */

typedef struct { int64_t src, tgt; } lo_pair;

static int lo_pair_cmp(const void *a, const void *b)
{
    const lo_pair *x = (const lo_pair *)a, *y = (const lo_pair *)b;
    return x->tgt < y->tgt ? -1 : x->tgt > y->tgt ? 1 :
           x->src < y->src ? -1 : x->src > y->src ? 1 : 0;
}
static int lo_cmp_i64(const void *a, const void *b)
{
    int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
    return x < y ? -1 : x > y ? 1 : 0;
}

/* Sort the pairs by (tgt, src) and fill the counts. Takes ownership of pairs. */
static int lo_finish(sesame_rowmap_t *m, lo_pair *pairs, int64_t npair, sesame_err_t *err)
{
    int64_t i, *tmp;
    qsort(pairs, (size_t)npair, sizeof *pairs, lo_pair_cmp);
    m->npair = npair;
    m->src = (int64_t *)malloc((size_t)(npair ? npair : 1) * sizeof(int64_t));
    m->tgt = (int64_t *)malloc((size_t)(npair ? npair : 1) * sizeof(int64_t));
    tmp    = (int64_t *)malloc((size_t)(npair ? npair : 1) * sizeof(int64_t));
    if (!m->src || !m->tgt || !tmp) {
        free(tmp); free(pairs); return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
    }
    for (i = 0; i < npair; i++) { m->src[i] = pairs[i].src; m->tgt[i] = pairs[i].tgt; }
    free(pairs);
    m->ntgt_covered = m->ntgt_multi = m->nsrc_mapped = 0;
    for (i = 0; i < npair; i++) {                 /* runs of equal target */
        if (i == 0 || m->tgt[i] != m->tgt[i-1]) m->ntgt_covered++;
        else if (i == 1 || m->tgt[i-1] != m->tgt[i-2]) m->ntgt_multi++;
    }
    memcpy(tmp, m->src, (size_t)npair * sizeof(int64_t));
    qsort(tmp, (size_t)npair, sizeof(int64_t), lo_cmp_i64);
    for (i = 0; i < npair; i++) if (i == 0 || tmp[i] != tmp[i-1]) m->nsrc_mapped++;
    free(tmp);
    return SESAME_OK;
}

void sesame_rowmap_free(sesame_rowmap_t *m)
{
    if (!m) return;
    free(m->src); free(m->tgt); free(m);
}

int sesame_rowmap_prefix(const char *src_platform, const sesame_index_t *src_ix,
                         const char *tgt_platform, const sesame_index_t *tgt_ix,
                         sesame_rowmap_t **out, sesame_err_t *err)
{
    int32_t nsrc, ntgt, i;
    int src_strip = 0, tgt_strip = 0;
    lo_ent *ent = NULL;
    lo_pair *pairs = NULL;
    int64_t npair = 0;
    sesame_rowmap_t *m = NULL;
    int rc;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!src_ix || !tgt_ix || !out) return sesame__fail(err, SESAME_ERR_IO, "null argument");
    nsrc = sesame_index_nprobes(src_ix);
    ntgt = sesame_index_nprobes(tgt_ix);
    if (nsrc <= 0 || ntgt <= 0) return sesame__fail(err, SESAME_ERR_FORMAT, "empty index");

    /* Direction: strip the suffix on the modern side of a modern<->legacy lift;
     * otherwise join on the full ID (same family, incl. identity). */
    if (is_legacy(tgt_platform) && is_modern(src_platform)) src_strip = 1;
    else if (is_modern(tgt_platform) && is_legacy(src_platform)) tgt_strip = 1;

    ent   = (lo_ent *)malloc((size_t)nsrc * sizeof *ent);
    pairs = (lo_pair *)malloc((size_t)ntgt * sizeof *pairs);   /* at most one per target */
    m     = (sesame_rowmap_t *)calloc(1, sizeof *m);
    if (!ent || !pairs || !m) { free(ent); free(pairs); free(m); return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }

    for (i = 0; i < nsrc; i++) { ent[i].id = sesame_index_probe_id(src_ix, i); ent[i].idx = i; }
    src_strip_g = src_strip;
    qsort(ent, (size_t)nsrc, sizeof *ent, lo_cmp);

    /* For each target probe, lower-bound its prefix in the sorted source keys. */
    for (i = 0; i < ntgt; i++) {
        const char *tid = sesame_index_probe_id(tgt_ix, i);
        int32_t lo = 0, hi = nsrc;             /* first ent with prefix >= tid */
        while (lo < hi) {
            int32_t mid = lo + ((hi - lo) >> 1);
            if (pcmp(ent[mid].id, src_strip, tid, tgt_strip) < 0) lo = mid + 1;
            else hi = mid;
        }
        if (lo < nsrc && pcmp(ent[lo].id, src_strip, tid, tgt_strip) == 0) {
            pairs[npair].src = ent[lo].idx; pairs[npair].tgt = i; npair++;
        }
    }
    free(ent);
    m->nsrc = nsrc; m->ntgt = ntgt;
    if ((rc = lo_finish(m, pairs, npair, err)) != SESAME_OK) { free(m); return rc; }
    *out = m;
    return SESAME_OK;
}

/* Load CpG_chrm + CpG_beg from the coord table, positional in the ordering.
 * chrm[i] strdup'd ("" if unmapped: "", "*" or "NA"), beg0[i] or -1. */
static int lo_load_coords(const char *path, int32_t np, char ***chrm_out,
                          long **beg_out, sesame_err_t *err)
{
    gzFile f = gzopen(path, "rb");
    char *line; size_t cap = 1 << 16;
    char **chrm; long *beg;
    int32_t row = 0;

    if (!f) return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", path);
    line = (char *)malloc(cap);
    chrm = (char **)calloc((size_t)np, sizeof(char *));
    beg  = (long *)malloc((size_t)np * sizeof(long));
    if (!line || !chrm || !beg) {
        free(line); free(chrm); free(beg); gzclose(f);
        return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
    }
    gzgets(f, line, (int)cap);                              /* header */
    while (gzgets(f, line, (int)cap)) {
        char *tab, *tab2;
        if (row >= np) { row++; continue; }                 /* count overflow */
        line[strcspn(line, "\r\n")] = '\0';
        tab = strchr(line, '\t'); if (tab) *tab = '\0';
        if (line[0] == '\0' || !strcmp(line, "*") || !strcmp(line, "NA")) {
            chrm[row] = strdup(""); beg[row] = -1;
        } else {
            chrm[row] = strdup(line);
            tab2 = tab ? strchr(tab + 1, '\t') : NULL;
            if (tab2) *tab2 = '\0';
            beg[row] = tab ? strtol(tab + 1, NULL, 10) : -1;
        }
        row++;
    }
    free(line); gzclose(f);
    if (row != np) {
        for (int32_t i = 0; i < row && i < np; i++) free(chrm[i]);
        free(chrm); free(beg);
        return sesame__fail(err, SESAME_ERR_FORMAT,
            "%s has %d data rows, ordering has %d -- lineage mismatch", path, row, np);
    }
    *chrm_out = chrm; *beg_out = beg;
    return SESAME_OK;
}

int sesame_rowmap_coords(const char *coords_path, int32_t nprobe,
                         const char *cr_path, sesame_rowmap_t **out,
                         sesame_err_t *err)
{
    sesame_rowmap_t *m = NULL;
    char **chrm = NULL; long *beg = NULL;
    char crbuf[4096];
    cfile_t cf;
    cdata_t cr;
    row_finder_t fdr;
    lo_pair *pairs = NULL;
    int64_t npair = 0;
    int32_t i, rc;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!coords_path || !cr_path || !out || nprobe <= 0)
        return sesame__fail(err, SESAME_ERR_IO, "null argument");
    if ((rc = lo_load_coords(coords_path, nprobe, &chrm, &beg, err)) != SESAME_OK)
        return rc;

    snprintf(crbuf, sizeof crbuf, "%s", cr_path);
    cf = open_cfile(crbuf);
    if (!cf.fh) { rc = sesame__fail(err, SESAME_ERR_IO, "cannot open %s", cr_path); goto out; }
    cr = read_cdata1(&cf);
    bgzf_close(cf.fh);
    if (cr.n == 0 || cr.fmt != '7') {
        free_cdata(&cr);
        rc = sesame__fail(err, SESAME_ERR_FORMAT,
            "%s is not a format-7 row-coordinate track", cr_path);
        goto out;
    }
    m     = (sesame_rowmap_t *)calloc(1, sizeof *m);
    pairs = (lo_pair *)malloc((size_t)nprobe * sizeof *pairs);
    if (!m || !pairs) { free_cdata(&cr); rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto out; }
    /* cr.n on a format-7 record is its payload length, not its row count;
     * fmt7_data_length() walks the rows. This is the universe size every
     * genome-indexed file must have -- 29,401,795 for hg38's cpg_nocontig. */
    m->nsrc = nprobe; m->ntgt = (int64_t)fmt7_data_length(&cr);

    fdr = init_finder(&cr);
    for (i = 0; i < nprobe; i++) {
        uint64_t r = 0;
        /* row_finder_search() exits the process on a chromosome the track does
         * not carry -- an alt/random-contig probe against cpg_nocontig.cr, which
         * drops contigs by design -- so check membership first. Such a probe
         * is simply unmapped, and counts as one. */
        if (chrm[i][0] && beg[i] >= 0 &&
            kh_get(str2int, fdr.h, chrm[i]) != kh_end(fdr.h))
            r = row_finder_search(chrm[i], (uint64_t)beg[i] + 1, &fdr, &cr);  /* beg1 */
        if (r) { pairs[npair].src = i; pairs[npair].tgt = (int64_t)r - 1; npair++; }
    }
    free_row_finder(&fdr);
    free_cdata(&cr);
    if ((rc = lo_finish(m, pairs, npair, err)) != SESAME_OK) { pairs = NULL; goto out; }
    pairs = NULL;
    *out = m; m = NULL;
out:
    for (i = 0; i < nprobe; i++) free(chrm[i]);
    free(chrm); free(beg); free(pairs); free(m);
    return rc;
}

int sesame_rowmap_invert(const sesame_rowmap_t *in, sesame_rowmap_t **out,
                         sesame_err_t *err)
{
    sesame_rowmap_t *m;
    lo_pair *pairs;
    int64_t i;
    int rc;
    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!in || !out) return sesame__fail(err, SESAME_ERR_IO, "null argument");
    m = (sesame_rowmap_t *)calloc(1, sizeof *m);
    pairs = (lo_pair *)malloc((size_t)(in->npair ? in->npair : 1) * sizeof *pairs);
    if (!m || !pairs) { free(m); free(pairs); return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }
    m->nsrc = in->ntgt; m->ntgt = in->nsrc;
    for (i = 0; i < in->npair; i++) { pairs[i].src = in->tgt[i]; pairs[i].tgt = in->src[i]; }
    if ((rc = lo_finish(m, pairs, in->npair, err)) != SESAME_OK) { free(m); return rc; }
    *out = m;
    return SESAME_OK;
}

/* ------------------------------------------------------------- apply --- */

/* Write the companion .idx (same shape cgwrite.c writes). */
static void lo_write_idx(const char *path, char *const *names,
                         const int64_t *offs, int32_t nsamp)
{
    size_t z = strlen(path) + 5;
    char *idxpath = (char *)malloc(z);
    FILE *ix;
    if (!idxpath) return;
    snprintf(idxpath, z, "%s.idx", path);
    if ((ix = fopen(idxpath, "w"))) {
        for (int32_t j = 0; j < nsamp; j++)
            fprintf(ix, "%s\t%lld\n", names[j], (long long)offs[j]);
        fclose(ix);
    }
    free(idxpath);
}

/* The read block and its inflated form: two objects normally, one when a mask
 * was inflated in place. */
static void lo_free2(cdata_t *c, cdata_t *d, int same)
{
    free_cdata(c);
    if (!same) free_cdata(d);
}

/* A source row's beta, for the numeric reductions: format 4 as stored (NA if
 * negative), format 3 as M/(M+U) (NA if 0,0). */
static double lo_beta(const cdata_t *d, int64_t s)
{
    if (d->fmt == '4') { float v = ((const float *)d->s)[s]; return v < 0.0f ? NAN : (double)v; }
    else {
        uint64_t mu = f3_get_mu((cdata_t *)d, (uint64_t)s);
        uint64_t M = mu >> 32, U = mu & 0xffffffffULL;
        return (M + U) ? (double)M / (double)(M + U) : NAN;
    }
}

int sesame_liftover_apply(const char *in_cx, const char *out_cx,
                          const sesame_rowmap_t *m, int depth, sesame_err_t *err)
{
    cfile_t cf;
    snames_t sn;
    BGZF *fp = NULL;
    char inbuf[4096];
    int64_t *offs = NULL, k;
    char **names = NULL;
    int32_t ns = 0, cap = 0, i, rc = SESAME_OK;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!in_cx || !out_cx || !m) return sesame__fail(err, SESAME_ERR_IO, "null argument");
    if (depth < 0) return sesame__fail(err, SESAME_ERR_IO, "depth must be >= 0");

    snprintf(inbuf, sizeof inbuf, "%s", in_cx);
    cf = open_cfile(inbuf);
    if (!cf.fh) return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", in_cx);
    sn = loadSampleNamesFromIndex(inbuf);
    if (!(fp = bgzf_open2(out_cx, "w"))) {
        rc = sesame__fail(err, SESAME_ERR_IO, "cannot open %s for writing", out_cx); goto out;
    }

    for (;;) {
        cdata_t c = read_cdata1(&cf), d, o;
        int in3, in4, numeric, is_mask;
        if (c.n == 0) break;                                   /* EOF */
        /* A mask block (format 0/1) is inflated the way mask.c does it --
         * prepare_mask -> convertToFmt0, in place, c.n then the bit count --
         * because decompress() hands back 0 rows for a compressed format-0
         * record. Everything else inflates with decompress(). */
        is_mask = c.fmt < '2';
        if (is_mask) { prepare_mask(&c); d = c; }
        else d = decompress(c);
        if ((int64_t)d.n != m->nsrc) {
            lo_free2(&c, &d, is_mask);
            rc = sesame__fail(err, SESAME_ERR_FORMAT,
                "%s has %llu rows but the source row space has %lld -- not indexed to it",
                in_cx, (unsigned long long)d.n, (long long)m->nsrc);
            goto out;
        }
        if (d.fmt == '7') {
            lo_free2(&c, &d, is_mask);
            rc = sesame__fail(err, SESAME_ERR_UNSUPPORTED,
                "%s is a row-coordinate track (format 7), which is the map, not data", in_cx);
            goto out;
        }
        in3 = d.fmt == '3'; in4 = d.fmt == '4'; numeric = in3 || in4;
        if (depth > 0 && !numeric) {
            lo_free2(&c, &d, is_mask);
            rc = sesame__fail(err, SESAME_ERR_UNSUPPORTED,
                "--simulated-depth needs betas (format 4) or M/U (format 3); %s is format %c",
                in_cx, d.fmt);
            goto out;
        }

        /* The output record, initialised to the format's own "absent" so a
         * target in no pair reads as missing, never as a value. */
        memset(&o, 0, sizeof o);
        o.compressed = 0; o.n = (uint64_t)m->ntgt;
        if (depth > 0 || in3) {                                /* M/U, 0,0 = missing */
            o.fmt = '3'; o.unit = 8; o.s = (uint8_t *)calloc((size_t)m->ntgt, 8);
        } else if (in4) {                                      /* betas, NA = -1.0 */
            float *s = (float *)malloc((size_t)m->ntgt * sizeof(float));
            o.fmt = '4'; o.unit = sizeof(float); o.s = (uint8_t *)s;
            if (s) for (int64_t r = 0; r < m->ntgt; r++) s[r] = -1.0f;
        } else if (d.fmt == '0') {                             /* 1 bit per row */
            o.fmt = '0'; o.unit = d.unit; o.s = (uint8_t *)calloc(((size_t)m->ntgt + 7) >> 3, 1);
        } else if (d.fmt == '6') {                             /* 2 bits per row */
            o.fmt = '6'; o.unit = d.unit; o.s = (uint8_t *)calloc(((size_t)m->ntgt + 3) >> 2, 1);
        } else if (d.fmt == '2') {                             /* state track */
            /* Layout is [keys...][\0][rows...]: the key table rides across
             * unchanged (as rowsub's slice does it) and the rows follow it.
             * calloc leaves an unmapped row at code 0 -- the track's FIRST
             * state, which for an MRMP is the null pattern (Pna). */
            uint64_t keys_nb = fmt2_get_keys_nbytes(&d);
            o.fmt = '2'; o.unit = d.unit;
            o.s = (uint8_t *)calloc(1, (size_t)keys_nb + 1 + (size_t)m->ntgt * d.unit);
            if (o.s) memcpy(o.s, d.s, (size_t)keys_nb + 1);
        } else {                                               /* '1', '5': unit-wide */
            o.fmt = d.fmt; o.unit = d.unit; o.s = (uint8_t *)calloc((size_t)m->ntgt, d.unit);
        }
        if (!o.s) { lo_free2(&c, &d, is_mask); rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto out; }
        /* where the rows start, for the unit-wide copy: after the key table
         * on a state track, at the buffer start otherwise */
        const uint8_t *ibase = d.fmt == '2' ? fmt2_get_data(&d) : d.s;
        uint8_t *obase = d.fmt == '2' ? o.s + fmt2_get_keys_nbytes(&d) + 1 : o.s;

        /* One pass over the pairs; each run of equal targets is one reduction. */
        for (k = 0; k < m->npair; ) {
            int64_t t = m->tgt[k], k2 = k;
            if (numeric && (depth > 0 || in4)) {               /* mean of non-NA betas */
                double sum = 0.0; int n = 0;
                for (; k2 < m->npair && m->tgt[k2] == t; k2++) {
                    double v = lo_beta(&d, m->src[k2]);
                    if (!isnan(v)) { sum += v; n++; }
                }
                if (n) {
                    double v = sum / n;
                    if (depth > 0) {
                        uint64_t M = (uint64_t)llround(v * depth);
                        f3_set_mu(&o, (uint64_t)t, M, (uint64_t)depth - M);
                    } else ((float *)o.s)[t] = (float)v;
                }
            } else if (in3) {                                  /* real counts: pool */
                uint64_t M = 0, U = 0;
                for (; k2 < m->npair && m->tgt[k2] == t; k2++) {
                    uint64_t mu = f3_get_mu(&d, (uint64_t)m->src[k2]);
                    M += mu >> 32; U += mu & 0xffffffffULL;
                }
                f3_set_mu(&o, (uint64_t)t, M, U);
            } else if (d.fmt == '0') {                         /* OR */
                int bit = 0;
                for (; k2 < m->npair && m->tgt[k2] == t; k2++) {
                    int64_t s = m->src[k2];
                    if (d.s[s >> 3] & (1u << (s & 7))) bit = 1;
                }
                if (bit) o.s[t >> 3] |= (uint8_t)(1u << (t & 7));
            } else if (d.fmt == '6') {                         /* max code */
                uint8_t v = 0;
                for (; k2 < m->npair && m->tgt[k2] == t; k2++) {
                    int64_t s = m->src[k2];
                    uint8_t x = (d.s[s >> 2] >> ((s & 3) * 2)) & 3;
                    if (x > v) v = x;
                }
                o.s[t >> 2] |= (uint8_t)(v << ((t & 3) * 2));
            } else {                                           /* first wins */
                memcpy(obase + (size_t)o.unit * (size_t)t, ibase + (size_t)d.unit * (size_t)m->src[k], o.unit);
                for (; k2 < m->npair && m->tgt[k2] == t; k2++) ;
            }
            k = k2;
        }

        if (ns >= cap) {
            cap = cap ? cap * 2 : 8;
            offs = (int64_t *)realloc(offs, (size_t)cap * sizeof(int64_t));
            names = (char **)realloc(names, (size_t)cap * sizeof(char *));
        }
        offs[ns] = bgzf_tell(fp);
        cdata_compress(&o);
        cdata_write1(fp, &o);
        free(o.s);
        names[ns] = strdup(ns < sn.n ? sn.s[ns] : "");
        ns++;
        lo_free2(&c, &d, is_mask);
    }
    bgzf_close(fp); fp = NULL;
    lo_write_idx(out_cx, names, offs, ns);
out:
    if (fp) bgzf_close(fp);
    bgzf_close(cf.fh);
    cleanSampleNames2(sn);
    for (i = 0; i < ns; i++) free(names[i]);
    free(names); free(offs);
    return rc;
}

/* --------------------------------------------- matrix form, for the API --- */

int sesame_liftover_betas(const char *src_platform, const sesame_index_t *src_ix,
                          const char *tgt_platform, const sesame_index_t *tgt_ix,
                          const double *mat_in, int32_t nsamp,
                          double **mat_out, sesame_err_t *err)
{
    sesame_rowmap_t *m = NULL;
    double *out;
    int64_t k;
    int32_t j, i, rc;

    if (!mat_in || !mat_out) return sesame__fail(err, SESAME_ERR_IO, "null argument");
    if ((rc = sesame_rowmap_prefix(src_platform, src_ix, tgt_platform, tgt_ix, &m, err)) != SESAME_OK)
        return rc;
    out = (double *)malloc((size_t)nsamp * (size_t)m->ntgt * sizeof(double));
    if (!out) { sesame_rowmap_free(m); return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }
    for (j = 0; j < nsamp; j++) {
        const double *si = mat_in + (size_t)j * (size_t)m->nsrc;
        double *ti = out + (size_t)j * (size_t)m->ntgt;
        for (i = 0; i < (int32_t)m->ntgt; i++) ti[i] = NAN;
        for (k = 0; k < m->npair; k++) ti[m->tgt[k]] = si[m->src[k]];   /* one per target */
    }
    sesame_rowmap_free(m);
    *mat_out = out;
    return SESAME_OK;
}
