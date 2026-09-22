/* vcf.c -- genotype the SNP (rs) probes on an Infinium array and write a VCF,
 * a port of sesame's formatVCF (R/vcf.R).
 *
 * Genotyping uses the RAW signal (no channel inference): for a SNP that flips an
 * Infinium-I probe's color channel, the out-of-band signal IS the alt-allele
 * evidence, so inferInfiniumIChannel would erase exactly what we measure. Per
 * SNP probe we form a variant-allele fraction (vaf):
 *   U == "ALT"       -> vaf = 1 - beta          (the U bead reports ALT)
 *   U == "REF_InfI"  -> vaf = out-of-band AF     (Type-I channel-switch SNP)
 *   otherwise (REF)  -> vaf = beta
 * then a 3-genotype binomial model (background 0.1, 40 beads) picks 0/0, 0/1, or
 * 1/1 and a phred-like score. Sites-only VCF, INFO carries PVF/GT/GS/Probe_ID/rs.
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

/* --- binomial genotyper (sesame R/vcf.R genotyper) --------------------------
 * GL_j = dbinom(round(x*40), 40, p_j) for p in {0.1, 0.5, 0.9}. The binomial
 * coefficient is the same across j (same k), so it cancels in the argmax and in
 * the posterior; we work in log space to avoid underflow. */
static void genotyper(double x, char gt[4], int *gs)
{
    static const double p[3] = { 0.1, 0.5, 0.9 };
    static const char  *g[3] = { "0/0", "0/1", "1/1" };
    const double n = 40.0;
    double k, lg[3], rest = 0.0, oneminus, score;
    int i, ind = 0;

    if (isnan(x)) { gt[0] = '.'; gt[1] = '/'; gt[2] = '.'; gt[3] = 0; *gs = 0; return; }
    k = round(x * n);
    if (k < 0) k = 0; else if (k > n) k = n;
    for (i = 0; i < 3; i++) lg[i] = k * log(p[i]) + (n - k) * log(1.0 - p[i]);
    for (i = 1; i < 3; i++) if (lg[i] > lg[ind]) ind = i;
    for (i = 0; i < 3; i++) if (i != ind) rest += exp(lg[i] - lg[ind]);
    oneminus = rest / (1.0 + rest);          /* = 1 - GL[ind]/sum(GL), precisely */
    score = floor(-log10(oneminus) * 10.0);
    if (!isfinite(score) || score > 255) score = 255;
    memcpy(gt, g[ind], 4);
    *gs = (int)score;
}

/* Type-I out-of-band allele fraction (getAFTypeIbySumAlleles): the share of a
 * probe's total signal that falls in the OTHER color channel. */
static double af_typeI(const sesame_sigdf_t *sdf, int32_t i)
{
    double mg = sdf->MG[i], ug = sdf->UG[i], mr = sdf->MR[i], ur = sdf->UR[i];
    if (isnan(mg) || isnan(ug) || isnan(mr) || isnan(ur)) return NAN;
    if (sdf->col[i] == SESAME_COL_G)         /* green in-band: OOB = red */
        return fmax(mr + ur, 1.0) / fmax(mr + ur + mg + ug, 2.0);
    if (sdf->col[i] == SESAME_COL_R)         /* red in-band: OOB = green */
        return fmax(mg + ug, 1.0) / fmax(mr + ur + mg + ug, 2.0);
    return NAN;                              /* Inf-II: not a Type-I SNP */
}

/* --- Probe_ID -> row index hash (open addressing, linear probe) ------------- */
typedef struct { char **key; int32_t *val; int32_t cap; } idmap_t;

static uint32_t hstr(const char *s) { uint32_t h = 2166136261u; while (*s) { h ^= (uint8_t)*s++; h *= 16777619u; } return h; }

static int idmap_build(idmap_t *m, const sesame_index_t *ix)
{
    int32_t n = sesame_index_nprobes(ix), i, cap = 1;
    while (cap < n * 2) cap <<= 1;
    m->cap = cap;
    m->key = (char **)calloc((size_t)cap, sizeof(char *));
    m->val = (int32_t *)malloc((size_t)cap * sizeof(int32_t));
    if (!m->key || !m->val) return -1;
    for (i = 0; i < n; i++) {
        const char *id = sesame_index_probe_id(ix, i);
        uint32_t h = hstr(id) & (uint32_t)(cap - 1);
        while (m->key[h]) h = (h + 1) & (uint32_t)(cap - 1);
        m->key[h] = (char *)id; m->val[h] = i;   /* id points into the index blob */
    }
    return 0;
}
static int32_t idmap_get(const idmap_t *m, const char *id)
{
    uint32_t h = hstr(id) & (uint32_t)(m->cap - 1);
    while (m->key[h]) { if (!strcmp(m->key[h], id)) return m->val[h]; h = (h + 1) & (uint32_t)(m->cap - 1); }
    return -1;
}
static void idmap_free(idmap_t *m) { free(m->key); free(m->val); }

/* --- VCF records ------------------------------------------------------------ */
typedef struct {
    char chrom[24]; long pos; char ref[8], alt[8], rs[80], probe[48];
    double vaf; char gt[4]; int gs;
} rec_t;

static int rec_cmp(const void *a, const void *b)          /* chrom lexicographic, then pos */
{
    const rec_t *x = (const rec_t *)a, *y = (const rec_t *)b;
    int c = strcmp(x->chrom, y->chrom);
    if (c) return c;
    return x->pos < y->pos ? -1 : (x->pos > y->pos ? 1 : 0);
}

/* ACT -> H, AGT -> D (IUPAC-ish encoding used by sesame for tri-allelic REF/ALT). */
static void enc_allele(char *a) { if (!strcmp(a, "ACT")) strcpy(a, "H"); else if (!strcmp(a, "AGT")) strcpy(a, "D"); }

/* one logical (possibly long) gz line into buf, newline stripped. 1/0/-1. */
static int gzline(gzFile f, char **buf, size_t *cap)
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

int sesame_format_vcf(const sesame_sigdf_t *sdf, const char *snp_path,
                      const char *genome, int variants_only, FILE *out,
                      sesame_err_t *err)
{
    idmap_t map; gzFile f;
    char *line = NULL; size_t cap = 1 << 16;
    double *betas = NULL;
    rec_t *rec = NULL; int32_t nrec = 0, crec = 0, i;
    int rc = SESAME_ERR_IO;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    memset(&map, 0, sizeof map);

    betas = (double *)malloc((size_t)sdf->n * sizeof(double));
    if (!betas) return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
    if (sesame_get_betas(sdf, 0, betas, err) != SESAME_OK) { free(betas); return err ? err->code : SESAME_ERR_IO; }
    if (idmap_build(&map, sdf->ix) != 0) { free(betas); sesame__fail(err, SESAME_ERR_NOMEM, "oom"); return SESAME_ERR_NOMEM; }

    if (!(f = gzopen(snp_path, "rb"))) { rc = sesame__fail(err, SESAME_ERR_IO, "cannot open %s", snp_path); goto done; }
    if (!(line = (char *)malloc(cap))) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto done; }
    gzline(f, &line, &cap);                                /* header: chrm beg end strand rs designType U REF ALT Probe_ID */

    while (gzline(f, &line, &cap) == 1) {
        char *tok[10], *s = line; int nf = 0;
        double beta, vaf; int32_t idx;
        for (;;) { char *t = strchr(s, '\t'); if (nf < 10) tok[nf] = s; nf++; if (!t) break; *t = 0; s = t + 1; }
        if (nf < 10) continue;
        /* tok: 0 chrm, 1 beg, 2 end, 3 strand, 4 rs, 5 designType, 6 U, 7 REF, 8 ALT, 9 Probe_ID */
        /* variants_only: drop the non-switching Infinium-I probes (REF_InfI with
         * no overlapping rs), which carry no known variant and almost never call
         * one; keep rs probes, channel-switching Inf-I, and Type-II SNP probes. */
        if (variants_only && !strcmp(tok[6], "REF_InfI") && !strcmp(tok[4], "NA")) continue;
        idx = idmap_get(&map, tok[9]);
        if (idx < 0) continue;                             /* SNP probe not on this array */
        beta = betas[idx];
        if (!strcmp(tok[6], "ALT"))          vaf = 1.0 - beta;
        else if (!strcmp(tok[6], "REF_InfI")) vaf = af_typeI(sdf, idx);
        else                                  vaf = beta;
        if (isnan(vaf)) continue;                          /* R yields NA; drop the site */

        if (nrec >= crec) { crec = crec ? crec * 2 : 4096; rec = (rec_t *)realloc(rec, (size_t)crec * sizeof(rec_t)); }
        {
            rec_t *r = &rec[nrec++];
            snprintf(r->chrom, sizeof r->chrom, "%s", tok[0]);
            r->pos = strtol(tok[2], NULL, 10);
            snprintf(r->ref, sizeof r->ref, "%s", tok[7]); enc_allele(r->ref);
            snprintf(r->alt, sizeof r->alt, "%s", tok[8]); enc_allele(r->alt);
            snprintf(r->rs, sizeof r->rs, "%s", tok[4]);
            snprintf(r->probe, sizeof r->probe, "%s", tok[9]);
            r->vaf = vaf;
            genotyper(vaf, r->gt, &r->gs);
        }
    }
    gzclose(f); f = NULL;

    qsort(rec, (size_t)nrec, sizeof(rec_t), rec_cmp);

    fprintf(out, "##fileformat=VCFv4.0\n##source=sesame-cli %s\n##reference=%s\n",
            SESAME_VERSION, genome ? genome : "hg38");
    fputs("##INFO=<ID=PVF,Number=1,Type=Float,Description=\"Pseudo Variant Frequency\">\n"
          "##INFO=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
          "##INFO=<ID=GS,Number=1,Type=Integer,Description=\"Genotyping score from 7 to 85\">\n"
          "##INFO=<ID=Probe_ID,Number=1,Type=String,Description=\"Infinium Probe ID\">\n"
          "##INFO=<ID=rs_ID,Number=1,Type=String,Description=\"Overlapping rs ID from dbSNP\">\n"
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", out);
    for (i = 0; i < nrec; i++) {
        rec_t *r = &rec[i];
        fprintf(out, "%s\t%ld\t.\t%s\t%s\t%d\t%s\tPVF=%.3f;GT=%s;GS=%d;Probe_ID=%s",
                r->chrom, r->pos, r->ref, r->alt, r->gs, r->gs > 20 ? "PASS" : "FAIL",
                r->vaf, r->gt, r->gs, r->probe);
        if (strcmp(r->rs, "NA") != 0 && r->rs[0]) fprintf(out, ";rs_ID=%s", r->rs);
        fputc('\n', out);
    }
    rc = SESAME_OK;

done:
    if (f) gzclose(f);
    free(line); free(betas); free(rec); idmap_free(&map);
    return rc;
}
