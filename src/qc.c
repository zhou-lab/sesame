/* qc.c -- the sesameQC panel (R sesameQC_calcStats), computed from a raw SigDF.
 *
 * Groups mirror R/QC.R: detection, numProbes, intensity, channel, dyeBias (incl.
 * the GCT bisulfite-conversion metric), betas.
 * Two faithful-but-noted differences live here (see NUMERICS.md):
 *   - detection runs pOOBAH internally (same mask lineage as P), and num_dtna
 *     counts probes with no signal in *either* channel -- the D2 fix means R's
 *     pOOBAH no longer returns NA there;
 *   - the beta group is computed on a D->B->P copy, exactly as R does
 *     (getBetas(pOOBAH(noob(dyeBiasNL(sdf))))).
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present The Children's Hospital of Philadelphia
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "sesame.h"
#include "internal.h"
#include <zlib.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------- TSV formatting --- */

typedef struct { const char *name; int is_int; } qc_desc_t;
#define QC_TAG_I 1
#define QC_TAG_D 0
static const qc_desc_t QC_DESC[] = {
#define QC_ROW(t, nm) { #nm, QC_TAG_##t },
    SESAME_QC_FIELDS(QC_ROW)
#undef QC_ROW
};
enum { QC_N = (int)(sizeof(QC_DESC) / sizeof(QC_DESC[0])) };

static void qc_values(const sesame_qc_t *q, double *v)
{
    int i = 0;
#define QC_ROW(t, nm) v[i++] = q->nm;
    SESAME_QC_FIELDS(QC_ROW)
#undef QC_ROW
}

#define QC_ROW(t, nm) "\t" #nm
static const char QC_HEADER_S[] = SESAME_QC_FIELDS(QC_ROW);
#undef QC_ROW

const char *sesame_qc_header(void) { return QC_HEADER_S + 1; } /* skip leading tab */

int sesame_qc_format_row(const sesame_qc_t *q, char *buf, size_t n)
{
    double v[QC_N];
    size_t off = 0;
    int i;
    qc_values(q, v);
    for (i = 0; i < QC_N; i++) {
        const char *sep = i ? "\t" : "";
        size_t room = n > off ? n - off : 0;
        int w;
        if (isnan(v[i]))            w = snprintf(buf + off, room, "%sNA", sep);
        else if (QC_DESC[i].is_int) w = snprintf(buf + off, room, "%s%lld", sep,
                                                 (long long)llround(v[i]));
        else                        w = snprintf(buf + off, room, "%s%.10g", sep, v[i]);
        if (w < 0) return -1;
        off += (size_t)w;
        if (off >= n) return -1;                       /* truncated */
    }
    return (int)off;
}

/* --------------------------------------------------------------- helpers --- */

static int starts2(const char *s, char a, char b) { return s[0] == a && s[1] == b; }
#define NA1(x) (isnan(x) ? 1 : 0)

#define sigdf_dup sesame_sigdf_dup   /* was a local static; now the public one */

/* Beta-distribution stats over the subset selected by pt (NULL = all probes,
 * else a 2-char probe-type prefix). scr holds >= n doubles. */
static void beta_stats(const double *b, const sesame_sigdf_t *s, const char *pt,
                       double *scr, double *mean, double *median,
                       double *funm, double *fmet, double *numna, double *fracna)
{
    int32_t i, len = 0, nn = 0, na = 0, lt = 0, gt = 0;
    long double sum = 0.0L;
    for (i = 0; i < s->n; i++) {
        if (pt && !starts2(sesame_index_probe_id(s->ix, i), pt[0], pt[1])) continue;
        len++;
        if (isnan(b[i])) { na++; continue; }
        scr[nn++] = b[i]; sum += (long double)b[i];
        if (b[i] < 0.3) lt++;
        if (b[i] > 0.7) gt++;
    }
    *numna  = na;
    *fracna = len ? (double)na / (double)len : NAN;
    *mean   = nn ? (double)(sum / (long double)nn) : NAN;
    if (nn) { sesame__sort(scr, nn); *median = sesame__median_sorted(scr, nn); }
    else    *median = NAN;
    *funm = nn ? (double)lt / (double)nn : NAN;
    *fmet = nn ? (double)gt / (double)nn : NAN;
}

/* ----------------------------------------------------------------- calc --- */

/* Load the positional Type-I extension-base codes (<plat>.typeI_ext.tsv.gz:
 * header then one token per probe, C/T/. -> 1/2/0) for the GCT metric. */
int sesame_load_ext_codes(const char *path, int32_t np, uint8_t **out,
                          sesame_err_t *err)
{
    gzFile f = gzopen(path, "rb");
    char line[256]; uint8_t *e; int32_t row = 0;
    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!f) return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", path);
    if (!(e = (uint8_t *)malloc((size_t)np))) { gzclose(f); return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }
    if (!gzgets(f, line, sizeof line)) { free(e); gzclose(f); return sesame__fail(err, SESAME_ERR_FORMAT, "empty %s", path); }
    while (gzgets(f, line, sizeof line)) {
        if (row < np) e[row] = (line[0]=='C') ? 1 : (line[0]=='T') ? 2 : 0;
        row++;
    }
    gzclose(f);
    if (row != np) { free(e); return sesame__fail(err, SESAME_ERR_FORMAT,
        "%s has %d rows, ordering has %d -- lineage mismatch", path, row, np); }
    *out = e;
    return SESAME_OK;
}

int sesame_qc_calc(const sesame_sigdf_t *s, const uint8_t *bgmask, int32_t bn,
                   const uint8_t *ext, int32_t extn,
                   sesame_qc_t *out, sesame_err_t *err)
{
    int32_t n, i;
    const sesame_index_t *ix;
    double *pv = NULL, *bgG = NULL, *bgR = NULL, *tmp = NULL;
    double *betas = NULL, *scr = NULL, *dR = NULL, *dG = NULL;
    sesame_sigdf_t *cp = NULL;
    int rc = SESAME_OK;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!s || !bgmask || !out) return sesame__fail(err, SESAME_ERR_IO, "null argument");
    n = s->n; ix = s->ix;
    if (bn != n)
        return sesame__fail(err, SESAME_ERR_FORMAT,
            "background mask length %d != probe count %d", bn, n);
    memset(out, 0, sizeof *out);

    /* ===== detection: pOOBAH p per probe (same pools as P) ===== */
    bgG = (double *)malloc((size_t)n * 3 * sizeof(double));
    bgR = (double *)malloc((size_t)n * 3 * sizeof(double));
    tmp = (double *)malloc((size_t)n * 3 * sizeof(double));
    pv  = (double *)malloc((size_t)n * sizeof(double));
    if (!bgG || !bgR || !tmp || !pv) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto done; }
    {
        int32_t nG = 0, nR = 0;
        for (i = 0; i < n; i++) {
            if (bgmask[i]) continue;
            if (s->col[i] == SESAME_COL_R)      { bgG[nG++] = s->MG[i]; bgG[nG++] = s->UG[i]; }
            else if (s->col[i] == SESAME_COL_G) { bgR[nR++] = s->MR[i]; bgR[nR++] = s->UR[i]; }
        }
        for (i = 0; i < n; i++)
            if (sesame__is_neg_control(sesame_index_probe_id(ix, i))) {
                bgG[nG++] = s->UG[i]; bgR[nR++] = s->UR[i];
            }
        nG = sesame__drop_na(bgG, nG, tmp); memcpy(bgG, tmp, (size_t)nG * sizeof(double));
        nR = sesame__drop_na(bgR, nR, tmp); memcpy(bgR, tmp, (size_t)nR * sizeof(double));
        if (nG <= 100) { for (i = 0; i < 1000; i++) bgG[i] = i + 1; nG = 1000; }
        if (nR <= 100) { for (i = 0; i < 1000; i++) bgR[i] = i + 1; nR = 1000; }
        sesame__sort(bgG, nG);
        sesame__sort(bgR, nR);

        for (i = 0; i < n; i++) {
            double rm = sesame__pmax2_narm(s->MR[i], s->UR[i]);
            double gm = sesame__pmax2_narm(s->MG[i], s->UG[i]);
            double pr = isnan(rm) ? NAN : 1.0 - (double)sesame__count_le(bgR, nR, rm) / (double)nR;
            double pg = isnan(gm) ? NAN : 1.0 - (double)sesame__count_le(bgG, nG, gm) / (double)nG;
            if (isnan(pr) && isnan(pg)) pv[i] = NAN;    /* no signal either channel */
            else if (isnan(pr))         pv[i] = pg;
            else if (isnan(pg))         pv[i] = pr;
            else                        pv[i] = pr < pg ? pr : pg;
        }
    }
    free(bgG); bgG = NULL; free(bgR); bgR = NULL; free(tmp); tmp = NULL;
    {
        int32_t dtna = 0, dt = 0, ndef, dt_mk = 0, nmk = 0;
        int32_t dt_t[3] = {0,0,0}, len_t[3] = {0,0,0};
        static const char pts[3][2] = { {'c','g'}, {'c','h'}, {'r','s'} };
        for (i = 0; i < n; i++) {
            const char *id;
            int t;
            if (isnan(pv[i])) { dtna++; continue; }
            if (pv[i] <= 0.05) dt++;
            if (!s->mask[i]) { nmk++; if (pv[i] <= 0.05) dt_mk++; }
            id = sesame_index_probe_id(ix, i);
            for (t = 0; t < 3; t++) if (starts2(id, pts[t][0], pts[t][1])) {
                len_t[t]++; if (pv[i] <= 0.05) dt_t[t]++;
            }
        }
        ndef = n - dtna;
        out->num_dtna = dtna;       out->frac_dtna = (double)dtna / (double)n;
        out->num_dt = dt;           out->frac_dt = ndef ? (double)dt / (double)ndef : NAN;
        out->num_dt_mk = dt_mk;     out->frac_dt_mk = nmk ? (double)dt_mk / (double)nmk : NAN;
        out->num_dt_cg = dt_t[0];   out->frac_dt_cg = len_t[0] ? (double)dt_t[0]/(double)len_t[0] : NAN;
        out->num_dt_ch = dt_t[1];   out->frac_dt_ch = len_t[1] ? (double)dt_t[1]/(double)len_t[1] : NAN;
        out->num_dt_rs = dt_t[2];   out->frac_dt_rs = len_t[2] ? (double)dt_t[2]/(double)len_t[2] : NAN;
    }
    free(pv); pv = NULL;

    /* ===== numProbes ===== */
    {
        int32_t nII = 0, nIR = 0, nIG = 0, ncg = 0, nch = 0, nrs = 0;
        for (i = 0; i < n; i++) {
            const char *id = sesame_index_probe_id(ix, i);
            if (s->col[i] == SESAME_COL_II) nII++;
            else if (s->col[i] == SESAME_COL_R) nIR++;
            else nIG++;
            if      (starts2(id, 'c', 'g')) ncg++;
            else if (starts2(id, 'c', 'h')) nch++;
            else if (starts2(id, 'r', 's')) nrs++;
        }
        out->num_probes = n; out->num_probes_II = nII;
        out->num_probes_IR = nIR; out->num_probes_IG = nIG;
        out->num_probes_cg = ncg; out->num_probes_ch = nch; out->num_probes_rs = nrs;
    }

    /* ===== intensity ===== */
    {
        long double smi=0,smu=0,sii=0,sibg=0,sibr=0,sobg=0,sobr=0;
        int32_t cmi=0,cmu=0,cii=0,cibg=0,cibr=0,cobg=0,cobr=0;
        int32_t naM=0,naU=0,naig=0,nair=0,naii=0;
        for (i = 0; i < n; i++) {
            double MG=s->MG[i],MR=s->MR[i],UG=s->UG[i],UR=s->UR[i];
            double M,U,t;
            if (s->col[i]==SESAME_COL_G)      { M=MG; U=UG; }
            else if (s->col[i]==SESAME_COL_R) { M=MR; U=UR; }
            else                              { M=UG; U=UR; }
            if (!s->mask[i]) {                              /* signalMU(mask=TRUE) */
                if (!isnan(M)) { smi+=M; cmi++; } else naM++;
                if (!isnan(U)) { smi+=U; cmi++; } else naU++;
            }
            t = M + U;
            if (!isnan(t)) { smu += t; cmu++; }             /* totalIntensities */
            if (s->col[i]==SESAME_COL_II) {
                if (!isnan(UG)) { sii+=UG; cii++; }
                if (!isnan(UR)) { sii+=UR; cii++; }
                naii += NA1(UG) + NA1(UR);
            } else if (s->col[i]==SESAME_COL_G) {
                if (!isnan(MG)) { sibg+=MG; cibg++; }
                if (!isnan(UG)) { sibg+=UG; cibg++; }
                if (!isnan(MR)) { sobr+=MR; cobr++; }
                if (!isnan(UR)) { sobr+=UR; cobr++; }
                naig += NA1(MG)+NA1(MR)+NA1(UG)+NA1(UR);
            } else {                                        /* col R */
                if (!isnan(MR)) { sibr+=MR; cibr++; }
                if (!isnan(UR)) { sibr+=UR; cibr++; }
                if (!isnan(MG)) { sobg+=MG; cobg++; }
                if (!isnan(UG)) { sobg+=UG; cobg++; }
                nair += NA1(MG)+NA1(MR)+NA1(UG)+NA1(UR);
            }
        }
        out->mean_intensity    = cmi ? (double)(smi/cmi) : NAN;
        out->mean_intensity_MU = cmu ? (double)(smu/cmu) : NAN;
        out->mean_ii      = cii  ? (double)(sii/cii)   : NAN;
        out->mean_inb_grn = cibg ? (double)(sibg/cibg) : NAN;
        out->mean_inb_red = cibr ? (double)(sibr/cibr) : NAN;
        out->mean_oob_grn = cobg ? (double)(sobg/cobg) : NAN;
        out->mean_oob_red = cobr ? (double)(sobr/cobr) : NAN;
        out->na_intensity_M = naM; out->na_intensity_U = naU;
        out->na_intensity_ig = naig; out->na_intensity_ir = nair; out->na_intensity_ii = naii;
    }

    /* ===== channel: run inferInfiniumIChannel on a copy, count transitions ===== */
    cp = sigdf_dup(s);
    if (!cp) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto done; }
    if ((rc = sesame_prep_infer_channel(cp, 0, 0, err)) != SESAME_OK) goto done;
    {
        int32_t r2r=0,g2g=0,r2g=0,g2r=0;
        for (i = 0; i < n; i++) {
            if (s->col[i] == SESAME_COL_II) continue;
            if (s->col[i]==SESAME_COL_R && cp->col[i]==SESAME_COL_R) r2r++;
            else if (s->col[i]==SESAME_COL_G && cp->col[i]==SESAME_COL_G) g2g++;
            else if (s->col[i]==SESAME_COL_R && cp->col[i]==SESAME_COL_G) r2g++;
            else if (s->col[i]==SESAME_COL_G && cp->col[i]==SESAME_COL_R) g2r++;
        }
        out->InfI_switch_R2R = r2r; out->InfI_switch_G2G = g2g;
        out->InfI_switch_R2G = r2g; out->InfI_switch_G2R = g2r;
    }
    sesame_sigdf_free(cp); cp = NULL;

    /* ===== dyeBias ===== */
    dR = (double *)malloc((size_t)n * sizeof(double));
    dG = (double *)malloc((size_t)n * sizeof(double));
    if (!dR || !dG) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto done; }
    {
        int32_t nR=0,nG=0;
        double medR,medG,topR,topG;
        for (i = 0; i < n; i++) {
            if (s->col[i]==SESAME_COL_R) { double v=s->MR[i]+s->UR[i]; if(!isnan(v)) dR[nR++]=v; }
            else if (s->col[i]==SESAME_COL_G) { double v=s->MG[i]+s->UG[i]; if(!isnan(v)) dG[nG++]=v; }
        }
        sesame__sort(dR, nR); sesame__sort(dG, nG);
        medR = nR ? sesame__median_sorted(dR, nR) : NAN;
        medG = nG ? sesame__median_sorted(dG, nG) : NAN;
        topR = nR ? sesame__median_sorted(dR + (nR>20?nR-20:0), nR>20?20:nR) : NAN;
        topG = nG ? sesame__median_sorted(dG + (nG>20?nG-20:0), nG>20?20:nG) : NAN;
        out->medR=medR; out->medG=medG; out->topR=topR; out->topG=topG;
        out->RGratio = medR/medG;
        out->RGdistort = (topR/topG)/(medR/medG);
    }
    free(dR); dR = NULL; free(dG); dG = NULL;

    /* ===== GCT (bisConversionControl): mean green of Type-I extension-C probes
     * over that of extension-A/T probes (R/sesame.R:562). NaN without ext. ===== */
    {
        double sumC=0.0, sumT=0.0; int32_t cC=0, cT=0;
        if (ext && extn == n) for (i = 0; i < n; i++) {
            if (ext[i]==1) { if(!isnan(s->MG[i])){sumC+=s->MG[i];cC++;} if(!isnan(s->UG[i])){sumC+=s->UG[i];cC++;} }
            else if (ext[i]==2) { if(!isnan(s->MG[i])){sumT+=s->MG[i];cT++;} if(!isnan(s->UG[i])){sumT+=s->UG[i];cT++;} }
        }
        out->GCT = (cC && cT) ? (sumC/(double)cC)/(sumT/(double)cT) : NAN;
    }

    /* ===== betas: getBetas(pOOBAH(noob(dyeBiasNL(copy)))) ===== */
    cp = sigdf_dup(s);
    betas = (double *)malloc((size_t)n * sizeof(double));
    scr   = (double *)malloc((size_t)n * sizeof(double));
    if (!cp || !betas || !scr) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto done; }
    if ((rc = sesame_prep_dye_bias_nl(cp, err)) != SESAME_OK) goto done;
    if ((rc = sesame_prep_noob(cp, bgmask, bn, 1, 15.0, err)) != SESAME_OK) goto done;
    if ((rc = sesame_prep_poobah(cp, bgmask, bn, 0.05, 1, err)) != SESAME_OK) goto done;
    if ((rc = sesame_get_betas(cp, 1, betas, err)) != SESAME_OK) goto done;
    beta_stats(betas, s, NULL, scr, &out->mean_beta, &out->median_beta,
               &out->frac_unmeth, &out->frac_meth, &out->num_na, &out->frac_na);
    beta_stats(betas, s, "cg", scr, &out->mean_beta_cg, &out->median_beta_cg,
               &out->frac_unmeth_cg, &out->frac_meth_cg, &out->num_na_cg, &out->frac_na_cg);
    beta_stats(betas, s, "ch", scr, &out->mean_beta_ch, &out->median_beta_ch,
               &out->frac_unmeth_ch, &out->frac_meth_ch, &out->num_na_ch, &out->frac_na_ch);
    beta_stats(betas, s, "rs", scr, &out->mean_beta_rs, &out->median_beta_rs,
               &out->frac_unmeth_rs, &out->frac_meth_rs, &out->num_na_rs, &out->frac_na_rs);

done:
    free(bgG); free(bgR); free(tmp); free(pv);
    free(dR); free(dG); free(betas); free(scr);
    sesame_sigdf_free(cp);
    return rc;
}
