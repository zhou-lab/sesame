/* index.c -- the per-platform ordering table.
 *
 * Reads the ordering table in either shape:
 *     Probe_ID <SEP> M <SEP> U <SEP> col [ <SEP> mask ]
 * where SEP is a tab or a comma, sniffed from the header line -- Illumina
 * manifests are CSV, and a custom array's ordering is usually derived from one,
 * so insisting on TSV would make everyone convert for no reason. The sniff is
 * per file, not per line: whichever separator the header uses is the one used
 * throughout, so a Probe_ID containing the other character stays safe.
 * with M/U as decimal bead addresses or "NA", col in {G,R,2}, mask in {0,1}.
 * The 5-column form is the legacy tools/export_ordering.R (sesameData) export
 * with the mask inline; the 4-column form is the published InfiniumAnnotation
 * ordering, where the mask has moved to the companion .cm file (mask -> 0 here).
 *
 * The whole (decompressed) file is slurped and parsed in place -- for ~937k
 * rows that is far cheaper than per-line I/O, and the blob doubles as the
 * Probe_ID string storage.
 *
 * The struct is opaque so this format stays an implementation detail; a
 * mmap-able binary can replace it without touching any caller.
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
#include <zlib.h>

#define MAX_INDEX_BYTES ((size_t)1 << 30)  /* 1 GiB slurp cap */

struct sesame_index_t {
    int32_t   n;
    char     *blob;    /* Probe_ID chars, NUL-terminated in place */
    uint32_t *id_off;  /* n offsets into blob */
    uint32_t *M;       /* 0 == NA (bead addresses are never 0) */
    uint32_t *U;
    uint8_t  *col;     /* SESAME_COL_* */
    uint8_t  *mask;
};

int32_t sesame_index_nprobes(const sesame_index_t *ix)
{
    return ix ? ix->n : 0;
}

const char *sesame_index_probe_id(const sesame_index_t *ix, int32_t i)
{
    if (!ix || i < 0 || i >= ix->n) return NULL;
    return ix->blob + ix->id_off[i];
}

/* Accessors for other TUs (sigdf.c) without exposing the layout publicly. */
const uint32_t *sesame__index_M(const sesame_index_t *ix)    { return ix->M; }
const uint32_t *sesame__index_U(const sesame_index_t *ix)    { return ix->U; }
const uint8_t  *sesame__index_col(const sesame_index_t *ix)  { return ix->col; }
const uint8_t  *sesame__index_mask(const sesame_index_t *ix) { return ix->mask; }

/* Slurp a possibly-gzipped file. Returns malloc'd buffer, sets *len. */
static char *slurp(const char *path, size_t *len, sesame_err_t *err)
{
    gzFile f = gzopen(path, "rb");
    size_t cap = 1 << 22, used = 0;
    char *buf;

    if (!f) { sesame__fail(err, SESAME_ERR_IO, "cannot open %s", path); return NULL; }
    buf = (char *)malloc(cap);
    if (!buf) { gzclose(f); sesame__fail(err, SESAME_ERR_NOMEM, "oom"); return NULL; }

    for (;;) {
        int got;
        if (used + (1 << 20) + 1 > cap) {
            char *nb;
            if (cap > MAX_INDEX_BYTES) {
                free(buf); gzclose(f);
                sesame__fail(err, SESAME_ERR_FORMAT, "index %s too large", path);
                return NULL;
            }
            cap *= 2;
            nb = (char *)realloc(buf, cap);
            if (!nb) { free(buf); gzclose(f);
                       sesame__fail(err, SESAME_ERR_NOMEM, "oom"); return NULL; }
            buf = nb;
        }
        got = gzread(f, buf + used, 1 << 20);
        if (got < 0) { free(buf); gzclose(f);
                       sesame__fail(err, SESAME_ERR_IO, "read error on %s", path);
                       return NULL; }
        if (got == 0) break;
        used += (size_t)got;
    }
    gzclose(f);
    buf[used] = '\0';
    *len = used;
    return buf;
}

/* Parse a decimal uint from [p,end); "NA" -> 0. Returns 0 on success. */
static int parse_addr(const char *p, const char *end, uint32_t *out)
{
    uint64_t v = 0;
    if (p == end) return -1;
    if (end - p == 2 && p[0] == 'N' && p[1] == 'A') { *out = 0; return 0; }
    for (; p < end; p++) {
        if (*p < '0' || *p > '9') return -1;
        v = v * 10 + (uint64_t)(*p - '0');
        if (v > 0xffffffffu) return -1;
    }
    *out = (uint32_t)v;
    return 0;
}

sesame_index_t *sesame_index_open(const char *path, sesame_err_t *err)
{
    sesame_index_t *ix = NULL;
    char *buf = NULL;
    size_t len = 0;
    int32_t n = 0, row = 0;
    char *p, *eof, sep;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!path) { sesame__fail(err, SESAME_ERR_IO, "null path"); return NULL; }

    buf = slurp(path, &len, err);
    if (!buf) return NULL;
    eof = buf + len;

    /* Count rows (lines) to size the arrays in one pass. */
    for (p = buf; p < eof; p++) if (*p == '\n') n++;
    if (len > 0 && eof[-1] != '\n') n++;   /* unterminated last line */
    n--;                                   /* header */
    if (n <= 0) {
        free(buf);
        sesame__fail(err, SESAME_ERR_FORMAT, "%s: no data rows", path);
        return NULL;
    }

    ix = (sesame_index_t *)calloc(1, sizeof(*ix));
    if (!ix) { free(buf); sesame__fail(err, SESAME_ERR_NOMEM, "oom"); return NULL; }
    ix->blob   = buf;
    ix->id_off = (uint32_t *)malloc((size_t)n * sizeof(uint32_t));
    ix->M      = (uint32_t *)malloc((size_t)n * sizeof(uint32_t));
    ix->U      = (uint32_t *)malloc((size_t)n * sizeof(uint32_t));
    ix->col    = (uint8_t  *)malloc((size_t)n);
    ix->mask   = (uint8_t  *)malloc((size_t)n);
    if (!ix->id_off || !ix->M || !ix->U || !ix->col || !ix->mask) {
        sesame_index_close(ix);
        sesame__fail(err, SESAME_ERR_NOMEM, "oom sizing %d probes", n);
        return NULL;
    }

    /* Sniff the separator off the header, then skip it. A comma only wins if
     * the header has no tab, so a tab-separated file with commas inside a field
     * is still read as TSV. */
    p = buf;
    while (p < eof && *p != '\n' && *p != '\t' && *p != ',') p++;
    sep = (p < eof && *p == ',') ? ',' : '\t';
    while (p < eof && *p != '\n') p++;
    if (p < eof) p++;

    while (p < eof && row < n) {
        char *line = p, *le, *f[5];
        int nf = 0;

        le = line;
        while (le < eof && *le != '\n') le++;
        if (le == line) { p = le + 1; continue; }   /* blank line */

        /* Split on the sniffed separator, in place. */
        f[nf++] = line;
        for (char *q = line; q < le && nf < 5; q++) {
            if (*q == sep) { *q = '\0'; f[nf++] = q + 1; }
        }
        /* 5 columns = legacy sesameData export (mask inline); 4 columns = the
         * published InfiniumAnnotation ordering, where the mask now lives in the
         * companion .cm file, so there is no mask column here (mask defaults 0). */
        if (nf != 4 && nf != 5) {
            sesame_index_close(ix);
            sesame__fail(err, SESAME_ERR_FORMAT,
                         "%s: row %d has %d fields, expected 4 or 5", path, row + 1, nf);
            return NULL;
        }
        /* Terminate the final field (overwrites '\n', or the NUL at EOF). */
        *le = '\0';

        ix->id_off[row] = (uint32_t)(f[0] - ix->blob);

        if (parse_addr(f[1], f[1] + strlen(f[1]), &ix->M[row]) ||
            parse_addr(f[2], f[2] + strlen(f[2]), &ix->U[row])) {
            sesame_index_close(ix);
            sesame__fail(err, SESAME_ERR_FORMAT,
                         "%s: row %d has a bad address", path, row + 1);
            return NULL;
        }

        if      (f[3][0] == 'G' && f[3][1] == '\0') ix->col[row] = SESAME_COL_G;
        else if (f[3][0] == 'R' && f[3][1] == '\0') ix->col[row] = SESAME_COL_R;
        else if (f[3][0] == '2' && f[3][1] == '\0') ix->col[row] = SESAME_COL_II;
        else {
            sesame_index_close(ix);
            sesame__fail(err, SESAME_ERR_FORMAT,
                         "%s: row %d has bad col '%s' (want G, R or 2)",
                         path, row + 1, f[3]);
            return NULL;
        }

        ix->mask[row] = (uint8_t)(nf == 5 && f[4][0] == '1');
        row++;
        p = le + 1;
    }

    if (row != n) {
        sesame_index_close(ix);
        sesame__fail(err, SESAME_ERR_FORMAT,
                     "%s: parsed %d rows, expected %d", path, row, n);
        return NULL;
    }

    ix->n = n;
    return ix;
}

void sesame_index_close(sesame_index_t *ix)
{
    if (!ix) return;
    free(ix->blob);
    free(ix->id_off);
    free(ix->M);
    free(ix->U);
    free(ix->col);
    free(ix->mask);
    free(ix);
}
