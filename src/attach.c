/* attach.c -- attach the ordering's Probe_IDs to a positional data file.
 *
 * A YAME .cg/.cm/.cx stores one value per probe in ordering order, with NO probe
 * id inside the container (row names live in the ordering, not the data). The
 * per-probe genomic coordinate tables (<platform>.hg38.coord.tsv.gz) are the
 * same: positional, no Probe_ID column. This turns either into a labeled TSV by
 * pairing row i with the ordering's i-th Probe_ID -- so the output is directly
 * greppable / joinable. The lineage must match (same platform + tag that
 * produced the file); a row-count mismatch is a hard error, not silent
 * misalignment.
 *
 * The YAME per-format rendering mirrors `yame unpack`'s print_cdata1 (fmt0 mask
 * bit, fmt3 M/U or beta, fmt4 float, fmt5 ternary, fmt1/2 raw), except floats
 * print at full precision rather than 3 decimals since the .cg already stores
 * float32. Links YAME directly, same as cgwrite.c.
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

#include "cfile.h"    /* open_cfile, read_cdata1 */
#include "cdata.h"    /* cdata_t, decompress, free_cdata, f3_get_mu, f2_get_string */
#include "index.h"    /* loadSampleNamesFromIndex, cleanSampleNames2, snames_t */

#include <ctype.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

static int has_suffix(const char *s, const char *suf)
{
    size_t ls = strlen(s), lf = strlen(suf);
    return ls >= lf && strcmp(s + ls - lf, suf) == 0;
}

/* A YAME container, by extension. Everything else is treated as text. */
static int is_yame(const char *path)
{
    return has_suffix(path, ".cg") || has_suffix(path, ".cm") ||
           has_suffix(path, ".cx") || has_suffix(path, ".cr");
}

/* Read one logical (possibly long) line into the grown buffer, stripping the
 * trailing newline. Returns 1 on a line, 0 at EOF, -1 on OOM. */
static int read_gzline(gzFile f, char **buf, size_t *cap)
{
    size_t len = 0;
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

/* The ordering columns requested with --with, in the ordering's own column
 * order, so the output reads like the ordering it came from. */
static void with_head(FILE *out, const sesame_attach_opt_t *opt)
{
    if (opt->with & SESAME_WITH_M)    fputs("\tM", out);
    if (opt->with & SESAME_WITH_U)    fputs("\tU", out);
    if (opt->with & SESAME_WITH_COL)  fputs("\tcol", out);
    if (opt->with & SESAME_WITH_MASK) fputs("\tmask", out);
}

static void with_row(FILE *out, const sesame_index_t *ix, int32_t i,
                     const sesame_attach_opt_t *opt)
{
    if (opt->with & SESAME_WITH_M) {
        uint32_t v = sesame__index_M(ix)[i];
        if (v) fprintf(out, "\t%u", v); else fputs("\tNA", out);
    }
    if (opt->with & SESAME_WITH_U) {
        uint32_t v = sesame__index_U(ix)[i];
        if (v) fprintf(out, "\t%u", v); else fputs("\tNA", out);
    }
    if (opt->with & SESAME_WITH_COL) {
        uint8_t c = sesame__index_col(ix)[i];
        fprintf(out, "\t%s", c == SESAME_COL_G ? "G" :
                              c == SESAME_COL_R ? "R" : "2");
    }
    if (opt->with & SESAME_WITH_MASK)
        fprintf(out, "\t%u", (unsigned)sesame__index_mask(ix)[i]);
}

/* Count data rows (total lines minus the header, unless no_header). */
static int count_text_rows(const char *path, int no_header, int32_t *nrow,
                           sesame_err_t *err)
{
    gzFile f = gzopen(path, "rb");
    char *buf; size_t cap = 1 << 16;
    long lines = 0; int r;
    if (!f) return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", path);
    if (!(buf = (char *)malloc(cap))) { gzclose(f);
        return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }
    while ((r = read_gzline(f, &buf, &cap)) == 1) lines++;
    free(buf); gzclose(f);
    if (r < 0) return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
    *nrow = (int32_t)(no_header ? lines : (lines > 0 ? lines - 1 : 0));
    return SESAME_OK;
}

/* ------------------------------------------------------------- text ------
 *
 * gzopen reads plain and gzipped text alike. The file's own header is kept and
 * prefixed with "Probe_ID"; each data line gets its positional Probe_ID. Row
 * count is checked up front so a lineage mismatch fails before any output. */
static int attach_text(const char *path, const sesame_index_t *ix,
                       const sesame_attach_opt_t *opt, FILE *out,
                       sesame_err_t *err)
{
    gzFile f;
    char *buf = NULL;
    size_t cap = 1 << 16;
    int32_t nid = sesame_index_nprobes(ix), nrow = 0, row = 0;
    int first = 1, r;

    if (count_text_rows(path, opt->no_header, &nrow, err) != SESAME_OK)
        return err ? err->code : SESAME_ERR_IO;
    if (nrow != nid)
        return sesame__fail(err, SESAME_ERR_FORMAT,
            "%s has %d data rows, ordering has %d -- lineage mismatch",
            path, nrow, nid);

    if (!(f = gzopen(path, "rb")))
        return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", path);
    if (!(buf = (char *)malloc(cap))) { gzclose(f);
        return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }

    while ((r = read_gzline(f, &buf, &cap)) == 1) {
        if (first && !opt->no_header) {
            fputs("Probe_ID", out); with_head(out, opt);
            fprintf(out, "\t%s\n", buf);
            first = 0;
            continue;
        }
        first = 0;
        fputs(sesame_index_probe_id(ix, row), out);
        with_row(out, ix, row, opt);
        fprintf(out, "\t%s\n", buf);
        row++;
    }
    free(buf);
    gzclose(f);
    if (r < 0) return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
    return SESAME_OK;
}

/* ------------------------------------------------------------- YAME ------ */

/* One value of decompressed record d at row i. Mirrors yame unpack. */
static void render1(FILE *out, cdata_t *d, uint64_t i,
                    const sesame_attach_opt_t *opt)
{
    switch (d->fmt) {
    case '0':                              /* mask bit */
        fputc(((d->s[i>>3] >> (i&0x7)) & 0x1) + '0', out);
        break;
    case '1':                              /* raw byte */
        fputc(d->s[i], out);
        break;
    case '2':                              /* state label */
        fputs(f2_get_string(d, i), out);
        break;
    case '3': {                            /* M/U counts */
        uint64_t mu = f3_get_mu(d, i);
        uint64_t M = mu >> 32, U = (mu << 32) >> 32;
        if (opt->beta) {
            if (M == 0 && U == 0) fputs("NA", out);
            else fprintf(out, "%.6g", (double)M / (double)(M + U));
        } else {
            fprintf(out, "%" PRIu64 "\t%" PRIu64, M, U);
        }
        break;
    }
    case '4': {                            /* float32, negative = NA */
        float v = ((float *)d->s)[i];
        if (v < 0) fputs("NA", out); else fprintf(out, "%.6g", (double)v);
        break;
    }
    case '5':                              /* ternary, 2 = NA */
        if (d->s[i] == 2) fputs("NA", out); else fputc(d->s[i] + '0', out);
        break;
    default:
        fputs("NA", out);
        break;
    }
}

/* Header cell(s) for a sample under format fmt (fmt3 M/U is two columns). */
static void render_head(FILE *out, char fmt, const char *name,
                        const sesame_attach_opt_t *opt)
{
    if (fmt == '3' && !opt->beta) fprintf(out, "%s_M\t%s_U", name, name);
    else                          fputs(name, out);
}

/* Several .cg files, side by side.
 *
 * They are all positional to the same ordering, so pairing them is a column
 * concatenation -- there is no key to match on. Doing it here rather than
 * making the caller write out two labelled TSVs and join them back on their
 * Probe_ID is both shorter and honest about what the operation is: the
 * files already line up, row i is probe i in every one of them. The row
 * counts must agree, and a mismatch is the lineage error it always was. */
static int attach_yame(const char *const *paths, int npath,
                       const sesame_index_t *ix,
                       const sesame_attach_opt_t *opt, FILE *out,
                       sesame_err_t *err)
{
    cfile_t cf;
    snames_t sn;
    cdata_t *recs = NULL;
    char **names = NULL;                     /* per record, owned */
    char pathbuf[4096];
    int32_t nid = sesame_index_nprobes(ix), nrec = 0, cap = 0, j;
    int64_t np = -1, i;
    int rc = SESAME_OK, k;

    for (k = 0; k < npath; k++) {
        int32_t got = 0;
        snprintf(pathbuf, sizeof pathbuf, "%s", paths[k]);
        cf = open_cfile(pathbuf);
        if (!cf.fh) {
            rc = sesame__fail(err, SESAME_ERR_IO, "cannot open %s", paths[k]);
            goto done;
        }
        sn = loadSampleNamesFromIndex(pathbuf);

        for (;;) {
            cdata_t c = read_cdata1(&cf), d;
            if (c.n == 0) break;                        /* EOF */
            d = decompress(c);
            free_cdata(&c);
            if (np < 0) np = (int64_t)d.n;
            else if ((int64_t)d.n != np) {
                free_cdata(&d);
                cleanSampleNames2(sn); bgzf_close(cf.fh);
                rc = sesame__fail(err, SESAME_ERR_FORMAT,
                    "%s has %" PRIu64 " rows, an earlier file has %" PRId64
                    " -- they are not the same ordering", paths[k],
                    (uint64_t)d.n, np);
                goto done;
            }
            if (nrec >= cap) {
                cap = cap ? cap * 2 : 8;
                recs  = (cdata_t *)realloc(recs, (size_t)cap * sizeof(cdata_t));
                names = (char **)realloc(names, (size_t)cap * sizeof(char *));
            }
            {   const char *nm = (got < sn.n) ? sn.s[got] : NULL;
                names[nrec] = strdup((nm && *nm) ? nm : "V"); }
            recs[nrec++] = d;
            got++;
            if (!opt->all) break;                       /* first record only */
        }
        cleanSampleNames2(sn);
        bgzf_close(cf.fh);
        if (got == 0) {
            rc = sesame__fail(err, SESAME_ERR_FORMAT, "no records in %s", paths[k]);
            goto done;
        }
    }

    if (nrec == 0) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT, "no records in %s", paths[0]);
        goto done;
    }
    if (np != nid) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT,
            "%s has %" PRId64 " probes, ordering has %d -- lineage mismatch",
            paths[0], np, nid);
        goto done;
    }

    if (!opt->no_header) {
        fputs("Probe_ID", out);
        with_head(out, opt);
        for (j = 0; j < nrec; j++) {
            fputc('\t', out);
            render_head(out, recs[j].fmt, names[j], opt);
        }
        fputc('\n', out);
    }
    for (i = 0; i < np; i++) {
        fputs(sesame_index_probe_id(ix, (int32_t)i), out);
        with_row(out, ix, (int32_t)i, opt);
        for (j = 0; j < nrec; j++) {
            fputc('\t', out);
            render1(out, &recs[j], (uint64_t)i, opt);
        }
        fputc('\n', out);
    }

done:
    for (j = 0; j < nrec; j++) { free_cdata(&recs[j]); free(names[j]); }
    free(recs); free(names);
    return rc;
}

int sesame_attach_probe(const char *path, const sesame_index_t *ix,
                        const sesame_attach_opt_t *opt, FILE *out,
                        sesame_err_t *err)
{
    return sesame_attach_probe_n(&path, 1, ix, opt, out, err);
}

int sesame_attach_probe_n(const char *const *paths, int npath,
                          const sesame_index_t *ix,
                          const sesame_attach_opt_t *opt, FILE *out,
                          sesame_err_t *err)
{
    static const sesame_attach_opt_t deflt = { 0, 0, 0, 0 };
    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!opt) opt = &deflt;
    if (npath < 1) return sesame__fail(err, SESAME_ERR_IO, "no input file");
    if (!is_yame(paths[0])) {
        if (npath > 1)
            return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
                "several inputs are only supported for YAME files");
        return attach_text(paths[0], ix, opt, out, err);
    }
    return attach_yame(paths, npath, ix, opt, out, err);
}
