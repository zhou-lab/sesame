/* idat.c -- Illumina IDAT v3 reader
 *
 * Ported from sesame's R/readIDAT.R (which was in turn adapted from
 * illuminaio with fixes for de-identified files).
 *
 * Divergences from the R implementation (see NUMERICS.md):
 *   D3 -- R parses the header twice: readIDAT() reads `version` as int32 to
 *         dispatch, then readIDAT_nonenc() reopens the file and re-reads it as
 *         int64 (R/readIDAT.R:22 vs :178). It only works because the value is
 *         small and little-endian. We parse once, as int64.
 *   D4 -- R's readLong() uses readBin(what="integer", size=8), which base R
 *         does not actually support (R/readIDAT.R:44-46). We read a real
 *         int64_t.
 *
 * This code parses untrusted input (IDATs are routinely downloaded from GEO).
 * nFields, byteOffset and nSNPsRead are all attacker-controlled, and the field
 * table is a 10-byte *unpadded* record. Every length is bounds-checked and
 * every allocation is capped before it is made.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

/* Sanity caps. Generous enough for any real array, tight enough that a corrupt
 * header cannot induce a huge allocation. Current largest platform (EPICv2) is
 * ~937k probes; IDATs carry ~1-2M beads. */
#define MAX_FIELDS   4096
#define MAX_PROBES   50000000

#define FC_NSNPSREAD 1000
#define FC_ILLUMINAID 102
#define FC_SD         103
#define FC_MEAN       104
#define FC_NBEADS     107

#define fail sesame__fail

/* gzread() is int-limited and may return short reads; loop until satisfied.
 * Returns 0 on success, -1 on short read/error. */
static int rd_full(gzFile f, void *buf, size_t len)
{
    unsigned char *p = (unsigned char *)buf;
    while (len > 0) {
        unsigned chunk = (len > 0x40000000u) ? 0x40000000u : (unsigned)len;
        int got = gzread(f, p, chunk);
        if (got <= 0) return -1;
        p += got;
        len -= (size_t)got;
    }
    return 0;
}

static int rd_u16(gzFile f, uint16_t *v)
{
    unsigned char b[2];
    if (rd_full(f, b, 2)) return -1;
    *v = (uint16_t)((uint16_t)b[0] | ((uint16_t)b[1] << 8));
    return 0;
}

static int rd_i32(gzFile f, int32_t *v)
{
    unsigned char b[4];
    if (rd_full(f, b, 4)) return -1;
    *v = (int32_t)((uint32_t)b[0] | ((uint32_t)b[1] << 8) |
                   ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24));
    return 0;
}

static int rd_i64(gzFile f, int64_t *v)
{
    unsigned char b[8];
    uint64_t x = 0;
    int i;
    if (rd_full(f, b, 8)) return -1;
    for (i = 7; i >= 0; i--) x = (x << 8) | (uint64_t)b[i];
    *v = (int64_t)x;
    return 0;
}

/* Decode a little-endian u16 block in place. No-op on LE hosts. */
static void fix_u16(uint16_t *a, int32_t n)
{
    const uint16_t probe = 1;
    if (*(const unsigned char *)&probe == 1) return; /* little-endian host */
    for (int32_t i = 0; i < n; i++)
        a[i] = (uint16_t)((a[i] >> 8) | (a[i] << 8));
}

static void fix_u32(uint32_t *a, int32_t n)
{
    const uint16_t probe = 1;
    if (*(const unsigned char *)&probe == 1) return;
    for (int32_t i = 0; i < n; i++) {
        uint32_t v = a[i];
        a[i] = ((v >> 24) & 0xffu) | ((v >> 8) & 0xff00u) |
               ((v << 8) & 0xff0000u) | ((v << 24) & 0xff000000u);
    }
}

/* gzseek from the start. Returns 0 on success. */
static int seek_to(gzFile f, int64_t off)
{
    if (off < 0) return -1;
    if (gzseek(f, (z_off_t)off, SEEK_SET) != (z_off_t)off) return -1;
    return 0;
}

int sesame_idat_read(const char *path, sesame_idat_t **out, sesame_err_t *err)
{
    gzFile f = NULL;
    sesame_idat_t *d = NULL;
    unsigned char magic[4];
    int64_t version = 0;
    int32_t n_fields = 0, n = 0, i;
    int64_t off_n = -1, off_id = -1, off_sd = -1, off_mean = -1, off_nb = -1;
    int rc;

    if (!path || !out) return fail(err, SESAME_ERR_IO, "null argument");
    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }

    /* zlib reads uncompressed files transparently, so one path serves both
     * .idat and .idat.gz -- same as R's gzfile(). */
    f = gzopen(path, "rb");
    if (!f) return fail(err, SESAME_ERR_IO, "cannot open %s", path);

    if (rd_full(f, magic, 4))
        { rc = fail(err, SESAME_ERR_IO, "short read on magic"); goto done; }
    if (memcmp(magic, "IDAT", 4) != 0)
        { rc = fail(err, SESAME_ERR_FORMAT,
                    "invalid IDAT magic (expected 'IDAT')"); goto done; }

    /* D3/D4: single parse, real int64. */
    if (rd_i64(f, &version))
        { rc = fail(err, SESAME_ERR_IO, "short read on version"); goto done; }
    if (version < 3)
        { rc = fail(err, SESAME_ERR_UNSUPPORTED,
                    "unsupported IDAT version %lld (encrypted/v<3 not supported)",
                    (long long)version); goto done; }

    if (rd_i32(f, &n_fields))
        { rc = fail(err, SESAME_ERR_IO, "short read on nFields"); goto done; }
    if (n_fields <= 0 || n_fields > MAX_FIELDS)
        { rc = fail(err, SESAME_ERR_FORMAT, "implausible nFields %d",
                    n_fields); goto done; }

    /* Field table: n_fields x { u16 code; i64 offset }, unpadded (10 bytes). */
    for (i = 0; i < n_fields; i++) {
        uint16_t code;
        int64_t  off;
        if (rd_u16(f, &code) || rd_i64(f, &off))
            { rc = fail(err, SESAME_ERR_IO,
                        "short read in field table at %d", i); goto done; }
        if (off < 0)
            { rc = fail(err, SESAME_ERR_FORMAT,
                        "negative byteOffset in field table at %d", i); goto done; }
        switch (code) {
        case FC_NSNPSREAD:  off_n    = off; break;
        case FC_ILLUMINAID: off_id   = off; break;
        case FC_SD:         off_sd   = off; break;
        case FC_MEAN:       off_mean = off; break;
        case FC_NBEADS:     off_nb   = off; break;
        default: break; /* strings/RunInfo/MidBlock deliberately ignored */
        }
    }

    if (off_n < 0 || off_id < 0 || off_sd < 0 || off_mean < 0 || off_nb < 0)
        { rc = fail(err, SESAME_ERR_FORMAT,
                    "missing required field(s): need 1000/102/103/104/107");
          goto done; }

    if (seek_to(f, off_n) || rd_i32(f, &n))
        { rc = fail(err, SESAME_ERR_IO, "cannot read nSNPsRead"); goto done; }
    if (n <= 0 || n > MAX_PROBES)
        { rc = fail(err, SESAME_ERR_FORMAT, "implausible nSNPsRead %d", n);
          goto done; }

    d = (sesame_idat_t *)calloc(1, sizeof(*d));
    if (!d) { rc = fail(err, SESAME_ERR_NOMEM, "oom"); goto done; }
    d->n = n;
    d->version = version;
    d->n_fields = n_fields;
    d->off_mean = off_mean;
    d->addr   = (uint32_t *)malloc((size_t)n * sizeof(uint32_t));
    d->mean   = (uint16_t *)malloc((size_t)n * sizeof(uint16_t));
    d->sd     = (uint16_t *)malloc((size_t)n * sizeof(uint16_t));
    d->nbeads = (uint8_t  *)malloc((size_t)n * sizeof(uint8_t));
    if (!d->addr || !d->mean || !d->sd || !d->nbeads)
        { rc = fail(err, SESAME_ERR_NOMEM, "oom allocating %d probes", n);
          goto done; }

    if (seek_to(f, off_id) || rd_full(f, d->addr, (size_t)n * 4))
        { rc = fail(err, SESAME_ERR_IO, "cannot read IlluminaID"); goto done; }
    fix_u32(d->addr, n);

    if (seek_to(f, off_sd) || rd_full(f, d->sd, (size_t)n * 2))
        { rc = fail(err, SESAME_ERR_IO, "cannot read SD"); goto done; }
    fix_u16(d->sd, n);

    if (seek_to(f, off_mean) || rd_full(f, d->mean, (size_t)n * 2))
        { rc = fail(err, SESAME_ERR_IO, "cannot read Mean"); goto done; }
    fix_u16(d->mean, n);

    if (seek_to(f, off_nb) || rd_full(f, d->nbeads, (size_t)n))
        { rc = fail(err, SESAME_ERR_IO, "cannot read NBeads"); goto done; }

    gzclose(f);
    *out = d;
    return SESAME_OK;

done:
    if (f) gzclose(f);
    sesame_idat_free(d);
    return rc;
}

void sesame_idat_free(sesame_idat_t *d)
{
    if (!d) return;
    free(d->addr);
    free(d->mean);
    free(d->sd);
    free(d->nbeads);
    free(d);
}
