/* deidentify.c -- strip the genetic fingerprint from an IDAT.
 *
 * Port of sesame's deIdentify/reIdentify (R/deidentify.R). The SNP (rs) probes
 * carry per-individual genotype; setting their bead Mean intensities to zero
 * (or reversibly scrambling them) makes the IDAT unusable for identity matching
 * while leaving every other byte -- and thus every non-SNP measurement --
 * intact. Only the Mean section is rewritten; the rest of the file is copied
 * verbatim, so a de-identified IDAT reads exactly like the original elsewhere.
 *
 * The scramble/restore path uses this tool's own PRNG (splitmix64 + Fisher-Yates
 * keyed by --seed), so a randomized file round-trips with `sesame reidentify
 * --seed`, not with R. Zeroing is fully R-equivalent (and irreversible).
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present The Children's Hospital of Philadelphia
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "sesame.h"
#include "internal.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

/* --- a tiny open-addressing set of bead addresses (0 is never a real one) --- */
typedef struct { uint32_t *slot; int32_t cap; } aset_t;

static int aset_init(aset_t *s, int32_t want)
{
    int32_t cap = 16;
    while (cap < want * 2) cap <<= 1;
    s->slot = (uint32_t *)calloc((size_t)cap, sizeof(uint32_t));
    s->cap = cap;
    return s->slot ? 0 : -1;
}
static void aset_add(aset_t *s, uint32_t a)
{
    uint32_t h;
    if (!a) return;                                   /* NA address */
    h = (a * 2654435761u) & (uint32_t)(s->cap - 1);
    while (s->slot[h]) { if (s->slot[h] == a) return; h = (h + 1) & (uint32_t)(s->cap - 1); }
    s->slot[h] = a;
}
static int aset_has(const aset_t *s, uint32_t a)
{
    uint32_t h;
    if (!a) return 0;
    h = (a * 2654435761u) & (uint32_t)(s->cap - 1);
    while (s->slot[h]) { if (s->slot[h] == a) return 1; h = (h + 1) & (uint32_t)(s->cap - 1); }
    return 0;
}
static void aset_free(aset_t *s) { free(s->slot); }

/* --- deterministic shuffle keyed by seed (self-consistent, not R's RNG) --- */
static uint64_t splitmix(uint64_t *s)
{
    uint64_t z = (*s += 0x9e3779b97f4a7c15ull);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
}
/* perm[] <- a permutation of [0,k) from seed (Fisher-Yates). */
static void make_perm(int32_t *perm, int32_t k, uint64_t seed)
{
    uint64_t s = seed ? seed : 0x9e3779b97f4a7c15ull;
    int32_t i;
    for (i = 0; i < k; i++) perm[i] = i;
    for (i = k - 1; i > 0; i--) {
        int32_t j = (int32_t)(splitmix(&s) % (uint64_t)(i + 1));
        int32_t t = perm[i]; perm[i] = perm[j]; perm[j] = t;
    }
}

/* Collect the addr[] indices of the SNP (rs) probes' M/U addresses, in bead
 * order. Returns count, or -1 on OOM; *out is a malloc'd list the caller frees. */
static int32_t collect_snp_idx(const sesame_idat_t *d, const sesame_index_t *ix,
                               int32_t **out)
{
    int32_t np = sesame_index_nprobes(ix), i, nrs = 0, k = 0;
    const uint32_t *M = sesame__index_M(ix), *U = sesame__index_U(ix);
    aset_t set;
    int32_t *idx;

    for (i = 0; i < np; i++)
        if (strncmp(sesame_index_probe_id(ix, i), "rs", 2) == 0) nrs++;
    if (aset_init(&set, nrs * 2 + 8)) return -1;
    for (i = 0; i < np; i++)
        if (strncmp(sesame_index_probe_id(ix, i), "rs", 2) == 0)
            { aset_add(&set, M[i]); aset_add(&set, U[i]); }

    idx = (int32_t *)malloc((size_t)(nrs * 2 + 8) * sizeof(int32_t));
    if (!idx) { aset_free(&set); return -1; }
    for (i = 0; i < d->n; i++)
        if (aset_has(&set, d->addr[i])) idx[k++] = i;
    aset_free(&set);
    *out = idx;
    return k;
}

/* Rewrite the IDAT with the (modified) mean[] in place of the Mean section,
 * copying the head and tail verbatim. gz input is read transparently; the
 * output is a plain .idat. */
static int rewrite_idat(const char *in_path, const char *out_path,
                        const uint16_t *mean, int32_t n, int64_t off_mean,
                        sesame_err_t *err)
{
    gzFile in = gzopen(in_path, "rb");
    FILE *out;
    char buf[65536];
    int64_t rem;
    int32_t i;

    if (!in) return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", in_path);
    if (!(out = fopen(out_path, "wb")))
        { gzclose(in); return sesame__fail(err, SESAME_ERR_IO, "cannot write %s", out_path); }

    for (rem = off_mean; rem > 0; ) {                 /* head: bytes before Mean */
        int want = rem < (int64_t)sizeof buf ? (int)rem : (int)sizeof buf;
        int got = gzread(in, buf, (unsigned)want);
        if (got <= 0) goto ioerr;
        if (fwrite(buf, 1, (size_t)got, out) != (size_t)got) goto ioerr;
        rem -= got;
    }
    for (i = 0; i < n; i++) {                          /* new Mean section (LE) */
        unsigned char b[2] = { (unsigned char)(mean[i] & 0xff),
                               (unsigned char)((mean[i] >> 8) & 0xff) };
        if (fwrite(b, 1, 2, out) != 2) goto ioerr;
    }
    for (rem = (int64_t)n * 2; rem > 0; ) {            /* skip old Mean in input */
        int want = rem < (int64_t)sizeof buf ? (int)rem : (int)sizeof buf;
        int got = gzread(in, buf, (unsigned)want);
        if (got <= 0) goto ioerr;
        rem -= got;
    }
    for (;;) {                                         /* tail: copy the rest */
        int got = gzread(in, buf, (unsigned)sizeof buf);
        if (got < 0) goto ioerr;
        if (got == 0) break;
        if (fwrite(buf, 1, (size_t)got, out) != (size_t)got) goto ioerr;
    }
    gzclose(in);
    if (fclose(out)) return sesame__fail(err, SESAME_ERR_IO, "write error on %s", out_path);
    return SESAME_OK;

ioerr:
    gzclose(in); fclose(out);
    return sesame__fail(err, SESAME_ERR_IO, "I/O error rewriting %s", out_path);
}

/* Shared core: read, collect SNP indices, let `mode` transform their means,
 * rewrite. mode 0 = zero, 1 = scramble(seed), 2 = restore(seed). */
static int deid_core(const char *in_path, const char *out_path,
                     const sesame_index_t *ix, int mode, uint64_t seed,
                     sesame_err_t *err)
{
    sesame_idat_t *d = NULL;
    int32_t *idx = NULL, k, j, rc = SESAME_ERR_IO;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (sesame_idat_read(in_path, &d, err) != SESAME_OK) return err ? err->code : SESAME_ERR_IO;

    k = collect_snp_idx(d, ix, &idx);
    if (k < 0) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto done; }

    if (mode == 0) {                                   /* zero the SNP means */
        for (j = 0; j < k; j++) d->mean[idx[j]] = 0;
    } else {                                           /* scramble / restore */
        uint16_t *vals = (uint16_t *)malloc((size_t)(k ? k : 1) * sizeof(uint16_t));
        int32_t *perm  = (int32_t  *)malloc((size_t)(k ? k : 1) * sizeof(int32_t));
        if (!vals || !perm) { free(vals); free(perm); rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto done; }
        for (j = 0; j < k; j++) vals[j] = d->mean[idx[j]];   /* current values */
        make_perm(perm, k, seed);
        if (mode == 1) for (j = 0; j < k; j++) d->mean[idx[j]]       = vals[perm[j]]; /* deid[j]=orig[perm[j]] */
        else           for (j = 0; j < k; j++) d->mean[idx[perm[j]]] = vals[j];       /* invert */
        free(vals); free(perm);
    }

    rc = rewrite_idat(in_path, out_path, d->mean, d->n, d->off_mean, err);
    if (rc == SESAME_OK)
        fprintf(stderr, "sesame: %s -> %s (%d SNP bead%s %s)\n", in_path, out_path,
                k, k == 1 ? "" : "s", mode == 0 ? "zeroed" : mode == 1 ? "scrambled" : "restored");
done:
    free(idx);
    sesame_idat_free(d);
    return rc;
}

int sesame_deidentify(const char *in_path, const char *out_path,
                      const sesame_index_t *ix, int randomize, uint64_t seed,
                      sesame_err_t *err)
{ return deid_core(in_path, out_path, ix, randomize ? 1 : 0, seed, err); }

int sesame_reidentify(const char *in_path, const char *out_path,
                      const sesame_index_t *ix, uint64_t seed, sesame_err_t *err)
{ return deid_core(in_path, out_path, ix, 2, seed, err); }
