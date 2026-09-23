/* fuzz_idat.c -- fuzz harness for the IDAT reader.
 *
 * The reader consumes untrusted files (IDATs are routinely pulled from GEO)
 * and the header is fully attacker-controlled: nFields, per-field byteOffset,
 * and nSNPsRead all steer allocation and seeking.
 *
 * Build (Linux/clang with libFuzzer):
 *     clang -std=c11 -g -O1 -Iinclude -fsanitize=fuzzer,address,undefined \
 *           fuzz/fuzz_idat.c src/idat.c -lz -o fuzz_idat
 *     ./fuzz_idat corpus/
 *
 * Build (standalone replayer -- works with Apple clang, which ships no
 * libFuzzer; run a corpus under ASan/UBSan):
 *     make fuzz-replay && ./fuzz_replay corpus/all-inputs
 *
 * sesame_idat_read() takes a path, so the harness materialises the input into
 * a temp file. Slower than an in-memory reader, but it exercises the real
 * zlib/seek path rather than a parallel one.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static void run_once(const uint8_t *data, size_t size)
{
    char tmpl[] = "/tmp/sesame_fuzz_XXXXXX";
    int fd = mkstemp(tmpl);
    if (fd < 0) return;

    if (size && write(fd, data, size) < 0) { close(fd); unlink(tmpl); return; }
    close(fd);

    sesame_idat_t *d = NULL;
    sesame_err_t e;
    if (sesame_idat_read(tmpl, &d, &e) == SESAME_OK) {
        /* Touch every buffer so ASan sees any bad bound. */
        volatile uint64_t acc = 0;
        for (int32_t i = 0; i < d->n; i++)
            acc += d->addr[i] + d->mean[i] + d->sd[i] + d->nbeads[i];
        (void)acc;
        sesame_idat_free(d);
    } else {
        /* On failure the reader must not have handed us anything to free. */
        if (d != NULL) abort();
    }
    unlink(tmpl);
}

#ifdef SESAME_FUZZ_STANDALONE
int main(int argc, char **argv)
{
    for (int i = 1; i < argc; i++) {
        FILE *f = fopen(argv[i], "rb");
        if (!f) continue;
        fseek(f, 0, SEEK_END);
        long n = ftell(f);
        fseek(f, 0, SEEK_SET);
        if (n < 0) { fclose(f); continue; }
        uint8_t *buf = (uint8_t *)malloc((size_t)n ? (size_t)n : 1);
        if (!buf) { fclose(f); continue; }
        size_t got = fread(buf, 1, (size_t)n, f);
        fclose(f);
        run_once(buf, got);
        free(buf);
        fprintf(stderr, "replayed %s (%zu bytes)\n", argv[i], got);
    }
    return 0;
}
#else
int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    run_once(data, size);
    return 0;
}
#endif
