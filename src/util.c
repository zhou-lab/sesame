/* util.c -- shared helpers.
 * SPDX-License-Identifier: LicenseRef-CHOP-Academic-BSD-2-Clause
 *
 * Copyright (C) 2026-present The Children's Hospital of Philadelphia
 *
 * Use of this software is available to academic and non-profit institutions
 * for research purposes under the 2-Clause BSD License; for use or transfers
 * to commercial entities, inquire with Dr. Wanding Zhou at zhouw3@chop.edu.
 * See the LICENSE file at the root of the repository for the full terms.
 */
#include "internal.h"

#include <stdarg.h>
#include <stdio.h>

int sesame__fail(sesame_err_t *err, int code, const char *fmt, ...)
{
    if (err) {
        va_list ap;
        err->code = code;
        va_start(ap, fmt);
        vsnprintf(err->msg, sizeof(err->msg), fmt, ap);
        va_end(ap);
    }
    return code;
}

const char *sesame_strerror(int code)
{
    switch (code) {
    case SESAME_OK:              return "ok";
    case SESAME_ERR_IO:          return "I/O error";
    case SESAME_ERR_FORMAT:      return "malformed input";
    case SESAME_ERR_UNSUPPORTED: return "unsupported input";
    case SESAME_ERR_NOMEM:       return "out of memory";
    default:                     return "unknown error";
    }
}
