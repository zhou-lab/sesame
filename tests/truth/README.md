# Frozen ground truth

Expected per-probe values for the tests that compare against the R package.
Committed, so `make test` grades against them **without needing R**, and so a
divergence from R is a row in a table rather than a threshold in a script.

## The rule that keeps this honest

A truth file is **not** a snapshot of what sesame happens to print. Every row
carries where its value came from:

| `src` | meaning |
|---|---|
| `R` | R and sesame agree; the value is R's. The overwhelming majority. |
| `C:<code>` | they differ, or R cannot compute it. The value is sesame's, and `<code>` names a documented reason in `../../DIVERGENCES.md`. |

A `C:` row may only be added together with its entry in `DIVERGENCES.md`, which
says what differs, why, and what evidence settled it. Regenerating must never
silently convert an `R` row into a `C:` row — `make truth` prints every such
flip and refuses to write until each has a reason.

This is the difference between freezing a result and documenting a
disagreement. A frozen result nobody can explain is a bug with a test around
it.

## Regenerating

    make truth RSCRIPT=<an R with the latest sesame/sesameData>

Needs the store, the test IDATs and the oracle. It rewrites the `.tsv.gz`
files and `PROVENANCE.tsv`, which records the R, sesame, sesameData and
annotation versions each value was produced under — without those four numbers
"matches R" means nothing. Review the diff before committing: a changed value
is either a fix or a regression, and the file cannot tell you which.

## Why not just run R in the test?

Three reasons, all met in practice rather than in theory:

1. **R is not always available.** CI has no sesame, so today every numeric gate
   skips there and only the local release run exercises them.
2. **R is not always able.** `detectionPnegEcdf` raises
   `ecdf(negctls$G): 'x' must have 1 or more non-missing values` on EPIC and
   HM450 with sesameData 1.29.10 — no oracle at all for half the platforms.
3. **R moves underneath the test.** The annotation inside `sesameData` and the
   annotation in the store are versioned separately, and when they disagree the
   suite goes red for a reason that has nothing to do with this code.
