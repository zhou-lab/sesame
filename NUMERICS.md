# NUMERICS.md — divergence register

## What these numbers were measured against

Every claim below is a measurement, and a measurement without its oracle is an
assertion. As of **2026-09-20** the numbers here come from:

| | |
|---|---|
| R | 4.6.0 (2026-04-24), `Rscript-4.6.0` from the lab's shared_Renv shims |
| sesame | 1.29.5, re-measured identically under 1.31.4 |
| sesameData | 1.29.10 |
| Bioconductor | 3.23.1 |
| annotation | the store at InfiniumAnnotation v8.1 |
| sesame-cli | 0.5, YAME v1.51 |

This block was missing until 2026-09-20, and its absence cost real time: when
`test-vcf` began failing, nothing in this file said which sesame the "exact"
claim had been measured against, so whether 22 disagreeing probes were new or
long-standing had to be re-derived from scratch rather than read off. The same
four versions are captured per run in `tests/truth/PROVENANCE.tsv`, and the
release SOP requires them in the log entry.

The oracle was upgraded to sesame **1.31.4** on 2026-09-22 and every number
above was re-measured: **all identical**, including the 22 excluded vcf probes
and the liftover set-difference. That is the expected result and a useful
control — the divergences live in `sesameData`, which is unchanged at 1.29.10,
not in sesame's code. They will only move when `sesameData` does, and 1.29.x
is the devel series for Bioc 3.23, so that is a next-cycle event.


sesame is an independent implementation of sesame's basic preprocessing. It is
**not** bit-exact with the R package by mandate: where the R code has a defect or
a numerically fragile construction, sesame fixes it and records the difference
here.

The R package (`zwdzwd/sesame`) remains the reference oracle and is never
modified by this project. Every difference below must be *justified*, not merely
reported — for numerical reformulations that means a proof or an
arbitrary-precision oracle showing sesame is closer to truth.

| id | site | R behavior | sesame behavior | status |
|----|------|-----------|------------------|--------|
| D1 | `R/sesame.R:293-298` `readIDAT1` | Merges Grn/Red by **positional `cbind`**, rownames taken from Grn. No join on IlluminaID. The assumption *holds for well-formed files* (verified: Grn/Red address vectors are identical, sorted and unique on every test array). It fails silently and array-wide only when mismatched IDATs are paired. | Verifies the address vectors match and errors otherwise. Confirmed: pairing an EPICv2 Grn with an HM450 Red is rejected, where R would silently emit garbage. | **done (P1)** |
| D2 | `R/detection.R:163-165` `pOOBAH` | `pmax(..., na.rm=TRUE)` but `pmin(...)` **without** `na.rm` → an NA in one channel yields `p=NA`, and `NA > 0.05` is NA, so the probe is **not** masked. Sibling `ELBAR` at `R/detection.R:60` does `pvals[is.na(pvals)] <- 1.0`, evidence this is an oversight. | Both channels NA → `p=1` (mask). One channel NA → use the other. | **done (P)** |
| D3 | `R/readIDAT.R:22` vs `:178` | Header parsed **twice**: `readIDAT()` reads `version` as int32 to dispatch, then `readIDAT_nonenc()` reopens the file and re-reads it as int64. Works only because the value is small and little-endian. | Parsed once, as int64. | **done (P0)** |
| D4 | `R/readIDAT.R:44-46` `readLong` | `readBin(what="integer", size=8)` — base R does not support 8-byte integers; this works by accident. | Real `int64_t`. | **done (P0)** |
| D5 | `R/background.R:143-167` `normExpSignal` | Computes `exp(dnorm(...,log=TRUE) - pnorm(...,lower.tail=FALSE,log.p=TRUE))`; a difference of logs, which loses precision in the left tail where the `1e-6` floor fires. | Evaluates the inverse Mills ratio `λ(t)=φ(t)/Φ(t)` directly (`erfc` for `t > -5`, Laplace continued fraction below). Mathematically identical, better conditioned, never differences logs. | **done (B)** |
| D6 | `MASS::huber` via `R/background.R:132,135` | **Errors** when `MAD == 0`. | Falls back to `s = max(IQR/1.349, 1e-8)` and sets `SESAME_STAT_NOOB_MAD0` — a single degenerate array must not abort a 100k-sample run. | **done (B)** |
| D7 | `R/dye_bias.R:122`, `R/background.R:98`, `R/detection.R:155`, `R/mask.R:173` | Escape hatches are **silent**: `dyeBiasNL`→`maskIG` when RGdistort is NA or >10; `noob` returns unchanged if background <100; `pOOBAH` substitutes a `1:1000` prior; `backgroundMask` returns NULL for unknown platforms (silently changing both `noob` and `pOOBAH` background composition). | Each sets a bit in a per-sample status word exposed by the API, so pipelines can count fallbacks instead of discovering them in the betas. Done for `dyeBiasNL` (`SESAME_STAT_DYEBIAS_FAILED`) and `noob` (`SESAME_STAT_NOOB_SKIPPED`); `pOOBAH` prior + unknown-platform mask still planned. | partial |
| D8 | `preprocessCore::normalize.quantiles.use.target` via `R/dye_bias.R:135,150` | preprocessCore's own compiled C (`qnorm.c`). | Clean-room reimplementation — see below. Algebraically identical; **not bit-identical**. | **done (D)** |

### D8 — the one difference I could not eliminate, and why

`normalize.quantiles.use.target` is the only compiled third-party dependency in
the QCDPB path. preprocessCore is `LGPL (>= 2)`; sesame-cli is AGPL-3.0-or-later.

**Since sesame-cli went AGPL (2026-07-17), vendoring `qnorm.c` is now legally
open:** LGPL is upward-compatible with GPL/AGPL (a recipient may elect LGPL-3,
which is GPL-3 plus permissions and thus absorbable into an AGPL-3 work), and
AGPL already requires shipping complete source, so LGPL's relink obligation is
moot. So this entry could be *erased* by vendoring preprocessCore's `qnorm.c`
for bit-exactness.

**It is deliberately kept clean-room anyway**, for engineering reasons rather
than licensing ones:

- the clean-room already agrees with R to ~2 ULP (measured below) — biologically
  and numerically meaningless;
- vendoring adds an LGPL component and pulls in preprocessCore's `include/` +
  stubs, which are `LinkingTo:` plumbing for R packages, not a clean library a
  standalone C program links against — i.e. real extraction work;
- the clean-room is ~40 lines with no dependency.

If exact reproduction of R's last bit ever matters more than that, vendoring is
now the sanctioned route. Until then the clean-room stands. It was characterized
purely by black-box probing:

| probe | observed |
|---|---|
| `n == m`, no ties | `sorted(target)[rank(x)]` |
| `n != m` | `quantile(target, (rank-1)/(n-1), type=7)` |
| ties in `x` | quantile at the **average** rank |
| NA in target | dropped before ranking |
| NA in `x` | preserved; ranking uses non-NA only |

The tie rule was pinned with a deliberately non-linear target, which separates
"quantile at the average rank" (actual: `100`) from "mean of the quantiles at
each tied rank" (`106.67`). The former wins.

**The algebra is settled; the last bit is not.** I could not find an arithmetic
form that reproduces `qnorm.c` bit-for-bit. Both `(1-h)*lo + h*hi` and
`lo + h*(hi-lo)`, across three index formulations, top out at ~92% bitwise
agreement — so preprocessCore evaluates it in some third order I would have to
read the source to learn, which is precisely what the licence forbids.

Measured consequence on end-to-end betas (`prep="CD"`, NA patterns identical
everywhere):

| platform | betas differing | max abs | max rel |
|---|---|---|---|
| HM450 | 1.05% | 5.6e-16 | 2.5e-15 |
| EPIC | 1.47% | 1.0e-15 | 4.0e-15 |
| EPICv2 | 1.82% | 7.8e-16 | 3.0e-15 |
| MSA | 4.39% | 3.3e-16 | 1.1e-15 |

~2 ULP on a value in [0,1]. The acceptance gate for `D` is max relative
difference ≤ 1e-12, which this clears by three orders of magnitude — and the
project's own stated gate was 1e-8. This is accepted as last-bit arithmetic
ordering, not an algorithmic difference: it is not attributable to any decision
in this implementation, and neither result is "more correct" than the other.

## Mask lineage — Q is not bit-identical to R, by data version

The `Q` step (qualityMask) reads the recommended mask sets from the platform's
YAME `.cm` in the store. sesame-cli links YAME directly (both are AGPL-3.0, same
author) and reads the `.cm` in-process via the YAME C API -- no `yame` binary at
runtime. A probe is masked iff it is set in any of `recommendedMaskNames(platform)` — `M_1baseSwitchSNPcommon_5pt`,
`M_2extBase_SNPcommon_5pt`, `M_mapping`, `M_nonuniq`, `M_SNPcommon_5pt` for MSA.

sesame-cli's `Q` is verified **self-consistent**: it masks exactly the yame
union applied to the ordering (14,494 probes on the MSA test array — an exact
gate, `tests/run_qmask.sh`).

It is **not** bit-identical to R's `qualityMask`, and cannot be, because the
published `.cm` is a *newer mask lineage* than the KYCG object sesameData
currently ships — the same situation as the ordering table (284,309 vs 284,317
probes). Measured against R (`KYCG.MSA.Mask.20260122`) on the MSA test array:

| | count |
|---|---|
| R qualityMask | 14,249 |
| sesame-cli Q (repo `.cm`, "B1 coherent") | 14,494 |
| shared | 14,240 |
| **Jaccard** | **0.982** |

The residual (9 R-only, incl. the 8 deleted probes; 254 C-only) is the
mask-version bump, not an implementation error — per-track the two agree ~90%+
(`M_mapping`: 1693 shared of 1870/1939). This converges when sesameData ships
the same mask version. sesame-cli always applies whatever `.cm` the pinned tag
publishes; the number above pins *which* lineage was measured.

## P (pOOBAH) — algorithm exact, residual is boundary + lineage

`P` reads the background mask (union of `backgroundMask` names) from the same
`.cm`, builds the per-channel ecdf of out-of-band + negative-control signal, and
masks probes whose detection p exceeds 0.05. The ecdf p is `1 - k/n` with
`k = #{bg ≤ x}`, matching R's `1 - ecdf(bg)(x)` arithmetic bit-for-bit.

Verified on the MSA test array: R masks 9,433, sesame-cli 9,445, **R-only 0**
(sesame is a superset here), 12 disagreements. Every disagreement is
*explainable* and the test (`tests/run_poobah.sh`) fails if any is not:

- **all 12 are boundary flips** — R's p is `0.04999…`, i.e. within 1e-3 of the
  0.05 cutoff, so the lineage-different background pool (`.cm` `M_nonuniq` etc.
  are the "B1" build, ordering is 284,309 not 284,317) shifts them across;
- **0 were D2 cases** here, but D2 is implemented: a probe whose channel is
  all-NA gets `p = NA` in R (silently unmasked) and `p = 1` in sesame (masked),
  matching the sibling `ELBAR`.

No disagreement had R p far from 0.05 — i.e. no algorithm error. The p-value
formula is exact; the residual is the same data lineage as Q/D.

## B (noob) — arithmetic proven in isolation; residual is lineage only

`B` parameterizes a Normal background and Exponential foreground per channel
(from the same non-`backgroundMask` out-of-band + negative-control pools as `P`)
and deconvolves every signal via the inverse Mills ratio (D5), with the fixed-
scale Huber location (`MASS::huber`) for `mu`/`sigma`/`alpha` and the D6 fallback
for a zero MAD. Because the fit pools are lineage-different, `B` cannot be
bit-identical to R — so the test (`tests/run_noob.sh`) separates the two things
that could move a beta:

**(1) the arithmetic, isolated.** `normExpSignal` and `MASS::huber` are compared
to the C primitives on *identical* inputs via the `normexp_test` harness:

- `normExpSignal` over a 480-point grid: max **absolute** error `3.8e-9`, bulk
  (p90) **relative** error `1.8e-11`. The single worst relative point is a deep-
  tail near-cancellation (`signal = mu.sf + sigma·λ` with `sigma·λ ≈ -mu.sf`),
  i.e. exactly where R's difference-of-logs loses precision and the Mills form
  does not — there **C is the more accurate of the two**, so R is not a valid
  oracle for that bit. The gate is therefore absolute (`< 1e-6`).
- `MASS::huber` on four vectors incl. a 568k real signal column: `|Δmu| ≤ 1.1e-8`
  (long-double accumulation matches R's `sum()`), `|Δs| ≤ 1.1e-13`. Note MASS
  returns the *pre-update* `mu` on convergence (the break precedes the
  assignment) — replicated exactly. D6 confirmed: a zero-MAD vector, which R
  **errors** on, yields `s = 1e-8` and the status bit in C.

**(2) the integration, on raw-identical probes.** The published ordering assigns
different M/U addresses than sesameData's manifest to **182** replicate probes
(e.g. `cg…_TO12`, `cg…_BC11`), so *their raw betas already differ by up to 0.48
before any prep* — a pre-`B` lineage artifact. Restricting to the 284,127 probes
whose raw betas agree (identical input to noob), betas after `B` match R with
**median `6.5e-6`, p99 `1.0e-3`, max `1.6e-3`**, NA agreement exact. That residual
is the lineage-level fit shift (`muG` 1592.75 vs 1592.72, `sigma` differing by
~1.5) propagating through `normExp`; the worst probe is invariably a negative
control (lowest signal, most sensitive). An earlier real bug — applying the wrong
channel's fit — moved these same betas by 0.4–0.7, so the gate cleanly separates
lineage from error.

## QC panel — matches sesameQC to lineage scale; two noted behaviors

`sesame preprocess --output qc` computes the `sesameQC_calcStats` panel
(detection, numProbes, intensity, channel, dyeBias, betas) from a raw SigDF. It is not bit-identical to
R for the same reason `P`/`B` are not — the published ordering/mask is a newer
lineage than sesameData (fewer probes, different background mask), so probe
counts, detection, and beta stats all shift by a handful of probes. The test
(`tests/run_qc.sh`) pins every metric to within that scale: on MSA all 65 metrics
agree, worst gated relative error `1.75e-2` on the small `G->R` channel-switch
category (an absolute difference of 4 probes); intensity means, probe counts, and
detection fractions match to `< 1e-3`. A formula bug would be far off, not off by
a dozen probes.

Two behaviors are faithful to R but worth stating:

- **`num_dtna` and the D2 fix.** `sesameQC_calcStats_detection` derives
  "N. Probes w/ Missing Raw Intensity" from `sum(is.na(pOOBAH(sdf, return.pval)))`.
  Our own D2 fix to R (`R/detection.R`) makes pOOBAH set those `NA`s to `1`, so in
  the current R package `num_dtna` is **0**. sesame-cli instead reports the metric
  as intended — probes with no usable signal in *either* channel — computed
  directly, independent of pOOBAH's NA convention. The QC test reports these two
  fields (`num_dtna`, `frac_dtna`) but does not gate them.

- **The beta group runs `D->B->P`.** `sesameQC_calcStats_betas` computes
  `getBetas(pOOBAH(noob(dyeBiasNL(sdf))))` — i.e. the beta-distribution stats are
  taken after dye-bias, noob, and detection masking, not on the raw signal like
  the other groups. sesame-cli reproduces that per-sample on a copy, so `num_na`
  equals the `P`-step mask count (9433 R / 9445 C on MSA — the same lineage split
  the `P` test reports).

## DML — the one step with no lineage caveat

Unlike everything above, differential methylation consumes a **betas matrix**, not
IDATs, so there is no ordering/mask lineage in play: the same matrix and metadata
give the same result as R. `sesame dml` fits each probe's betas across samples by
Householder QR least squares (R's `lm` uses LINPACK QR), takes per-coefficient
t-tests and a holdout F-test per categorical variable, and the p-values come from
the regularized incomplete beta `I_x(a,b)` (Lentz continued fraction) via the t
and F tail CDFs.

The test (`tests/run_dml.sh`) generates a small dataset and fits it both ways with
`~ group + age` (including a few NA betas to exercise per-probe dropping). It
matches R's `DML` / `summaryExtractTest` to **~5e-10** on estimates, t and F
p-values, effect sizes, and the BH adjustment — i.e. exact up to QR-vs-LINPACK
rounding. `sesame__betai`/`pt`/`pf` agree with R's `pbeta`/`pt`/`pf` to ~1e-14
across the usual range.

The design is the deliberately-simple part: a main-effects formula
(categorical auto-dummied with treatment contrasts and alphabetical levels, to
match `model.matrix`; continuous as-is) or an explicit numeric design matrix
(`--design`) for anything R-formula-complex. Region calling (DMR) is separate and
needs the per-probe genomic-coordinate annotation, now published as
`<platform>.hg38.coord.tsv.gz`.

## CNV — arithmetic proven against R's lm; residual is channel-lineage leverage

`sesame cnv` reproduces `cnSegmentation` end to end — the bin log2 ratios **and**
their circular binary segmentation (see the CBS note below). For a target sample
it regresses per-probe total intensity on a panel of normals' totals (the **same
Householder QR as DML**,
`sesame__ols`), then `log2(target / max(fitted, 1))` per probe, and bins along the
genome (tile at 50 kb, subtract assembly gaps, left/right-merge to ≥ 20 probes,
median per bin — `getBinCoordinates` + `leftRightMerge1` + `binSignals`).

**The fit + log2 is numerically exact.** Fed R and C the *identical* per-probe
totals (extracted from C's own `.cg` with `attach-probe`, so no data-lineage or
channel ambiguity enters), C's OLS+log2 matches R's `lm` to **max 1.0e-5, median
5e-8** over 937k EPICv2 probes (`tests/run_cnv.sh`; the 1e-5 tail is float32 `.cg`
storage of the totals, not the fit). An independent pure-Python normal-equations
solve on the same inputs agrees with C to ~3e-8.

**Segmentation (CBS).** The bin log2 ratios are segmented by a clean-room
deterministic port of DNAcopy's circular binary segmentation (`src/cbs.c`), the
routine `cnSegmentation` calls. DNAcopy assesses each split by **permutation
against R's RNG**, which `cnSegmentation` never seeds — so its borderline splits
are not reproducible even run-to-run. We keep the deterministic parts and replace
the RNG with analysis:
- the change-point **location** is the arc maximising the F-standardised CBS
  statistic `BSS/((tss-BSS)/(n-2))`, `BSS = n(S_j-S_i)^2/(L(n-L))` — exact
  (a naive O(m²) scan finds the same argmax as DNAcopy's block-decomposed
  `tmaxo`, since real data has no ties);
- **significance** is the Siegmund (1988) analytic tail probability `tailp` that
  DNAcopy itself uses for large segments (`p.method="hybrid"`), applied over the
  full valid arc range so it also covers the small arcs DNAcopy would permute;
- the 1-vs-2 change-point edge test uses DNAcopy's `tpermp` two-sample statistic
  scored against the t distribution instead of a permutation;
- `sdundo` (merge adjacent segments whose median differs by < `undo.SD·trimmed.SD`)
  and `trimmed.variance` (with the compiled-in `inflfact(0.025)`) are already
  deterministic in DNAcopy and ported verbatim.

On real K562 copy-number bins (`tests/run_cbs.sh`, oracle = `DNAcopy::segment` with
sesame's parameters and a fixed seed) this recovers **97% of DNAcopy's segment
breakpoints within one bin**, with exact locations. It calls ~30% more segments:
those extra breakpoints have a median |Δmean| of 0.23 log2 — real local amplitude
changes DNAcopy's conservative permutation rejects, not noise. A clean step signal
segments exactly. Since DNAcopy is itself non-deterministic here, breakpoint recall
(not an exact match) is the meaningful gate.

**Residual vs R's `cnSegmentation` on real IDATs is a data-lineage leverage
effect, not an arithmetic one.** Running the full R path (`readIDATpair` target,
`EPICv2.8.SigDF` normals, its own `probeCoords`) diverges by ~2.6e-4 median.
Traced to the root: ~25 EPICv2 Infinium-I probes carry a **different in-band
channel (`col`) assignment** between the published ordering and the manifest
`readIDATpair` uses (e.g. `cg19844224_TC11`: green in-band gives total 4731, red
gives 53548). These are not a channel-inference difference — running C's `C` step
makes it *worse* (643 probes), confirming `readIDATpair` is not inferring here.
Because the CNV fit is **global**, those ~25 high-intensity probes act as
regression leverage and tilt every bin's log2 by ~1e-4. Same class as the Q/P/B
ordering-lineage divergence: pick one lineage for target, normals, and coords
(all positional in the same ordering, as the store guarantees) and the profile is
internally exact. Match the CNV inputs' processing state to the normal reference's
(the store `cnvnormals.cg` is raw signal, so the target should be raw too).

## VCF (formatVCF) — genotype and fraction exact; only the deep-tail score drifts

`sesame vcf` ports `formatVCF` (`R/vcf.R`): from the RAW SigDF (channel inference
would erase the Type-I channel-switch genotype signal) it forms a per-probe
variant fraction — `1-beta` if the U bead reports ALT, the out-of-band Type-I
allele fraction (`getAFTypeIbySumAlleles`) for `REF_InfI`, else `beta` — then a
3-genotype binomial model (background 0.1, 40 beads) picks 0/0, 0/1, or 1/1 and a
phred-like score `GS = floor(-log10(1 - GL_max/ΣGL)·10)`.

Validated against R on **127,550** of the 127,572 EPICv2 SNP probes (GM12878,
raw; `tests/run_vcf.sh`): the **genotype call (GT) and variant fraction (PVF)
are exact** — 0 GT mismatches, PVF max diff 0.

**The 22 excluded probes are a manifest disagreement, not a divergence**
(measured 2026-09-20). The two sides read different manifests: the C side the
store's `EPICv2.ordering.tsv.gz` at InfiniumAnnotation v8.1, R the one inside
`sesameData` (1.29.10 here). On 22 Type-I probes the two disagree about the
probe's **colour channel** — `G` in one, `R` in the other, all 22. The
out-of-band allele fraction is then computed from the opposite channel, so
`getAFTypeIbySumAlleles` returns the exact complement: C and R PVF sum to
1.0000 on every one of the 22, and 11 of them cross a genotype boundary. The
formula is identical on both sides (`R/sesame.R:242` vs `src/vcf.c:54`,
`pmax(other,1)/pmax(total,2)`); only the partition into `InfIG`/`InfIR`
differs. `tests/run_vcf.sh` therefore grades only probes where the two
manifests agree on the channel, decides that from the manifests alone (never
from whether the answer matched), and **fails if the excluded count moves off
22** — so annotation drift surfaces as a changed number rather than a widening
blind spot. Which manifest is right is **deliberately not decided**: this gate asks whether
the C port reproduces `formatVCF`, and that question does not need an answer to
whose channel assignment is correct. Deciding it would bind the test to
re-decide on every annotation release. The binomial coefficient cancels in the
posterior `GL_max/ΣGL`, so the ratio is computed in log space without it. `GS`
matches on all but **11 probes** (≤10 apart), all at extreme scores (> ~80, beyond
the nominal 7–85 range) where `1 - GL_max/ΣGL` is ~1e-40: there R's `dbinom`
(Loader's saddle-point algorithm) is more accurate in the deep tail than a naive
log-pmf, so `-log10` of a ~1e-40 quantity lands a few integers apart. A
quality-score-only effect on 0.009% of probes, never on a call.

`--variants` drops the non-switching Infinium-I probes (`REF_InfI` with no
overlapping `rs`) — the ~99% bulk that carry no known variant and essentially
always call 0/0 — keeping the rs and channel-switching probes (883 of 127,572 on
EPICv2). The default emits the full set, matching R.

## 0.2.0 ports — fidelity and the one divergence

Each was validated against R with a differential test (`tests/run_*.sh`):

| feature | interface | R oracle | result |
|---|---|---|---|
| `dyeBiasL` | prep code `E` | `dyeBiasL(inferInfiniumIChannel(sdf))` | bit-identical (≤1 ULP), EPICv2+MSA. Only Inf-II betas move (per-channel scaling leaves Inf-I `M/(M+U)` invariant). |
| `detectionPnegEcdf` | `preprocess --detection pneg` | `detectionPnegEcdf(sdf, return.pval=TRUE)` | ≤1 ULP on every real probe. C additionally emits `p=1` for all-NA control probes R's SigDF omits. |
| `betasCollapseToPfx` | `preprocess --collapse` | `betasCollapseToPfx(betas)` | exact (prefix set, NA pattern, values) — validated on C's own betas so the `.cg` float32 rounding is the only residual. |
| `mLiftOver` | `sesame mliftover` | `mLiftOver(betas, target, source)` | exact both directions (EPICv2↔EPIC): identical target set, order, NA pattern, and first-match choice. |
| `imputeBetasMatrixByMean` | `impute --method mean` | `imputeBetasMatrixByMean(mx, axis)` | exact (axis=probe 0.0; axis=sample float32). |
| `bisConversionControl` (GCT) | `qc.tsv` `GCT` column | the GCT formula on the store's ext lists | arithmetic exact (~1e-10) on all platforms; **not** bit-vs-R's `probeInfo` — see note. |

**imputeBetasByGenomicNeighbors** (`impute --method neighbors`) is the one with a
noted divergence. The algorithm matches R exactly — strand-aware, start-anchored
`resize` + **strand-specific** `findOverlaps` + `slice_min` with ties — verified
against a reference replication of R's algorithm on identical coordinates
(`tests/run_neighbors.sh`, 50,024 mapped EPICv2 probes, ≤ float32). Two points:

- **Coordinate source.** R reads `sesameData_getManifestGRanges`; the CLI reads the
  store's `<plat>.<genome>.coord.tsv.gz`. These are the same mapping for **every
  mapQ≥1 probe** and disagree only on **mapQ=0** (multi-/non-mapping) probes, which
  impute to noise regardless. So mapped, confidently-placed probes match R; the
  ~4% of imputed probes that don't all trace to a mapQ=0 probe in the probe or its
  neighbourhood.
- **Unmapped probes are left NA.** R places all unmapped probes at a shared
  position 0 so they "neighbour" each other and receive the mean of all unmapped
  betas — a meaningless artifact the CLI declines to reproduce.

**GCT extension lists (`typeI_ext`) — a curation-lineage difference, like the Q
mask.** The GCT ratio is `mean(green of extension-C Type-I probes) / mean(green of
extension-A/T)`. R's `bisConversionControl` reads the extension classes from
sesameData's `probeInfo$typeI.extC/extT` — a curated subset available only for
EPIC/HM450. sesame-cli reads the annotation pipeline's `<plat>.typeI_ext.tsv.gz`,
generated **uniformly for every platform** and deliberately slightly broader (on
EPIC it classifies ~3,463 more Type-I extension probes than `probeInfo`, almost
all additive). The **arithmetic is exact** — validated to ~1e-10 against the GCT
formula recomputed in R on the *same* lists (`tests/run_gct.sh`, all four
platforms) — but because the input set differs, GCT is **not** bit-identical to
`bisConversionControl` (≈7.6e-3 on EPIC/HM450). This is the same class of
data-version/curation difference as the Q mask above, not an arithmetic one, and
it is what lets GCT be defined for EPICv2/MSA at all (where R has no built-in
lists).

## Observations about the R implementation (not divergences)

**Ordering files are sorted by Probe_ID under R's *locale* collation, not byte
order.** All four platforms are in `order()` sequence under `en_US.UTF-8`.
HM450/EPIC/EPICv2 happen to coincide with byte order; **MSA does not**. It is the
only platform carrying `rs` IDs where one number prefixes another, e.g.

```
locale : rs12051_TC21 < rs12051548_TC21     <- what the files use
radix  : rs12051548_TC21 < rs12051_TC21
```

ICU ignores the `_` and compares `rs12051TC21` vs `rs12051548TC21` (`5` < `T`);
byte order compares `_` (0x5F) vs `5` (0x35) and flips it. Exactly 11 such pairs
in MSA.

This does not affect sesame-cli: Probe_IDs are consumed in file order and the
only binary search is over integer addresses. Two hazards it *does* create:

1. **Rebuild reproducibility.** `order()` is locale-dependent, so regenerating
   the `.address` object under `LC_COLLATE=C` (Docker, CI, minimal images) would
   reorder MSA, change its sha256, and change the output row order — which is
   load-bearing (`R/sesame.R:504`). `tools/export_ordering.R` preserves the
   object's existing order and is therefore safe; anyone rebuilding *upstream*
   must pin the locale.
2. **Never binary-search Probe_IDs with `strcmp`.** It would silently fail on
   those 11 MSA probes. Use a hash if ID lookup is ever needed.

**The ordering table's `mask` column is dead data.** `.address$ordering` carries
a `mask` column built by `sesameAnno_buildAddressFile` from
`create_default_mask()$ref_issue` (`R/sesameAnno.R:170`). Nothing ever reads it:
`chipAddressToSignal` explicitly sets `mask=FALSE` (`R/sesame.R:475,498`) and a
grep over `R/` finds no consumer anywhere in the pipeline. sesame parses the
column but deliberately ignores it, matching R. Seeding the SigDF mask from it
would wrongly mask ~257 (HM450) / ~2263 (EPICv2) / ~3120 (MSA) probes that R
keeps. Design masking is `qualityMask`'s job (the `Q` step), which uses a
different, curated vocabulary. *This is worth a look from the sesame side: the
column is either vestigial or an intent that was never wired up.*

> The published InfiniumAnnotation ordering (tag v7) has since **dropped the
> `mask` column** — it is now a 4-column table (`Probe_ID M U col`), with the
> curated masks living in the companion `.cm` file. `sesame_index_open` accepts
> either shape (4 columns → mask defaults 0), so nothing downstream changes; the
> point above stands, the dead column is simply gone from the source now.

**R's text parser is not correctly rounded — do not compare betas via text.**
`as.numeric("0.96236179722418258")` returns a double one ULP away from the value
that produced that string, even though C's own `strtod` round-trips it exactly.
R uses its own `R_strtod` rather than a correctly-rounded parser. Comparing
through `read.table` manufactured a ~2.2e-16 "divergence" across ~19% of probes
that did not exist in the computation: with raw float64 the betas are
**bit-identical**. The golden harness therefore compares binary (`--f64`), never
text.

## Deliberately *not* divergent

- **Row order.** The SigDF row order follows the manifest ordering table
  (`R/sesame.R:504`). This is load-bearing and is preserved exactly.
- **`approx(ties=mean)`** in `dyeBiasNL` (`R/dye_bias.R:143,158`). The tie
  collapse is load-bearing — intensities are integer-valued with massive
  duplication and tied `x` map to different `y`. Replicated exactly.
- **`quantile` type 7** at `R/channel_inference.R:32` and inside `IQR` at
  `R/background.R:103-104`. Replicated exactly.
- **`getBetas` clamps** `pmax(M,1)/pmax(M+U,2)` — note sesame uses **no** `+100`
  Illumina offset. Replicated exactly.
