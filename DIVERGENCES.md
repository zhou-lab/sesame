# Where sesame's answer differs from the R package's

R is sesame-cli's oracle: almost everything here is validated by running the R
package on the same input and demanding the same number, to the ULP where the
arithmetic allows it (`NUMERICS.md` has the per-step gates). This file is the
short list of places where the two **do not** agree, each with what differs,
why, and what evidence settled it.

It exists because the alternative is worse. A test that quietly excludes the
probes it cannot explain, or a claim of exactness with an asterisk nobody can
follow, hides exactly the thing a user needs when a number looks wrong. Every
code below is referenced from the frozen ground truth in `tests/truth/`, and
`tests/make_truth.sh` refuses to freeze a value whose code is not documented
here.

**None of these is a tolerance.** They are places where the two programs are
answering slightly different questions, usually because they read different
annotation.

## Which side anchors the truth

The frozen ground truth in `tests/truth/` is anchored on **sesame's own output
against the store's annotation**, not on R's. The reason is the annotation
rather than the code: the store carries InfiniumAnnotation v8.1, which is a
newer generation than the copy inside `sesameData`, and where the two disagree
it is the older one that is behind — visibly so in the control-probe names,
where the store's `ctl_10609447_NEGATIVE` carries the control's type and
`sesameData`'s `ctl_10609447` does not.

R remains the oracle, and that is not a contradiction. Every value the two
sides agree on is recorded as `R`, which is the overwhelming majority and is
what validates the port; the disagreements are recorded as `C:<code>` with an
entry below. So R still has to agree with us everywhere it can, and where it
cannot, the reason is written down rather than averaged away.

What this does **not** claim is that sesame is right about any particular
probe in Illumina's terms. It says which annotation generation this project
builds against, which is a statement about our inputs. Settling individual
probes against Illumina's own manifest would bind this file to re-settle them
on every annotation release.

---

### `C:vcf-channel`

**What differs.** On 22 of 127,572 EPICv2 SNP probes, `sesame vcf` and R's
`formatVCF` report complementary variant fractions — the two PVFs sum to
1.0000 — and 11 of those cross a genotype boundary.

**Why.** The out-of-band allele fraction of a Type-I probe is the share of
signal in the *other* colour channel, so it depends on which channel the probe
is assigned to. The two sides read different manifests: sesame reads the
store's `EPICv2.ordering.tsv.gz` (InfiniumAnnotation v8.1), R reads the one
inside `sesameData` (1.29.10 when this was measured). Those two disagree about
the colour channel of exactly these 22 probes — `G` in one, `R` in the other,
all 22 of them. The formula is identical on both sides (`R/sesame.R:242`
versus `src/vcf.c:54`, `pmax(other,1)/pmax(total,2)`); only the partition into
`InfIG`/`InfIR` differs.

**Evidence.** Measured 2026-09-20 on GM12878 EPICv2
(`206909630042_R08C01`): all 22 disagreeing probes have opposite `col` in the
two manifests, and C+R PVF sums to 1.0000 on every one. On the remaining
127,550 probes, genotype and fraction are exact.

**Which manifest is right is deliberately not decided.** sesame's job is to be
correct *given the annotation it was handed*, and deciding whose channel
assignment is correct would bind this file to re-decide on every annotation
release — and would put a claim about Illumina's manifest inside a fidelity
gate. The frozen truth therefore records sesame's value **for the store's
ordering**, which is a statement about our inputs, not about the array.

### `C:vcf-extra-probe`

**What differs.** sesame genotypes a probe R's output does not contain.

**Why.** A `.cg` is positional against the ordering, so every ordering row
exists; R's SigDF drops probes its manifest does not carry. The value is
sesame's because R has none.

### `C:pneg-no-negctl` (resolved 2026-09-22 — kept for the frozen truth)

**What differed.** On EPIC and HM450, R could not compute
`detectionPnegEcdf` at all: it raised *"'x' must have 1 or more non-missing
values"* from `ecdf(negctls$G)`, because the SigDF it built carried no usable
negative controls.

**Why.** Not an annotation difference at all — a bug in the R package.
`readControls()` built the control data.frame by positional assumption rather
than matching the control addresses against the decoded matrix, so on
platforms whose control block is not laid out the way it assumed, every
control signal came back `NA`. Fixed in `R/sesame.R` (sesame 1.31.5); with the
fix R answers on all four platforms and agrees with sesame.

**Status.** No longer a divergence. The code is still recognised by
`tests/merge_truth.py` so that an older oracle can be frozen against, and so
that this entry explains any archived truth row carrying it.

### `C:pneg-control-probe`

**What differs.** sesame reports `p = 1` for all-NA control probes that R's
SigDF omits entirely.

**Why.** Positional storage carries every ordering row, including controls
with no signal. Documented already in `NUMERICS.md`; listed here because the
frozen truth has to source those rows from somewhere.

### `C:pneg-design-type`

**What differs.** On EPIC and HM450, two probes — `cg07162498` and
`cg09334382` — get a detection p-value from sesame that R does not reproduce
(EPIC: 0 versus 0.00243309, i.e. 1/411).

**Why.** The two orderings disagree about these probes' **design type**, not
merely their colour. The store's v8.1 ordering calls `cg07162498` Infinium-I
red; `sesameData` calls it Infinium-II. A Type-I probe is read from two
addresses in one channel, a Type-II probe from one address in both, so the two
sides feed genuinely different intensities into the same ECDF. There is no
common value to compare.

A pure `G`↔`R` channel flip does **not** have this effect here, and is not
excluded: the comparison runs on the raw signal, where the p-value is
`pmin` over `pmax(M,U)` in each colour, so swapping the colours of a Type-I
probe leaves the pair of queries unchanged. MSA has 189 such flips and still
agrees at `max|diff| = 0.00e+00`.

**Evidence.** Measured 2026-09-22 against InfiniumAnnotation v8.1 and
sesameData 1.29.10. Excluding just these two probes takes EPIC from
`max|diff| = 2.43e-03` to exactly `0`. Counts of design-type disagreements per
platform: EPICv2 0, MSA 0, EPIC 2, HM450 2 — pinned in `tests/run_pneg.sh`,
which fails loudly if a manifest moves them.

**A measurement trap worth recording.** `sesameData` stores `col` as a factor
with levels `G`/`R` only, so Infinium-II probes are `NA`. Comparing the two
orderings with `col.R != col.C` therefore returns `NA` for exactly these rows
and `sum(..., na.rm=TRUE)` drops them — which is how an earlier pass counted
"5 channel-disagreeing probes" on EPIC and concluded, wrongly, that excluding
them changed nothing. Normalise `NA` to `"2"` before comparing.

### `C:pneg-negctl-type`

**What differs.** On HM450, about 4,500 probes differ from R by ~1.4e-3 —
every one of them by exactly the same pair of values, 82/614 versus 81/613.

**Why.** The two sides run the ECDF over negative-control pools of different
size. `sesameData`'s HM450 control table has **614 NEGATIVE controls and no
RESTORATION category at all**; the store's v8.1 ordering has **613 NEGATIVE
plus one RESTORATION**, at address 41636384. Both tables hold 850 controls, so
one control has simply been relabelled: `sesameData` lists 41636384 as
"Negative 604" (with the sentinel `Color_Channel = -99`), the store as
`ctl_41636384_RESTORATION`. The restoration control is a real HM450 control in
Illumina's own manifest, so the store's labelling is the one that matches the
array.

**Evidence.** Measured 2026-09-22. The negative-control address sets are
otherwise identical — EPIC 411 versus 411 with no difference either way, HM450
614 versus 613 differing only in 41636384. `tests/run_pneg.sh` drops that one
control from R's pool, after which HM450 agrees at exactly `0` over every
probe outside `C:pneg-design-type`.

**A second trap, in how the control is located.** `readControls()` drops
controls whose signal is NA, so `attr(sdf, "controls")` is shorter than the
`sesameData` control table — 848 rows versus 850 on this chip. Dropping the
control by its position in that table therefore removes a *different* negative
control: the pool is the right size, the addresses are wrong, and the residual
shrinks from 4.9e-2 to 1/613 = 1.6e-3 instead of going to zero. The harness
keys on the row name (`make.names("Negative 604")` = `Negative.604`) and
cross-checks its `G`/`R` against the SigDF row for address 41636384, so a
rename in a future `sesameData` fails the gate rather than passing quietly.

### `C:liftover-ctl-name`

**What differs.** On EPIC, 635 control probes are named `ctl_<address>_<TYPE>`
by sesame (`ctl_10609447_NEGATIVE`) and `ctl_<address>` by R
(`ctl_10609447`). Nothing else differs: same addresses, same values.

**Why.** The two sides take their target ordering from different places —
sesame from the store's `EPIC.ordering.tsv.gz`, R from
`sesameDataGet("EPIC.address")$ordering` (`R/mLiftOver.R`, `convertProbeID`).
The store's newer annotation adds the control's type to its name.

**Evidence.** Measured 2026-09-20: both orderings hold exactly 866,553 probes;
the 635 names that differ are all controls, and with the type suffix stripped
the two sets are identical. Compared as raw strings this reported 1,270 set
differences — 635 each way — for what is one naming change.

`tests/run_liftover.sh` normalises control names before comparing, so the R
cross-check stays exact on all 866,553 probes, and the frozen truth keeps
sesame's spelling, which is the one the store publishes.

