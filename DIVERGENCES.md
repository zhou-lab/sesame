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

### `C:pneg-no-negctl`

**What differs.** On EPIC and HM450, the whole `detectionPnegEcdf` p-value
vector is sesame's, with no R comparison at all.

**Why.** R cannot compute it on those platforms with sesameData 1.29.10:
`detectionPnegEcdf` calls `ecdf(negctls$G)` and raises *"'x' must have 1 or
more non-missing values"* because the SigDF it builds carries no usable
negative controls. This is not a disagreement about a number — the oracle does
not produce one.

**Evidence.** 2026-09-20: EPICv2 and MSA agree with R at `max|diff| = 0.00e+00`
over 284,309 and 937,690 probes; EPIC and HM450 raise the error above. The
frozen values keep those two platforms under a regression gate that would
otherwise not exist, and they will be re-compared to R the moment R can answer.

### `C:pneg-control-probe`

**What differs.** sesame reports `p = 1` for all-NA control probes that R's
SigDF omits entirely.

**Why.** Positional storage carries every ordering row, including controls
with no signal. Documented already in `NUMERICS.md`; listed here because the
frozen truth has to source those rows from somewhere.

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
