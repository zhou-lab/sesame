# Sesame2

[![build](https://github.com/zwdzwd/sesame-cli/actions/workflows/conda-build.yml/badge.svg)](https://github.com/zwdzwd/sesame-cli/actions/workflows/conda-build.yml)
[![conda](https://img.shields.io/conda/vn/zhou-lab/sesame?label=conda)](https://anaconda.org/zhou-lab/sesame)
[![license](https://img.shields.io/badge/license-BSD--2--Clause%20(academic)%20%2F%20commercial-blue.svg)](LICENSE)
[![coverage](https://img.shields.io/endpoint?url=https%3A%2F%2Fzwdzwd.github.io%2Fsesame-cli%2Fcoverage.json)](scripts/coverage.sh)
[![docs](https://img.shields.io/badge/docs-online-blueviolet)](https://zwdzwd.github.io/sesame-cli/)
[![language](https://img.shields.io/badge/C-C11-00599C)](include/sesame.h)
[![arrays](https://img.shields.io/badge/arrays-EPIC%20%7C%20EPICv2%20%7C%20HM450%20%7C%20MSA-brightgreen)](README.md#data-files-and-the-store)
[![vs R](https://img.shields.io/badge/betas%20vs%20R-bit--identical-success)](NUMERICS.md)

Infinium DNA-methylation analysis as a single C binary — **IDAT → betas → QC,
differential methylation, and copy number** — with no R and no Bioconductor.
Sesame2 is the second implementation of [sesame](https://github.com/zwdzwd/sesame) /
`openSesame`, validated against the R package as a permanent oracle.

> **Sesame2 (2.x, this C command line) and sesame (1.x, the R/Bioconductor
> package) are parallel, not sequential.** Neither is an upgrade of the other:
> one is an R package you `library()`, the other a binary you run in a shell.
> They share a name because they are one project and one method, and they are
> installed from different places — `BiocManager::install("sesame")` for the R
> package, `conda install -c zhou-lab sesame` for this one. The R package
> stays on the 1.x series.

<img src="docs/figures/pipeline.svg" alt="sesame-cli workflow" width="100%">

📖 **[Documentation site](https://zwdzwd.github.io/sesame-cli/)** ·
🔬 [Fidelity & numerics](NUMERICS.md)

## What it does

- **Read IDATs** — `readIDATpair`, **bit-identical** to R on 33/33 test files
  (~20.4M bead records, 9 platforms, plain and gzipped).
- **Preprocess to betas** — the full `QCDPB` pipeline (qualityMask,
  inferInfiniumIChannel, dyeBiasNL, pOOBAH, noob); `prep=""` betas are
  bit-identical to R, the rest match to documented tolerances.
- **QC** — the 66-metric `sesameQC` panel, one row per sample.
- **Differential methylation** — per-probe OLS with t/F tests and BH adjustment
  (`DML`), matching R's `lm` to ~1e-9.
- **Copy number** — log2 ratios vs a normal panel, genome binning, and a
  deterministic port of DNAcopy's **circular binary segmentation**.
- **SNP genotyping** — genotype the array's `rs` probes into a VCF
  (`formatVCF`); genotype calls and variant fractions are exact vs R.
- **Scale** — a cohort runs in one process across a thread pool, streaming to
  compact YAME `.cg` files; a bad IDAT never aborts the run.

Arrays: **EPIC, EPICv2, HM450, MSA**. Outputs are [YAME](https://github.com/zhou-lab/YAME)
`.cg` files, readable with `yame` or `libsesame`.

## Install

```sh
# conda (recommended) — the package is `sesame`; deps are zlib, libcurl and yame
conda install -c zhou-lab -c conda-forge sesame yame
```

Or build from source (a C compiler + `make`, `zlib`, `libcurl`, and the bundled
**YAME** submodule):

```sh
git clone --recurse-submodules https://github.com/zwdzwd/sesame-cli
cd sesame-cli && make          # builds YAME/htslib from the submodule, then ./sesame
cp sesame ~/.local/bin/         # put it on PATH
```

The binary is `sesame`. At runtime it needs only zlib and libcurl (YAME/htslib
are linked statically). You also want the [`yame`](https://github.com/zhou-lab/YAME)
binary on `PATH`: sesame does not download anything, so `yame fetch` is how the
annotation it reads gets there. `sesame version` names the yame release this
build expects.

## Quickstart

```sh
# 1. one-time: fetch a platform's annotation (ordering + mask + coords + SNP)
#    with yame, which owns downloading for the whole tool suite
yame fetch -y EPICv2

# 2. preprocess a cohort (default QCDPB) -> one indexed .cg per output + qc.tsv
#    each arg is an IDAT prefix, or a directory searched recursively for pairs
sesame preprocess --out out/ idats/
#    -> out/beta.cg  out/intensity.cg  out/pval.cg  out/qc.tsv

# 3a. differential methylation from the betas
sesame dml --betas out/beta.cg --index "$STORE/EPICv2/EPICv2.ordering.tsv.gz" \
           --meta samples.tsv --formula '~ group + age' > dml.tsv

# 3b. copy number for a tumor sample (writes segments + bins)
yame fetch -y hg38
sesame preprocess --prep "" --raw-signal --output total_intensity --out t/ tumor
sesame cnv --platform EPICv2 --normals EPICv2.cnvnormals.cg \
           t/total_intensity.cg out/segments.tsv out/bins.tsv
```

sesame is entirely offline -- it has no network code at all. `yame fetch` is
the one download path, and it fills a store every tool in the suite reads, so a
file fetched once serves sesame, kycg and yame alike. `sesame version` names
the yame this binary was built against and the store it reads.

## Commands

```
sesame preprocess   [options] <prefix|dir> [...]       # IDAT -> YAME .cg (+ qc.tsv)
sesame dml          --betas <beta.cg|matrix.tsv> (--formula .. --meta .. | --design ..)
sesame cnv          [--probes] [--platform P] <target.cg> <segments.tsv> <bins|probes.tsv>
sesame vcf          <prefix> [--platform P] --snp <snp.tsv.gz>     # SNP genotypes -> VCF
sesame mliftover    --to <platform> --platform <src> <in.cg> <out.cg>  # cross-platform betas (mLiftOver)
sesame impute       --method <mean|neighbors> <in.cg> <out.cg>    # fill NA betas
sesame region       (<chr:beg-end>|--gene NAME) --betas <beta.cg> # region betas -> long-form TSV
sesame deidentify   [--randomize|-r --seed N] <prefix|idat>        # strip SNP fingerprint from IDAT
sesame attach-probe [--platform P] <file.cg|.cm|.tsv>  # label a positional file with Probe_IDs
sesame idat-dump    [--head N] [--tsv] <file.idat[.gz]># inspect a raw IDAT
sesame version                                         # version, the yame it links, the store
```

---

## Preprocess — IDAT → betas

`preprocess` applies `--prep` (default `QCDPB`, `openSesame`'s default) to each
sample and writes one **YAME `.cg`** per requested output over the whole cohort —
one indexed file, one block per sample — into `--out DIR` (default `.`), plus a
`qc.tsv`.

<img src="docs/figures/beta_dist.svg" alt="beta-value distribution" width="66%">

| `--output` (comma list) | file | format |
|---|---|---|
| `beta` | `beta.cg` | 4 — float32 beta; NA = masked |
| `intensity` | `intensity.cg` | 3 — M/U integers (`yame` derives beta **and** coverage) |
| `total_intensity` | `total_intensity.cg` | 4 — M+U float |
| `pval` | `pval.cg` | 4 — pOOBAH detection p |
| `qc` | `qc.tsv` | the 66-metric `sesameQC` panel, one row per sample |

Default `--output` is `beta,intensity,pval,qc`.

```sh
sesame preprocess --out out/ idats/                           # a folder, searched recursively
sesame preprocess --output beta --out out/ pfx1 pfx2          # explicit prefixes -> just beta.cg
sesame preprocess --output total_intensity --raw-signal --out cnv/ tumor  # CNV input
sesame preprocess --prep "" --out raw/ pfx1                   # raw (bit-identical to R) betas
```

Signal outputs (`intensity`, `total_intensity`) reflect the prep; `--raw-signal`
takes them from the raw signal (what CNV wants). Other flags: `--index`/`--platform`
(else auto-detected from bead count), `--min-beads N`, `--threads N`, `--tmp DIR`.
Samples run in parallel; a failed sample becomes an NA block and the exit code is 1.
See [Batch mode](#batch-mode) and [The QCDPB pipeline](#the-qcdpb-pipeline).

## Differential methylation — `dml`

Per-probe **differential methylation**: for each probe, an OLS of its betas across
samples on a design, giving per-coefficient estimates and t-tests, a holdout
F-test per categorical variable, effect sizes, and BH-adjusted p-values — sesame's
`DML` / `summaryExtractTest`, as a TSV (one row per probe).

```sh
sesame dml --betas out/beta.cg --index "$STORE/EPICv2/EPICv2.ordering.tsv.gz" \
           --meta samples.tsv --formula '~ group + age' > dml.tsv
```

`--betas` is a `preprocess` `beta.cg` (with `--index <ordering>` to resolve probe
IDs, since a `.cg` is positional; `dml` takes no `--platform` and does not look
in the store, so give the ordering by path — and it still calls `--index`
deprecated, which is the only flag it has) or a `Probe_ID` matrix TSV; `--meta`'s first
column matches the sample names. `--formula` takes **main effects** (categorical
columns auto-dummied to match R's `model.matrix` — treatment contrasts,
alphabetical levels; continuous as-is, with an intercept). For interactions or
splines, build the design elsewhere and pass `--design <numeric.tsv>` (first
column = sample id). Runs in parallel (`--threads`). Because it consumes a betas
matrix, there is no data-lineage caveat — it matches R's `lm` to ~1e-9.

## Copy number — `cnv`

For each sample in the target `total_intensity.cg`, regress its per-probe total
intensity on a panel of normals (OLS), then `log2(target / max(fitted, 1))` per
probe — sesame's `cnSegmentation`. Probes are binned along the genome (50 kb tiles
minus assembly gaps, merged to ≥ `--min-probes`), each bin takes the median probe
signal, and the bins are segmented by **circular binary segmentation** (a
deterministic port of DNAcopy's CBS).

<img src="docs/figures/cnv_k562.png" alt="K562 copy-number profile" width="100%">

```sh
# <target.cg> <segments.tsv> <bins.tsv> are all required, positional:
sesame cnv --platform EPICv2 t/total_intensity.cg out/segments.tsv out/bins.tsv
# --probes writes per-probe log2 to the detail file instead of per-bin:
sesame cnv --platform EPICv2 --probes t/total_intensity.cg out/segments.tsv out/probes.tsv
```

`--coords` and the ordering default to the store for `--platform`
(`<store>/<P>/<P>.hg38.coord.tsv.gz`); genome tiling comes
from `yame fetch genomes/<--genome>` (default hg38). `--normals` is
user-supplied -- the normal panel is not published with the annotation. Both output
files carry a leading `sample` column — segments: chrom/start/end/nbin/seg.mean;
bins: chrom/start/end/nprobes/log2ratio; probes: Probe_ID/chrom/pos/log2ratio.
The fit + log2 matches R's `lm` to ~1e-5 on identical inputs; CBS recovers 97% of
DNAcopy's breakpoints with exact locations (DNAcopy's permutation significance is
itself unseeded, so an exact match is impossible — see `NUMERICS.md`). The figure
above is real `sesame cnv` output for K562.

## SNP genotyping — `vcf`

Genotype the array's SNP (`rs`) probes and write a sites-only VCF — sesame's
`formatVCF`. Per SNP probe it forms a variant fraction (`1-beta` for an ALT bead,
the out-of-band Type-I allele fraction for a channel-switch probe, else `beta`),
then a binomial model calls `0/0`, `0/1`, or `1/1` with a phred-like score in the
`INFO` field.

```sh
sesame vcf tumor --platform EPICv2 > geno.vcf              # SNP table from the store
sesame vcf tumor --platform EPICv2 --variants > geno.vcf   # informative probes only
```

It uses the **raw** signal (channel inference would erase the Type-I channel-switch
genotype). `--snp` is the platform's `<platform>.<genome>.snp.tsv.gz` (or the store
default). `--variants` keeps only the informative probes — rs and channel-switching
— dropping the non-switching Infinium-I bulk (~99% of rows) that almost never call
a variant. Genotype calls and variant fractions are **exact** vs R's `formatVCF` on every
probe whose colour channel the store's ordering and `sesameData`'s manifest
agree on — 127,550 of 127,572 on EPICv2. The other 22 are an annotation
disagreement rather than a divergence, and the quality score differs on ~0.01%
of probes in the deep tail; both are quantified in `NUMERICS.md`.

## Utilities

**`attach-probe`** — a YAME `.cg`/`.cm` and the per-probe coordinate tables are
**positional** (one row per probe in ordering order, no `Probe_ID` column inside).
This prepends the ordering's `Probe_ID` and prints a labeled TSV — greppable,
joinable:

```sh
sesame attach-probe --platform MSA MSA.hg38.coord.tsv.gz         # coord + Probe_ID
sesame attach-probe --platform MSA --all MSA.hg38.mask.cm        # every mask set, per probe
sesame attach-probe --index ord.tsv.gz --beta out/intensity.cg  # M/U -> beta
```

The row count **must** match the ordering — a mismatch is a hard error, never
silent misalignment.

**`idat-dump`** — inspect a raw `.idat`/`.idat.gz`: a summary header, or `--tsv`
for `addr<TAB>mean<TAB>sd<TAB>nbeads` (`--head N` limits rows).

**`version`** — the sesame version, the yame it was built against, the
annotation tags that come with that yame, and the store path in use. There is
no `fetch` subcommand: annotation is downloaded by `yame fetch`. See
[the store](#data-files-and-the-store).

---

## The QCDPB pipeline

`--prep` applies steps in the **order given**; `QCDPB` is the `openSesame` default.
Each is an independent, self-contained port of the R function (see `NUMERICS.md`).

| code | step | what it does |
|---|---|---|
| `Q` | `qualityMask` | OR the platform's recommended probe mask (from the `.cm`) into the SigDF mask |
| `C` | `inferInfiniumIChannel` | reassign each Infinium-I probe to its brighter color channel |
| `D` | `dyeBiasNL` | nonlinear dye-bias correction — pull the red and green Inf-I distributions to a common curve |
| `E` | `dyeBiasL` | linear dye-bias correction — scale each channel to a common reference intensity |
| `P` | `pOOBAH` | detection-p masking from the out-of-band + negative-control background; mask probes with `p > 0.05` |
| `B` | `noob` | normal-exponential background subtraction |

Extra `preprocess` options: `--collapse` averages EPICv2/MSA replicate probes to
their cg-prefix in the beta output (`betasCollapseToPfx`; axis → `beta.cg.probes`);
`--detection pneg` emits the negative-control-ECDF detection p (`detectionPnegEcdf`)
in `pval.cg` instead of pOOBAH. The `qc.tsv` panel adds `GCT`
(`bisConversionControl`) on every platform (annotation v8.1+).

Pass any subset in any order (`--prep CD` for channel + dye bias, `--prep ""` for
raw betas). `Q`, `P`, and `B` read the `.cm` mask, so the platform must be in the
store.

## Batch mode

`preprocess` over **multiple prefixes** runs them in one process. The index and
masks are parsed **once** and shared read-only; per-sample work runs across a
pthread pool (`--threads`, default = online CPUs).

- **One platform per batch** (the row space of the `.cg`).
- Each output is **one indexed `.cg`** — one block per sample, named by basename.
- **Deterministic:** a sample's block is its command-line position, so
  `--threads 1` is byte-identical to `--threads N`.
- **Robust:** a failed sample (missing/corrupt IDAT, wrong platform) becomes an
  all-`NA` block with a warning and sets [exit status](#exit-status) 1.
- **Scales past RAM:** each result store is an unlinked temp file the OS pages to
  disk (`--tmp DIR`), so a cohort larger than memory still works.

```sh
sesame preprocess --prep QCDPB --threads 8 --out out/ run/*prefixes
```

## Output formats

Numeric outputs are YAME `.cg` files (a BGZF container of per-sample blocks) plus a
`<file>.idx` mapping sample names to block offsets — read with `yame unpack` /
`yame info`, or `libsesame`'s `sesame_read_cg`.

- **format 3** (`intensity.cg`): M and U per probe as integers; `yame` derives
  both beta (`MU2beta`) and coverage (`MU2cov`), exact for raw IDAT signal.
- **format 4** (`beta.cg`, `total_intensity.cg`, `pval.cg`): one float32 per
  probe, NA as a negative value.

`qc.tsv` is a plain TSV, one row per sample. Betas are float32 in the `.cg`
(biologically lossless); the validation gates run through the double-precision
library path so they stay bit-exact against R.

## Data files and the store

Ordering tables, masks, per-probe coordinates and SNP tables live in
[`zhou-lab/InfiniumAnnotation`](https://github.com/zhou-lab/InfiniumAnnotation),
one folder per platform, versioned by **git tag**. **sesame does not download
them** — [`yame`](https://github.com/zhou-lab/YAME) owns fetching for the whole
tool suite, and `yame fetch -y EPICv2` mirrors that folder into
the shared **store**:

```
<store>/EPICv2/SHA256SUMS          byte-identical copy of the remote
<store>/EPICv2/EPICv2.ordering.tsv.gz
<store>/EPICv2/EPICv2.hg38.mask.cm(.idx)
<store>/EPICv2/EPICv2.hg38.coord.tsv.gz
<store>/EPICv2/EPICv2.hg38.snp.tsv.gz
<store>/EPICv2/KYCG/...            kycg's per-platform sets
<store>/hg38/seqinfo.tsv.gz  gaps  cytoband  genes.bed.gz(.tbi)
```

so `cd <store>/EPICv2 && shasum -a 256 -c SHA256SUMS` verifies it by hand. The
path is keyed by the browser hierarchy — species, then platform or genome
build — rather than by which tool asked, which is what makes one download
serve sesame, kycg and yame at once. Note it is not the same string as the
fetch address: `yame fetch` names a unit `<source>/<target>`, because one
store directory can be filled from more than one upstream (`<store>/hg38`
holds the genome annotation, and `<store>/hg38/KYCG` kycg's sets). Genome-level annotation comes the same way from
[`zhou-lab/genomes`](https://github.com/zhou-lab/genomes) via `yame fetch -y
hg38` (hg38, mm10, mm39 published); it drives CNV binning, the CNV
ideogram, and `region --gene`.

**Store location:** `$YAME_DATA_HOME`, else
`${XDG_DATA_HOME:-~/.local/share}/yame` — one variable for the whole suite,
because a shared store only dedupes if every tool agrees where it is. It is the
*data* tier rather than a cache tier on purpose: these are large references that
should survive somebody clearing `~/.cache`. **Asset lookup:** an explicit `--index` /
`--coords` / `--normals` / `--snp` is a **path**, used as given — there is no
store lookup behind it, so a bare `--snp EPICv2.hg38.snp.tsv.gz` fails unless
that file is in the working directory. Omit the flag and the chain
`<store>/<platform>/<file>` → `./<file>` resolves it.

**sesame has no network code at all** — it never downloads, never prompts (a
prompt would hang a Nextflow job or Docker build silently), and never falls back
to the network when an asset is absent. You get an error naming the exact `yame
fetch` to run and the annotation tag this build expects.

One asset is deliberately not published: the CNV normal panel
(`<platform>.cnvnormals.cg`). Build one with `make cnv-normals`, or pass your
own with `cnv --normals`.

### Custom arrays — the ordering table

*(the same ground with worked examples:
[docs → Custom array](https://zwdzwd.github.io/sesame-cli/#custom-array))*

`--platform` and `--index` take **the same thing**. `--platform EPICv2` is only a
lookup shortcut: it finds `<store>/EPICv2/EPICv2.ordering.tsv.gz` for you. So a
custom or pre-release array needs no new flag — hand `--platform` the ordering
that describes it, and a path is used as a path:

```sh
sesame preprocess  --platform MyArray.ordering.tsv.gz --out out/ idats/
sesame attach-probe --platform MyArray.ordering.tsv.gz out/beta.cg
```

**Tab or comma, either way.** The separator is sniffed from the header line, so
an ordering exported straight out of a CSV manifest needs no conversion. A comma
only wins if the header has no tab.

Every companion asset has its own flag, inferred when `--platform` is a name and
named explicitly otherwise — so a custom array is not second-class:

| flag | asset | inferred? |
|---|---|---|
| `--mask` | `.cm` quality/background mask | yes, for a named platform |
| `--coords` | per-probe coordinates | yes |
| `--snp` | SNP table for `vcf` | yes |
| `--normals` | CNV normal panel | no — not published, build with `make cnv-normals` |

```sh
# a custom array, full QCDPB: the ordering names the array, the mask comes with it
sesame preprocess --platform MyArray.ordering.csv --mask MyArray.mask.cm \
    --prep QCDPB --out out/ idats/
```

**Making the mask.** A `.cm` is a bitset with one bit per ordering row, in the
same order — 1 means drop the probe:

```sh
# one 0/1 per ordering row
# replace the condition -- an undefined name is 0 to awk, so this writes
# a mask of all zeros that packs and indexes fine and masks nothing
awk -F'\t' 'NR>1{print ($1 ~ /^ctl_/) ? 1 : 0}' MyArray.ordering.csv > m.txt

yame pack -f b m.txt MyArray.mask.cm     # -f b = binary/bitset
echo M_custom > names.txt
yame index -s names.txt MyArray.mask.cm  # name the track
```

Given directly, `--mask` wants a **single-track** `.cm`. A published platform's
mask holds many (EPICv2 has 25) and sesame selects the recommended subset *by
platform name*, so a file handed over without a name is ambiguous — it would
have to union all 25, masking 111k probes where the recommended 5 mask 50k.
sesame refuses and says so instead of guessing.

Without `--mask`, a custom array can still run `--prep C` or `""` — with
`--output beta`. The prep code is not the whole rule: the default `--output`
also asks for `pval` and `qc`, which need a mask whatever the prep, so
`--prep C` alone still stops with "P/B/pval/qc need a mask". `Q`, `P` and `B`
need one too. When a named platform ships more than one `.cm` (MM285 has mm10
and mm39) the lexicographically first is used — deterministic, and mm10 for
MM285 — but name `--mask` when the build matters.

> `--index` was the old name for `--platform <ordering>` and still works, with a
> deprecation notice. "index" read like a search structure; the ordering is the
> array's identity, which is what `--platform` means.

**It is not the Illumina/sesame manifest with tabs instead of commas.** It is a
four-column *ordering* derived from a manifest — everything sesame needs to turn
bead addresses into probes, and nothing else. No sequences, no coordinates:
those live in the separate `coord`/`mask`/`snp` files, which a custom array
simply does without.

| column | meaning |
|---|---|
| `Probe_ID` | any string; what `attach-probe` prints and what joins to your annotation |
| `M` | the **methylated** bead address, decimal — Illumina's `AddressB_ID`. `NA` for Infinium-II |
| `U` | the **unmethylated** bead address — Illumina's `AddressA_ID`. Always present |
| `col` | design and in-band channel: `G` Infinium-I green, `R` Infinium-I red, `2` Infinium-II |

Rules the parser enforces:

- **Tab-separated, with a header line.** The first line is skipped unconditionally,
  so it must be there; its text is not checked.
- **Four columns**, or five with a trailing `mask` (`0`/`1`). The five-column form
  is the legacy sesameData export; the published orderings use four and keep the
  mask in a companion `.cm`.
- `col` must be exactly `G`, `R` or `2` — anything else is a hard error naming the row.
- Addresses are decimal or the literal `NA`; `NA` becomes "no address".
- Plain or gzipped, both fine — `MyArray.ordering.tsv` and `.tsv.gz` read the same.

**Row order is load-bearing, and permanent.** A `.cg` stores one value per row in
this order with no Probe_IDs inside, so row *i* is probe *i* forever. Sorting or
inserting a row silently re-labels every `.cg` ever written against the old
file — which is what the row-count check catches when it says "lineage mismatch".
Freeze the order when you first publish an ordering, and version it if it changes.

From an Illumina manifest, that is a projection:

```
Probe_ID  <- IlmnID
M         <- AddressB_ID   (empty -> NA, which is every Infinium-II probe)
U         <- AddressA_ID
col       <- Infinium_Design_Type II            -> 2
             Infinium_Design_Type I + Color_Channel Grn -> G
             Infinium_Design_Type I + Color_Channel Red -> R
```

Sanity-check a new file before shipping it: every `col=G|R` row should have both
addresses, every `col=2` row should have `M=NA`, and the row count should equal
the number of probes you expect. `sesame attach-probe --index MyArray.ordering.tsv.gz`
on any `.cg` you built with it is the end-to-end check — it fails loudly if the
lengths disagree.

**Versions travel together.** `sesame version` (and the top of `sesame help`)
names the yame this binary was built against and the store it will read:

```
sesame <version>
    built against  YAME <version>
    store          /home/you/.local/share/yame   ($YAME_DATA_HOME unset; yame fetch fills it)
```

sesame links `libyame.a` from a pinned submodule, so that version is a fixed
fact of the build — and with it the store layout and the annotation tags, which
come from yame's shared catalog (`YAME/tools/registry/TAGS`) rather than from
sesame. Every tool in the suite therefore pins the same annotation generation;
`src/registry.h` is regenerated from that catalog with `make registry`. Use a
`yame` at least as new as the one named here.

## Configuration

Flags plus a few environment variables; no config file.

| variable | effect |
|---|---|
| `YAME_DATA_HOME` | the shared store — where `yame fetch` writes and sesame looks up orderings, masks, coordinates and SNP tables |
| `XDG_DATA_HOME` | fallback store root (`$XDG_DATA_HOME/yame`) when `YAME_DATA_HOME` is unset |
| `TMPDIR` | where a batch's temp-file-backed result matrix is created |
| `SESAME_TEST_IDATS` | test IDAT directory for `make test` (default `~/repo/InfiniumTestIDATs`) |

## Exit status

| code | meaning |
|---|---|
| `0` | success |
| `1` | a runtime error, or (in a batch) at least one sample failed — the matrix still prints, with `NA` for the failures |
| `2` | usage error (bad flags/arguments) |

## Fidelity and validation

R is the oracle, permanently. Precision-sensitive gates run through a harness
(`tests/pipeline_dump`) that dumps the library pipeline at **double precision**, so
they stay bit-identical / ULP-exact even though the `.cg` product stores betas as
float32. The golden ladder (`make test`):

| level | gate | status |
|---|---|---|
| IDAT reader | `IlluminaID`/`Mean`/`SD`/`NBeads` **bit-identical** | ✅ 33/33 files, ~20.4M records |
| `C` channel | channel calls identical | ✅ 5/5 |
| `D` dyeBias | max relative Δ ≤ 1e-12 (clean-room qnorm) | ✅ ~2 ULP |
| `Q` mask | recommended-track union of the `.cm` | ✅ self-consistent (mask lineage) |
| `P` pOOBAH | every C-vs-R flip is 0.05-boundary/D2 | ✅ 12 flips, all at the cutoff |
| `B` noob | `normExpSignal`/`huber` on identical inputs | ✅ ≤ few ULP; betas lineage-only |
| Betas `prep=""` | **bit-identical** | ✅ 6/6, 4.4M betas |
| Batch | `--threads 1 == N`; bad sample → NA block | ✅ ThreadSanitizer-clean |
| QC panel | every `sesameQC` metric within lineage scale | ✅ 65/65 |
| DML | vs R `DML`/`summaryExtractTest` | ✅ ~5e-10 |
| CNV | fit+log2 vs R `lm`; CBS breakpoints vs DNAcopy | ✅ ~1e-5; 97% recall, exact locations |
| VCF | genotypes vs R `formatVCF` | ✅ 127572 probes, GT + fraction exact |

"Lineage" means the published ordering/mask is a newer data version than the
installed `sesameData`, so a few probes differ for reasons that predate any
arithmetic — each quantified in `NUMERICS.md`, which carries one row per
intentional divergence from R. The IDAT parser (untrusted GEO input) is fuzzed and
run under ASan/UBSan.

## Development

```sh
make                # build ./sesame
make test           # the full golden ladder vs the R oracle (+ the docs check)
make test-docs      # run every documented example against the built binary
make test-docs-check# versions/subcommands/flags/metric counts vs the binary
scripts/coverage.sh # gcov line coverage of `make test` over src/ and cli/
make asan           # rebuild under ASan + UBSan
make fuzz-replay    # replay the IDAT-parser corpus under ASan/UBSan
make clean
```

`make test-docs` runs `docs/examples/*.sh` — the page's own commands, in a
sandbox, on this checkout's binary. It needs the store and the test IDATs, so
it is not part of `make test` and never runs in CI. `make test-docs-check`
needs neither and does: it is what fails when a doc quotes a flag the binary
does not have. Releases follow the
[release SOP](https://github.com/zhou-lab/labjournal) entry beside the project
org (`20210109_sesame_cli_RELEASE_SOP.md`).

`make test` needs `Rscript` with the sesame R package (the oracle) and test IDATs
at `$SESAME_TEST_IDATS`. Layout: `cli/` (the command), `src/` (`libsesame` — IDAT
reader, index, SigDF, prep steps, numerics, CNV/CBS, fetch), `include/sesame.h`
(public API), `tests/` (shell drivers + R oracles), `docs/` (the Pages site),
`YAME/` (submodule). `NUMERICS.md` documents every intentional divergence from R
and is required reading before changing a prep step.

## License

**2-Clause BSD for academic and non-profit research use**; for commercial use
or transfer, inquire with Dr. Wanding Zhou at zhouw3@chop.edu. © 2026-present
The Children's Hospital of Philadelphia — see `LICENSE`. The bundled
[YAME](https://github.com/zhou-lab/YAME) carries the same terms, which is what
lets sesame-cli link it directly to read `.cm` masks in-process. The R package
`sesame` is a *separate program* distributed through Bioconductor under
AGPL-3; invoking either way is unaffected. The dye-bias quantile normalization
is a clean-room reimplementation rather than `preprocessCore` (`LGPL`); it
agrees to ~2 ULP (see `NUMERICS.md`, D8).
