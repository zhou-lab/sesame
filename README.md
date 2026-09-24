<h1 align="center">SeSAMe2</h1>

<p align="center">
  <a href="https://github.com/zhou-lab/sesame/actions/workflows/conda-build.yml"><img alt="build" src="https://github.com/zhou-lab/sesame/actions/workflows/conda-build.yml/badge.svg"></a>
  <a href="https://anaconda.org/zhou-lab/sesame"><img alt="conda" src="https://img.shields.io/conda/vn/zhou-lab/sesame?label=conda"></a>
  <a href="LICENSE"><img alt="license" src="https://img.shields.io/badge/license-BSD--2--Clause%20(academic)%20%2F%20commercial-blue.svg"></a>
  <a href="scripts/coverage.sh"><img alt="coverage" src="https://img.shields.io/endpoint?url=https%3A%2F%2Fzhou-lab.github.io%2Fsesame%2Fcoverage.json"></a>
  <a href="https://zhou-lab.github.io/sesame/"><img alt="docs" src="https://img.shields.io/badge/docs-online-blueviolet"></a>
  <a href="include/sesame.h"><img alt="language" src="https://img.shields.io/badge/C-C11-00599C"></a>
  <a href="https://zhou-lab.github.io/sesame/"><img alt="arrays" src="https://img.shields.io/badge/arrays-EPIC%20%7C%20EPICv2%20%7C%20HM450%20%7C%20MSA-brightgreen"></a>
  <a href="NUMERICS.md"><img alt="betas vs R" src="https://img.shields.io/badge/betas%20vs%20R-bit--identical-success"></a>
</p>

<p align="center">
  <b>Infinium DNA-methylation analysis as a single C binary</b><br>
  IDAT &rarr; betas &rarr; QC, differential methylation, copy number, SNP
  genotyping. No R, no Bioconductor, no network.
</p>

<p align="center">
  &#128214; <b><a href="https://zhou-lab.github.io/sesame/">Documentation</a></b> &middot;
  &#128300; <a href="NUMERICS.md">Fidelity &amp; numerics</a> &middot;
  &#128269; <a href="DIVERGENCES.md">Where we differ from R</a>
</p>

SeSAMe2 is the second implementation of
[sesame](https://github.com/zwdzwd/sesame) / `openSesame`, validated against the
R package as a permanent oracle.

> **SeSAMe2 (2.x, this command line) and sesame (1.x, the R/Bioconductor
> package) are parallel, not sequential** — one method, two implementations.
> R: `BiocManager::install("sesame")`; shell: `conda install -c zhou-lab sesame`.

## Install

```sh
conda install -c zhou-lab -c conda-forge sesame yame
```

Or from source — a C compiler, `make`, `zlib`, `libcurl`, and the bundled YAME
submodule:

```sh
git clone --recurse-submodules https://github.com/zhou-lab/sesame
cd sesame && make
```

The binary is `sesame`. Keep [`yame`](https://github.com/zhou-lab/YAME) on
`PATH` too: sesame has no network code, so `yame fetch` is how the annotation it
reads gets there, into a store the whole tool suite shares. `sesame version`
names the yame release this build expects and the store it reads.

## Use

```sh
yame fetch -y EPICv2                        # once, per platform
sesame preprocess --out out/ idats/         # -> beta.cg intensity.cg pval.cg qc.tsv
```

`sesame help` lists every command, and `sesame <command> -h` is the authority on
its flags. Worked examples for each — preprocessing, DML, CNV, SNP genotyping,
region views, custom arrays, and reading the outputs back — are on the
**[documentation site](https://zhou-lab.github.io/sesame/)**.

## Fidelity

R is the permanent oracle. Raw betas are bit-identical; every prep step is gated
against R to a tolerance recorded in [NUMERICS.md](NUMERICS.md), and the handful
of places where the two genuinely answer different questions — all of them
annotation-lineage differences, none a tolerance — are written up in
[DIVERGENCES.md](DIVERGENCES.md).

## Development

```sh
make                 # build ./sesame
make test            # the golden ladder vs the R oracle
make test-docs       # run every documented example against this binary
make test-docs-check # versions/subcommands/flags vs the binary (no oracle needed)
scripts/coverage.sh  # gcov line coverage
make asan            # rebuild under ASan + UBSan
```

`make test` needs `Rscript` with the sesame R package and test IDATs at
`$SESAME_TEST_IDATS`; CI runs neither, so the numeric gates are the release
manager's job. Layout: `cli/` the command, `src/` the library, `include/sesame.h`
the public API, `tests/` shell drivers with R oracles, `docs/` the Pages site,
`YAME/` the submodule. Read [NUMERICS.md](NUMERICS.md) before changing a prep
step. Releases follow `20210109_sesame_cli_RELEASE_SOP.md` in the lab journal.

## License

**2-Clause BSD for academic and non-profit research use**; for commercial use or
transfer, inquire with Dr. Wanding Zhou at zhouw3@chop.edu. © 2026-present The
Children's Hospital of Philadelphia — see [LICENSE](LICENSE). The bundled
[YAME](https://github.com/zhou-lab/YAME) carries the same terms. The R package
`sesame` is a *separate program* distributed through Bioconductor under AGPL-3.
