# Conda recipe for sesame

Builds the `sesame` binary and publishes it to the **zhou-lab** channel.

## Install (end users)

```sh
conda install -c zhou-lab -c conda-forge sesame
```

This installs the `sesame` command. The package was called `sesame-cli` up to
and including 0.3, to keep its distance from the Bioconductor R package. That
distance is no longer wanted: the C implementation is the one we build on, and
the package, the binary and the method should all be spelled the same. There is
no actual name clash to avoid -- the R package is `bioconductor-sesame` on
bioconda, and nothing named `sesame` exists on conda-forge or bioconda.

**Renaming orphans the old name.** Anyone who installed `sesame-cli` keeps it
and will not see `sesame` updates; `conda install sesame-cli` still resolves,
but to 0.3 forever. If that matters, publish a `sesame-cli` metapackage that
depends on `sesame` so the old name keeps tracking the new one.

## Build & upload (maintainers)

Requires `conda-build` and `anaconda-client`:

```sh
conda install -n base conda-build anaconda-client

# Build (Linux + macOS; Windows is skipped). conda-forge supplies zlib + libcurl.
conda build -c conda-forge conda-recipe/

# Upload the built package to the zhou-lab channel
anaconda login
anaconda upload -u zhou-lab $(conda build -c conda-forge conda-recipe/ --output)
```

To upload automatically on every successful build:

```sh
conda config --set anaconda_upload yes
```

## Automated builds (CI)

`.github/workflows/conda-build.yml` runs `conda build` on every push/PR to `main`
across **linux-64** and **osx-arm64** (Intel `osx-64` is dropped — GitHub is
retiring its Intel macOS runners, so those jobs queue indefinitely), uploading
each package as a build artifact, and publishes to the `zhou-lab` channel on a
`vX.Y.Z` tag (guarded so the tag must match the recipe version). Publishing needs
an `ANACONDA_TOKEN` repository secret with upload rights to the org. Cutting a
release:

```sh
# bump `version` in conda-recipe/meta.yaml and SESAME_VERSION in include/sesame.h,
# commit, then:
git tag -a v2.0.0 -m "sesame 2.0.0" && git push origin v2.0.0
```

## Notes

- **In-repo build.** The recipe builds the checked-out tree (`source: path: ..`),
  including the vendored `YAME/` + `htslib` submodule, so a release build is just:
  check out the `vX.Y.Z` tag, then `conda build`. Keep `version` in `meta.yaml`
  and `SESAME_VERSION` in `include/sesame.h` in sync with the tag.
- **Submodule.** `YAME/` (and its nested `htslib/`) must be checked out before
  building (`git submodule update --init --recursive`); `path:` copies whatever
  is on disk. `htslib`'s `config.h` and `version.h` are committed, so no autoconf
  / git step is needed at build time.
- **Rebuilds the vendored static libs.** `YAME/libyame.a` and
  `YAME/htslib/libhts.a` are committed for the author's platform; `build.sh`
  deletes them (and the stale `.o`) so they rebuild for the target arch.
- **Flag handling.** The vendored Makefiles have rigid `CFLAGS` and ignore
  `CPPFLAGS`, so `build.sh` bakes conda's compile flags into `CC` (which every
  sub-make uses verbatim). `LDFLAGS` reaches the final link via the environment.
- **Dependencies.** Only zlib and libcurl (both from conda-forge, with
  `run_exports` that pin the runtime automatically). Make sure conda-forge is on
  your channel list when building.
- **glibc floor (do not remove).** `conda_build_config.yaml` pins
  `c_stdlib_version` to **2.17** on Linux (11.0 on macOS), and `meta.yaml` lists
  `{{ stdlib('c') }}`. Without the pin, `{{ compiler('c') }}` pulls
  `gcc_linux-64`, whose `sysroot_linux-64` dependency is unversioned, so the
  solver takes the newest sysroot (2.39) and the binary ends up needing symbols
  no ordinary machine has — that is how 0.2.0 build 0 shipped a linux-64 package
  that died with ``version `GLIBC_2.38' not found`` on RHEL 9 and Ubuntu 22.04.
  CI now fails the build if any required symbol is above `GLIBC_2.17`. Check a
  binary by hand with:

  ```sh
  objdump -T $(which sesame) | grep -oE 'GLIBC_[0-9.]+' | sort -uV | tail -1
  ```

- **Republishing a fixed build of an existing version.** Bump `number:` under
  `build:` in `meta.yaml` (not `version:`), then move the tag so CI rebuilds and
  uploads — the publish step's guard only requires tag == recipe version, and it
  uploads with `--force`:

  ```sh
  git tag -f -a v0.2.0 -m "sesame-cli 0.2.0 (build 1)" && git push -f origin v0.2.0
  ```

  Conda prefers the higher build number, so clients get the fixed package. Also
  delete the broken artifact from anaconda.org so it can never be selected:
  `anaconda remove zhou-lab/sesame-cli/0.2.0/linux-64/sesame-cli-0.2.0-h9bf148f_0.conda`
