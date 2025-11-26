---
title: ARM64 on Bioconda
subtitle: Enable existing Bioconda and Conda-forge packages to work on ARM64 architecture
shortTitle: ARM64 on Bioconda
weight: 1
---

ARM64 support is increasingly important as cloud providers offer ARM-based instances (like AWS Graviton) that are 20-40% cheaper than x86 alternatives with comparable performance. Enabling ARM64 compatibility allows users to leverage these cost-effective platforms.

This page explains how to enable existing Bioconda and Conda-forge packages to work on ARM64 architecture (`linux/arm64`).

## Bioconda and Conda-forge

Bioconda is a conda channel for bioinformatics packages, while conda-forge provides general-purpose recipes. Both systems differ in CI infrastructure and recipe structure, though they share commonalities.

This documentation addresses adapting existing bioconda/conda-forge packages to support `linux/arm64` (ARM processors like AWS Graviton) alongside the standard `linux/amd64` (Intel architecture).

## What you need to do

The goal is to enable users to run `conda install <package>` and receive the package on ARM systems.

## Set up your environment

### Install Bioconda

To use Bioconda on ARM systems:

1. Install Conda:

   ```bash
   wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh
   sh ./Miniforge3-Linux-aarch64.sh
   ```

   Select "Yes" at installation completion.

2. Install Bioconda:

   ```bash
   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda config --set channel_priority strict
   ```

### Create a development environment (optional)

For local development on ARM machines:

1. Create a development environment:

   ```bash
   conda create -n bioconda-builds -c conda-forge -c bioconda bioconda-utils
   conda activate bioconda-builds
   ```

:::note
Local testing is not always reliable. CI builds often succeed when local attempts fail.
:::

## Enable packages

### Test locally

On ARM machines, test whether a package installs correctly:

```bash
conda install <name>
```

If the package fails to install, it needs ARM64 support.

### Test with Wave

If you do not have access to ARM hardware, use Wave through [Seqera Containers](https://seqera.io/containers/) or the [Wave CLI](https://github.com/seqeralabs/wave-cli):

```bash
wave --conda <name> --platform linux/arm64 --freeze --await
```

- `--freeze` persists generated images to Seqera Containers
- `--await` maintains execution until completion

Build failures display error messages with links to detailed build pages.

### Interpret test output

Build failures typically involve two categories:

1. "Not found" errors that require package fixes
2. Missing dependencies that need investigation

Most fixable packages reside in Bioconda. However, conda-forge packages occupy individual GitHub repositories under the [conda-forge organisation](https://github.com/conda-forge/), formatted as `<package-name>-feedstock`.

## Bioconda recipes

### Generic recipes

Generic packages may have missing non-generic dependencies that prevent functionality until dependencies are resolved.

Generic package indicators in `recipes/<name>/meta.yaml`:

```yaml
build:
  ..
  noarch: generic
  ..
```

When generic packages fail, dependencies require fixing. This is often resolved through simple build number increments and pull requests. Since generic packages lock to specific build versions that pre-date bioconda ARM support, they may need rebuilding.

### Non-generic recipes

Non-generic packages compile native code or download binaries (typically x86-only). Examine `meta.yaml` and adjacent `build.sh` files.

Enable ARM builds by adding to `meta.yaml`:

```yaml
extra:
  additional-platforms:
    - linux-aarch64
    - osx-arm64
```

Bump the build number by adding 1 to the existing value. Add a package versioning pin if missing:

```yaml
build:
  number: 2
  noarch: generic
  run_exports:
    - { { pin_subpackage(name, max_pin='x.x') } }
```

### Test in the CI

Open a pull request and monitor CI results. Common solutions include:

- **Source file patching**: `meta.yaml` handles patches
- **Platform-specific logic**: `build.sh` supports conditional bash logic. Use `CFLAGS`, `CC`, and `${PREFIX}` (dependency base directory containing `bin`, `lib`, `include`)
- **Build isolation**: CI output shows conda packages built in isolated environments

### Test builds locally

Build locally for faster iteration:

```bash
bioconda-utils build --packages <package>
```

Consult the [Bioconda dev instructions](https://bioconda.github.io/contributor/building-locally.html) for detailed guidance.

## Conda-forge

### Quick fixes

For conda-forge packages (for example, a package named `perl-nonsense`):

1. Navigate to `https://github.com/conda-forge/<package-name>-feedstock`
2. Check for existing ARM migration pull requests. If a stale PR exists, comment: `@conda-forge-admin please rerender`
3. If no outstanding PR exists, edit the [migrations list](https://github.com/conda-forge/conda-forge-pinning-feedstock/blob/main/recipe/migrations/arch_rebuild.txt) alphabetically. Bots will automatically initiate a PR.
4. When builds succeed, maintainers merge and binaries become available within a few days.
5. For unresponsive maintainers, comment: `@conda-forge-admin, please ping team` or `@conda-forge/core please help review and merge this PR`

### Manual fixes

Clone the feedstock:

```bash
git clone git@github.com:<your-github-username>/<package-name>-feedstock.git
cd <package-name>-feedstock
git remote add bot https://github.com/regro-cf-autotick-bot/<package-name>-feedstock
git fetch bot
git checkout -b aarch64-fixes bot/bot-pr_arch_<branch-suffix>
```

Replace `<branch-suffix>` with the actual migration branch name. Use tab completion after typing `bot/bot-pr_arch_` to find available branches.

Fix issues by editing `meta.yaml` or `conda-forge.yml`, then:

```bash
git commit --all
git push --set-upstream origin aarch64-fixes
```

Create a new PR from your branch. Request rerendering by commenting: `@conda-forge-admin please rerender`

:::note
Current delays affect builds on `linux-ppc64le` and `linux-aarch64` via emulation or Travis ARM resources. Modify `conda-forge.yml` to remove problematic platforms if necessary.
:::
