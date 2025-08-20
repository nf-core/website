---
title: "Installation"
subtitle: Setting up nf-test in your development environment
weight: 10
---

## Installing nf-test

nf-test must be installed separately for nf-core pipeline development and testing.
The nf-core tools only provide test commands when working within the `nf-core/modules` repository context.

### Prerequisites

You will need the following to get started:

- **Java 11 or higher**
- **[Nextflow](https://www.nextflow.io/docs/latest/install.html)** (latest stable version recommended)
- **[nf-core/tools](https://nf-co.re/docs/nf-core-tools/installation)**

### Installation

nf-test is distributed via conda/mamba and as a self-installing package.

#### Conda/Mamba

To install nf-test with conda (or mamba):

1. Create an `nf-test` environment:

```bash
conda create -n nf-test -c bioconda nf-test
```

2. Activate your environment:

```bash
conda activate nf-test
```

3. Verify your install:

```bash
nf-test version
```

:::note
You can replace `conda` with `mamba` in the commands above if you prefer mamba's CLI, though modern conda (≥22.11) uses the same fast libmamba solver internally.
:::

#### Standalone binary

To install nf-test as a standalone binary:

1. Download the binary:

```bash
curl -fsSL https://get.nf-test.com | bash
```

2. Make `nf-test` executable:

```bash
chmod +x nf-test
```

3. Add `nf-test` to your `$PATH` variable. For example:

```bash
mkdir -p $HOME/.local/bin/
mv nf-test $HOME/.local/bin/
```

:::note
Ensure the directory `$HOME/.local/bin/` is included in your `PATH` variable.
Set `export PATH="$PATH:$HOME/.local/bin` to temporarily add this directory to `PATH`.
Add the export command to your shell configuration file, such as `~/.bashrc` or `~/.zshrc`, to add the directory to `PATH` permanently.
Alternatively, move the `nf-test` executable to a directory already in your `PATH`.
:::

## Next Steps

Once you have nf-test installed, proceed to [Repository Setup](./02_project_setup.md) to configure your testing environment.
