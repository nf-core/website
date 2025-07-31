---
title: "1. Installation"
subtitle: Setting up nf-test in your development environment
weight: 10
---

## Installing nf-test

**For nf-core pipeline development and testing, you need to install nf-test separately.** nf-core tools only provides test commands for the `nf-core/modules` repository context.

### Prerequisites

Before installing nf-test, ensure you have:

- **Java 11 or higher**
- **Nextflow** (latest stable version recommended)

### Recommended Installation: Using Conda/Mamba

```bash
conda create -n nf-core -c bioconda nf-core nf-test
conda activate nf-core
# or
mamba create -n nf-core -c bioconda nf-core nf-test
mamba activate nf-core

### Alternative Installation Methods

For a standalone binary install: 

```bash
curl -fsSL https://get.nf-test.com | bash
```

Move the `nf-test` file to a directory accessible by your `$PATH` variable:

```bash
# Make it executable and move to PATH
chmod +x nf-test
sudo mv nf-test /usr/local/bin/
```

#### Install specific version

For installing a specific version of nf-test, you can specify this as an option to the install script.

```bash
curl -fsSL https://get.nf-test.com | bash -s 0.9.0
```

> **For complete installation instructions and troubleshooting**, visit the [official nf-test installation documentation](https://www.nf-test.com/docs/getting-started/).

### Verification

Verify your installation worked correctly, you can check the version is printed.

```bash
nf-test version
```

You can also list all available tests in a repository already with:

```bash
nf-test list
```

## Next Steps

Once you have nf-test installed, proceed to [nf-test Commands & Integration](./02_commands_integration.md) to learn the essential commands.
