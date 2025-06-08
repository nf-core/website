---
title: "1. Installation"
subtitle: Setting up nf-test in your development environment
weight: 10
---

## Installation Options

### Option 1: Using Homebrew (macOS/Linux)

```bash
brew install nf-test
```

### Option 2: Using Conda/Mamba

```bash
conda install -c bioconda nf-test
# or
mamba install -c bioconda nf-test
```

### Option 3: Manual Installation

Download the latest version from the [nf-test releases page](https://github.com/askimed/nf-test/releases):

```bash
# Download and extract
curl -fsSL https://github.com/askimed/nf-test/releases/latest/download/nf-test.tar.gz | tar -xzf -

# Make executable and move to PATH
chmod +x nf-test
sudo mv nf-test /usr/local/bin/
```

## Verification

Verify your installation:

```bash
nf-test --version
```

## Essential Plugins

For nf-core pipeline testing, you'll commonly need these plugins (they are loaded in your `nf-test.config`):

- **nft-utils**: Provides utility functions like `getAllFilesFromDir()` and `removeNextflowVersion()`
- **nft-bam**: Provides BAM file testing utilities for pipelines that work with sequencing data

These plugins are automatically loaded when specified in your `nf-test.config`:

```groovy
config {
    plugins {
        load "nft-bam@0.5.0"
        load "nft-utils@0.0.3"
    }
}
```

## Prerequisites

Before installing nf-test, ensure you have:

- **Java 11 or higher**
- **Nextflow** (latest stable version recommended)

## Next Steps

Once installed, proceed to [nf-test Commands & Integration](./02_commands_integration.md) to learn the essential commands. 