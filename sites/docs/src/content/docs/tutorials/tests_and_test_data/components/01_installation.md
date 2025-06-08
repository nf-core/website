---
title: "1. Installation"
subtitle: Setting up nf-test in your development environment
weight: 10
---

## Installation Options

### Option 1: Official Installation Script (Recommended)

The official and recommended way to install nf-test:

```bash
curl -fsSL https://get.nf-test.com | bash
```

If you don't have curl installed, you can use wget:

```bash
wget -qO- https://get.nf-test.com | bash
```

This will create the `nf-test` executable file in the current directory. Optionally, move the `nf-test` file to a directory accessible by your `$PATH` variable:

```bash
# Make it executable and move to PATH
chmod +x nf-test
sudo mv nf-test /usr/local/bin/
```

#### Install a specific version

To install a specific version, pass it to the install script:

```bash
curl -fsSL https://get.nf-test.com | bash -s 0.9.0
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
nf-test version
```

You should see output similar to:

```
🚀 nf-test 0.9.0
https://code.askimed.com/nf-test
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Nextflow Runtime:

      N E X T F L O W
      version 23.10.1 build 5891
      created 12-01-2024 22:01 UTC (23:01 CET)
      cite doi:10.1038/nbt.3820
      http://nextflow.io
```

## Troubleshooting Installation

### Nextflow Binary not found?

If you get an error message about Nextflow binary not being found, you have two options:

1. **Move Nextflow to PATH**: Ensure your Nextflow binary is in a directory accessible by your `$PATH` variable
2. **Set NEXTFLOW_HOME**: Set the environment variable `NEXTFLOW_HOME` to the directory containing the Nextflow binary:

```bash
export NEXTFLOW_HOME=/path/to/nextflow/directory
```

### Updating nf-test

To update an existing nf-test installation to the latest version:

```bash
nf-test update
```

## Common Plugins

For nf-core pipeline testing, you may need these plugins (loaded in your `nf-test.config`):

### Available nf-test Plugins

| Plugin Name | Description/Use |
|-------------|-----------------|
| **nft-utils** | Essential utility functions like `getAllFilesFromDir()` and `removeNextflowVersion()` - widely used across nf-core pipelines |
| **nft-bam** | Helper functionality for handling SAM/BAM/CRAM files during tests - critical for genomics pipelines |
| **nft-vcf** | Support for VCF files - utilities for testing variant call format files |
| **nft-fasta** | Support for FASTA files - enables validation and testing of FASTA file formats |
| **nft-fastq** | Support for FASTQ files - validation and testing utilities for sequencing data files |
| **nft-csv** | Support for CSV files - utilities for testing comma-separated value files |
| **nft-compress** | Support for ZIP files - enables testing of compressed file handling |
| **nft-anndata** | Support for AnnData (h5ad) files - utilities for testing single-cell analysis data formats |
| **nft-tiff** | Support for TIFF files - enables testing of Tagged Image File Format files |

### Plugin Configuration

These plugins are automatically loaded when specified in your `nf-test.config`:

```groovy
config {
    plugins {
        load "nft-utils@0.0.3"
        load "nft-bam@0.6.0"
        // Add additional plugins as needed:
        // load "nft-vcf@1.0.7"
        // load "nft-fastq@0.0.1"
    }
}
```

> **Note:** For the complete list of available plugins and their latest versions, visit the [nf-test plugins registry](https://plugins.nf-test.com/).

## Prerequisites

Before installing nf-test, ensure you have:

- **Java 11 or higher**
- **Nextflow** (latest stable version recommended)

## Next Steps

Once installed, proceed to [nf-test Commands & Integration](./02_commands_integration.md) to learn the essential commands. 