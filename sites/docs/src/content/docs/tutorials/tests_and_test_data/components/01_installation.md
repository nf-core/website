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
conda install -c bioconda nf-core
conda install -c bioconda nf-test
# or
mamba install -c bioconda nf-core
mamba install -c bioconda nf-test
```

### Alternative Installation Methods

#### Official Installation Script

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

```bash
curl -fsSL https://get.nf-test.com | bash -s 0.9.0
```

> **For complete installation instructions and troubleshooting**, visit the [official nf-test installation documentation](https://www.nf-test.com/docs/getting-started/).

### Verification

Verify your installation:

```bash
nf-test version
```

List available tests:
```bash
nf-test list
```

---


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

## Next Steps

Once you have nf-test installed, proceed to [nf-test Commands & Integration](./02_commands_integration.md) to learn the essential commands.
