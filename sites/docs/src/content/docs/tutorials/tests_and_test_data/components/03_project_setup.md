---
title: "3. Repository Setup"
subtitle: Configuring your nf-core pipeline repository for testing with nf-test
weight: 30
---

## Example nf-core Repository Structure

Here is an example of an nf-core pipeline project with nf-test:

```
my-pipeline/
├── main.nf                          # Main pipeline script
├── nextflow.config                  # Main pipeline configuration
├── nextflow_schema.json            # Pipeline parameter schema
├── nf-test.config                  # nf-test configuration
├── modules.json                    # Module dependencies
├── tests/                          # Pipeline-level tests
│   ├── default.nf.test            # Main pipeline test
│   ├── params-1.nf.test           # Pipeline param-1 test
|   ├── params-2.nf.test           # Pipeline param-2 test
│   ├── nextflow.config            # Test-specific config
│   └── .nftignore                 # Files to ignore in snapshots
├── conf/                           # Configuration profiles
│   ├── test.config                # Test profile
│   ├── test_full.config           # Full test profile
│   ├── base.config                # Base configuration
│   └── igenomes.config            # Reference genome configurations
├── assets/                         # Pipeline assets
│   ├── samplesheet.csv            # Test samplesheet
│   └── schema_input.json          # Input validation schema
├── modules/nf-core/               # Module tests
│   └── */tests/main.nf.test      # Individual module tests
├── subworkflows/nf-core/          # Subworkflow tests
│   └── */tests/main.nf.test      # Individual subworkflow tests
└── .nf-test/                      # nf-test working directory
```

### 1. Key Configuration Files

#### `nf-test.config` (Main Configuration)

**Purpose**: Controls nf-test global settings

```groovy
config {
    testsDir "."
    workDir System.getenv("NFT_WORKDIR") ?: ".nf-test"
    configFile "tests/nextflow.config"
    ignore 'modules/nf-core/**/*', 'subworkflows/nf-core/**/*'
    profile "test"
    triggers 'nextflow.config', 'nf-test.config', 'conf/test.config', 'tests/nextflow.config', 'tests/.nftignore'
    plugins {
        load "nft-utils@0.0.3"
        // Add plugins based on your pipeline needs:
        // load "nft-bam@0.6.0"    // For BAM file testing
        // load "nft-vcf@1.0.7"    // For VCF file testing
        // load "nft-fastq@0.0.1"  // For FASTQ file testing
    }
}
```

> **For complete configuration options**, see the [official nf-test configuration documentation](https://www.nf-test.com/docs/configuration/).

### Available nf-test Plugins

For nf-core pipeline testing, plugins provide additional functionality to help validate different output file types and ensure data integrity during testing.

Plugins can be customized - like nft-utils, which has been tailored by the nf-core community for automating the capture and validation of pipeline-level outputs.

| Plugin Name      | Description/Use                                                                                                              |
| ---------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| **nft-utils**    | Essential utility functions like `getAllFilesFromDir()` and `removeNextflowVersion()` - widely used across nf-core pipelines |
| **nft-bam**      | Helper functionality for handling SAM/BAM/CRAM files during tests - critical for genomics pipelines                          |
| **nft-vcf**      | Support for VCF files - utilities for testing variant call format files                                                      |
| **nft-fasta**    | Support for FASTA files - enables validation and testing of FASTA file formats                                               |
| **nft-fastq**    | Support for FASTQ files - validation and testing utilities for sequencing data files                                         |
| **nft-csv**      | Support for CSV files - utilities for testing comma-separated value files                                                    |
| **nft-compress** | Support for ZIP files - enables testing of compressed file handling                                                          |
| **nft-anndata**  | Support for AnnData (h5ad) files - utilities for testing single-cell analysis data formats                                   |
| **nft-tiff**     | Support for TIFF files - enables testing of Tagged Image File Format files                                                   |

> **Note:** For the complete list of available plugins and their latest versions, visit the [nf-test plugins registry](https://plugins.nf-test.com/).

#### `tests/nextflow.config` (Test-Specific Settings)

**Purpose**: Defines parameters and settings specifically for test execution
This is very similar to your typical nf-core test profiles.
However one change is the `aws.client.anonymous` option that <XYZ>

```groovy
params {
    modules_testdata_base_path   = 's3://ngi-igenomes/testdata/nf-core/modules/'
    pipelines_testdata_base_path = 'https://raw.githubusercontent.com/nf-core/test-datasets/PIPELINE_NAME/'
    input       = "${projectDir}/assets/samplesheet.csv"
    outdir      = 'results'
}

process {
    resourceLimits = [
        cpus: 2,
        memory: '3.GB',
        time: '2.h'
    ]
}

aws.client.anonymous = true
```

## Test Data Integration

nf-core tools 3.3+ includes a command for discovering test datasets, making it easier to refer to these within your tests:

```bash
# List all available test datasets
nf-core test-datasets list

# List datasets for a specific branch
nf-core test-datasets list --branch mag

# Search for specific datasets
nf-core test-datasets search --branch mag minigut_reads

# Generate download URLs for datasets
nf-core test-datasets list --branch mag --generate-dl-url

# Generate nf-test compatible paths
nf-core test-datasets list --branch mag --generate-nf-path
```

> **Note:** This feature requires [nf-core/tools 3.3+](https://nf-co.re/blog/2025/tools-3_3#new-nf-core-test-datasets-command).

## Generating Tests (Optional)

Most nf-core pipelines already have tests generated. If needed:

```bash
# Generate test template for a process
nf-test generate process path/to/main.nf

# Generate test for workflow
nf-test generate workflow path/to/main.nf
```

## Next Steps

With your project properly configured, proceed to [Testing Modules](./04_testing_modules.md) to start writing comprehensive module tests.
