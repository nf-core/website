---
title: "Repository Setup"
subtitle: Configuring your nf-core pipeline repository for testing with nf-test
weight: 20
---

## nf-test repository architecture

Understanding how nf-test integrates with your nf-core pipeline repository is crucial for effective testing.
The following diagram shows the key components and their relationships:

## Example nf-core repository structure

The following example shows an nf-core pipeline project with nf-test:

```
my-pipeline/
├── main.nf                          # Main pipeline script
├── nextflow.config                  # Main pipeline configuration
├── nextflow_schema.json            # Pipeline parameter schema
├── nf-test.config                  # nf-test configuration
├── modules.json                    # Module dependencies
├── tests/                          # Pipeline-level tests
│   ├── default.nf.test            # Main pipeline test
│   ├── test_alternatives.nf.test           # Pipeline test alternative tools
|   ├── test_minimal.nf.test           # Pipeline test most basic run
│   ├── nextflow.config            # Test-specific config
│   └── .nftignore                 # Files to ignore in snapshots
├── conf/                           # Configuration profiles
│   ├── test.config                # Test profile
│   ├── test_full.config           # Full test profile
│   ├── base.config                # Base configuration
│   ├── test_alternatives.config                # Parameter configuration for test_alternative test
│   ├── test_minimal.config                # Parameter configuration for test_minimal test
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

## nf-test files overview

The following files are required for nf-test configuration:

**Core configuration files**

You'll edit these:

- `nf-test.config` - Main nf-test configuration (plugins, profiles, global settings)
- `tests/nextflow.config` - Test-specific Nextflow parameters and settings
- `tests/.nftignore` - Files to exclude from nf-test snapshots

**Test files**

You'll create/edit these:

- `tests/*.nf.test` - Pipeline-level test files (e.g., `default.nf.test` and `test_minimal.nf.test`)
- `modules/*/tests/main.nf.test` - Individual module test files
- `subworkflows/*/tests/main.nf.test` - Individual subworkflow test files

**Generated files**

nf-test creates these automatically:

- `tests/*.nf.test.snap` - Snapshot files containing expected test outputs
- `.nf-test/` - Working directory for nf-test (i.e., cache, temporary files)

### Key configuration files

#### `nf-test.config` (main configuration)

**Purpose**: Controls nf-test global settings and behavior across your entire pipeline

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
        load "nft-bam@0.6.0"    // For BAM file testing
        load "nft-vcf@1.0.7"    // For VCF file testing
        load "nft-fastq@0.0.1"  // For FASTQ file testing
    }
}
```

**Key configuration options**

nf-core specific:

- `profile "test"`: Sets `test` as the default Nextflow profile for all tests. This means every test will automatically use your pipeline's test configuration by default (resource limits, test data paths, etc.) without needing to specify `-profile test` each time.

- `ignore`: Excludes nf-core modules and subworkflows from testing in your pipeline repository. These components have their own tests in the nf-core/modules repository. For local/custom components, consider adding: `'modules/local/**/tests/*', 'subworkflows/local/**/tests/*'`

- `triggers`: Files that automatically invalidate test caches when changed. Critical for nf-core pipelines because changes to configuration files should trigger test re-runs to ensure consistency.

- `configFile "tests/nextflow.config"`: Points to your test-specific Nextflow configuration. This allows separation of test settings from production pipeline settings.

- `testsDir "."`: Sets the root directory for test discovery. Using "." means nf-test will recursively find all `*.nf.test` files in your repository (but ignore what is defined in the `ignore` option above).

- `workDir`: Defines where nf-test stores its working files and caches. The environment variable `NFT_WORKDIR` allows CI systems to customize this location.

- `plugins`: Essential for nf-core testing. The `nft-utils` plugin provides helper functions like `getAllFilesFromDir()` that are used across many nf-core pipeline tests.

:::note
For complete configuration options, see the [official nf-test configuration documentation](https://www.nf-test.com/docs/configuration/).
:::

### Available nf-test plugins

For nf-core pipeline testing, plugins provide additional functionality to help validate different output file types and ensure data integrity during testing.

Plugins can be customized, like nft-utils, which the nf-core community has tailored to automate the capture and validation of pipeline-level outputs.

| Plugin name      | Description/Use                                                                                                              |
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

:::note
For the complete list of available plugins and their latest versions, visit the [nf-test plugins registry](https://plugins.nf-test.com/).
:::

#### `tests/nextflow.config`

Defines base parameters and settings across all test executions
This is very similar to your typical nf-core test profiles.
However one change is the `aws.client.anonymous` option that ensures anonymous access to AWS S3 buckets for test data, which is useful when your test data is publicly accessible and you don't want to configure AWS credentials for CI/CD.

:::note
Base configuration like below can be overridden for other tests in specific test\_\*.config files
:::

```groovy
params {
    modules_testdata_base_path   = 's3://ngi-igenomes/testdata/nf-core/modules/'
    pipelines_testdata_base_path = 'https://raw.githubusercontent.com/nf-core/test-datasets/<PIPELINE_NAME>/' // Replace <PIPELINE_NAME> with your pipeline's name
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

## Test data integration

From nf-core tools v3.3 onwards, there is a command for discovering test datasets, making it easier to refer to these within your tests:

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

:::note
For more information, see the nf-core/tools v3.3 [blog post](https://nf-co.re/blog/2025/tools-3_3#new-nf-core-test-datasets-command).
:::

## Pipeline-level testing setup

### nf-core/tools 3.3+ pipeline template

Starting with nf-core/tools v3.3, pipeline-level nf-tests are included in the nf-core pipeline template.
The template provides:

- `tests/default.nf.test`: A default pipeline test that mirrors the setup in `config/test.config`
- `tests/nextflow.config`: Base test-specific Nextflow configuration
- `tests/.nftignore`: Files to ignore in nf-test snapshots
- `nf-test.config`: Pipeline-level nf-test configuration

### Initial snapshot generation ⚠️

:::warning
The pipeline template includes `tests/default.nf.test` but **does not include a snapshot file**.
:::

This means:

1. **The initial nf-test CI run will fail** because there's no snapshot for `default.nf.test`
2. **You must generate the snapshot manually** before the tests will pass

To generate the required snapshot:

```bash
# Generate snapshot for the default test (run twice to ensure stability)
nf-test test tests/ --profile=+docker
nf-test test tests/ --profile=+docker

# Or with other execution profiles:
nf-test test tests/ --profile=+singularity
nf-test test tests/ --profile=+singularity

nf-test test tests/ --profile=+conda
nf-test test tests/ --profile=+conda
```

:::Note

- The `=+` notation extends the Nextflow `-profile test` option rather than overwriting it.
- Running the test **twice is required** to ensure snapshot stability: the first run generates the snapshot, the second run validates it's consistent.
  :::

After running this command:

1. **Commit the generated snapshot**: `tests/default.nf.test.snap`
2. **Push to your repository** to fix the failing CI

### Creating additional pipeline tests

You can create additional pipeline tests by copying and renaming the default test:

```bash
# Copy the default test for different scenarios
cp tests/default.nf.test tests/test_alternative.nf.test
cp tests/default.nf.test tests/test_full.nf.test
cp tests/default.nf.test tests/test_minimal.nf.test
```

Generally you should have one `conf/test_*.config` per `test_*.nf.test` file.

Each test config should include the parameter and/or input dataset combination you want to test.

You then modify each `test_*.nf.test` file to describe what you want to test in the output of the pipeline.

Finally, each test will need its own snapshot generated using the same `nf-test test` command.

## Generating tests (optional)

For (local) module and subworkflow tests, you can generate test templates:

```bash
# Generate test template for a process
nf-test generate process path/to/main.nf

# Generate test for workflow
nf-test generate workflow path/to/main.nf
```

## Next steps

With your project properly configured, proceed to [Testing Modules](./03_testing_modules.md) to start writing comprehensive module tests.
