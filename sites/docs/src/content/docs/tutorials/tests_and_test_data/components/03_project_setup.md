---
title: "3. Project Setup"
subtitle: Configuring your nf-core pipeline project for testing with nf-test
weight: 30
---

## nf-core Project Structure

Here's the standard structure of an nf-core pipeline project with nf-test:

```
my-pipeline/
├── main.nf                          # Main pipeline script
├── nextflow.config                  # Main pipeline configuration
├── nextflow_schema.json            # Pipeline parameter schema
├── nf-test.config                  # nf-test configuration
├── modules.json                    # Module dependencies
├── tests/                          # Pipeline-level tests
│   ├── default.nf.test            # Main pipeline test
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

## Real nf-test Configuration

### Standard nf-test.config for nf-core pipelines

```groovy
config {
    // location for all nf-test tests
    testsDir "."

    // nf-test directory including temporary files for each test
    workDir System.getenv("NFT_WORKDIR") ?: ".nf-test"

    // location of an optional nextflow.config file specific for executing tests
    configFile "tests/nextflow.config"

    // ignore tests coming from the nf-core/modules repo
    ignore 'modules/nf-core/**/*', 'subworkflows/nf-core/**/*'

    // run all test with defined profile(s) from the main nextflow.config
    profile "test"

    // list of filenames or patterns that should trigger a full test run
    triggers 'nextflow.config', 'nf-test.config', 'conf/test.config', 'tests/nextflow.config', 'tests/.nftignore'

    // load commonly used plugins
    plugins {
        load "nft-utils@0.0.3"
        // Add other plugins as needed for your pipeline
        // load "nft-bam@0.6.0"    // For genomics pipelines working with BAM files
        // load "nft-vcf@1.0.7"    // For variant calling pipelines
        // load "nft-fastq@0.0.1"  // For quality control and preprocessing
    }
}
```

### Key Configuration Elements

#### 1. **Test Directory Structure**
- `testsDir "."` - Tests are located throughout the project
- Ignores nf-core modules/subworkflows to test only local components
- Uses environment variable `NFT_WORKDIR` for flexible work directory

#### 2. **Common Plugins**
- **nft-utils@0.0.3**: For file management and utility functions (widely used)
- **nft-bam@0.5.0**: For BAM file testing utilities (when working with genomics data)

#### 3. **Trigger Files**
Files that trigger full test runs when changed (meaning all tests in the project will run, not just affected ones):
- `nextflow.config` - Main pipeline configuration
- `nf-test.config` - Test configuration
- `conf/test.config` - Test profile configuration
- `tests/nextflow.config` - Test-specific configuration
- `tests/.nftignore` - Snapshot ignore patterns

## Test-Specific Configuration

### Standard tests/nextflow.config for nf-core pipelines

```groovy
/*
========================================================================================
    Nextflow config file for running tests
========================================================================================
*/

params {
    // Base directory for nf-core/modules test data
    modules_testdata_base_path = 's3://ngi-igenomes/testdata/nf-core/modules/'
    
    // Base directory for pipeline-specific test data
    pipelines_testdata_base_path = 'https://raw.githubusercontent.com/nf-core/test-datasets/PIPELINE_NAME/'
    
    // Input data - modify these paths for your pipeline's test data
    input  = "${projectDir}/assets/samplesheet.csv"
    outdir = 'results'
    
    // Add pipeline-specific test parameters here
    // Examples for different pipeline types:
    // genome = 'R64-1-1'                    // For genomics pipelines
    // fasta  = 'test_reference.fa'          // Reference genome
    // gtf    = 'test_annotation.gtf'        // Gene annotations
    // design = 'design.csv'                 // Experimental design file
}

// Resource limits for testing
process {
    resourceLimits = [
        cpus: 2,
        memory: '6.GB',
        time: '6.h'
    ]
}

// CI/CD optimizations
aws.client.anonymous = true

// Disable Nextflow reporting options for faster tests
timeline { enabled = false }
report   { enabled = false }
trace    { enabled = false }
dag      { enabled = false }
```

### Key Test Configuration Features

#### 1. **Parameter Management**
- **modules_testdata_base_path**: Centralized nf-core modules test data
- **Pipeline-specific params**: Add your pipeline's required parameters
- **Remote URLs**: Use for consistency across environments

#### 2. **Resource Management**
- **resourceLimits**: Conservative limits suitable for CI/CD environments
- **Standard limits**: `cpus: 2, memory: '6.GB', time: '6.h'` (nf-core defaults)
- **Scalable**: Adjust based on your pipeline's requirements

#### 3. **CI/CD Optimizations**
- **aws.client.anonymous = true**: Enables S3 access in CI environments
- **Disabled reporting**: Removes timeline, report, trace, dag for faster execution
- **Profile-based configuration**: Uses `test` profile for consistent testing

## Profile Configuration

### Test Profile (conf/test.config)

```groovy
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow config file for running minimal tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Defines input files and everything required to run a fast and simple pipeline test.

    Use as follows:
        nextflow run nf-core/PIPELINE_NAME -profile test,<docker/singularity> --outdir <OUTDIR>
----------------------------------------------------------------------------------------
*/

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Input data
    input = "${projectDir}/assets/samplesheet.csv"

    // Genome references
    fasta       = "${params.pipelines_testdata_base_path}/reference/genome.fa.gz"
    fasta_index = "${params.pipelines_testdata_base_path}/reference/genome.fa.fai"
}

process {
    resourceLimits = [
        cpus: 4,
        memory: '15.GB',
        time: '1.h'
    ]
}
```

## Test Data Structure

### Example Samplesheet (assets/samplesheet.csv)

```csv
sample,fastq_1,fastq_2
test_sample1,https://github.com/nf-core/test-datasets/raw/PIPELINE_NAME/testdata/sample1_R1.fastq.gz,https://github.com/nf-core/test-datasets/raw/PIPELINE_NAME/testdata/sample1_R2.fastq.gz
test_sample2,https://github.com/nf-core/test-datasets/raw/PIPELINE_NAME/testdata/sample2.fastq.gz,
```

### Test Samplesheet Best Practices

1. **Mixed data types**: Include both single-end and paired-end samples (if applicable)
2. **Remote URLs**: Use direct links to test datasets for consistency
3. **Minimal size**: Keep datasets small for fast testing
4. **Real data**: Use actual data when possible, synthetic when necessary

## Input Validation Schema

### Schema Structure (assets/schema_input.json)

```json
{
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://raw.githubusercontent.com/nf-core/PIPELINE_NAME/master/assets/schema_input.json",
    "title": "nf-core/PIPELINE_NAME pipeline - params.input schema",
    "description": "Schema for the file provided with params.input",
    "type": "array",
    "items": {
        "type": "object",
        "properties": {
            "sample": {
                "type": "string",
                "pattern": "^\\S+$",
                "errorMessage": "Sample name must be provided and cannot contain spaces",
                "meta": ["id"]
            },
            "fastq_1": {
                "type": "string",
                "format": "file-path",
                "exists": true,
                "pattern": "^\\S+\\.f(ast)?q\\.gz$",
                "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
            },
            "fastq_2": {
                "type": "string",
                "format": "file-path",
                "exists": true,
                "pattern": "^\\S+\\.f(ast)?q\\.gz$",
                "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
            }
        },
        "required": ["sample", "fastq_1"]
    }
}
```

## Project Initialization Steps

### 1. Initialize nf-test (if not present)

```bash
# Initialize nf-test in the project
nf-test init

# This creates:
# - nf-test.config
# - .nf-test/ directory
```

### 2. Set Up Test Configuration

Create or verify these key files using the standard nf-core structure:

#### nf-test.config
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

#### tests/nextflow.config
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
manifest { nextflowVersion = '!>=24.10.2' }
timeline { enabled = false }
report   { enabled = false }
trace    { enabled = false }
dag      { enabled = false }
```

### 3. Create Test Assets

#### Minimal samplesheet (assets/samplesheet.csv)
```csv
sample,fastq_1,fastq_2
test_sample1,https://github.com/nf-core/test-datasets/raw/PIPELINE_NAME/testdata/sample1_R1.fastq.gz,https://github.com/nf-core/test-datasets/raw/PIPELINE_NAME/testdata/sample1_R2.fastq.gz
test_sample2,https://github.com/nf-core/test-datasets/raw/PIPELINE_NAME/testdata/sample2.fastq.gz,
```

### 4. Set Up Profile Configuration

#### conf/test.config
```groovy
params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'
    input = "${projectDir}/assets/samplesheet.csv"
    // Add your specific test parameters
}

process {
    resourceLimits = [
        cpus: 4,
        memory: '15.GB',
        time: '1.h'
    ]
}
```

## Environment Setup

### Development Environment

```bash
# Verify nf-test installation
nf-test --version

# Test basic functionality
nf-test test --tag pipeline --verbose
```

### CI/CD Environment Variables

For GitHub Actions, set these in your repository:

```yaml
env:
  NFT_WORKDIR: /tmp/.nf-test
  NFT_DIFF: "pdiff"
  NFT_DIFF_ARGS: "--line-numbers --expand-tabs=2"
```

## Verification Steps

### 1. Test Configuration

```bash
# Verify nf-test configuration
nf-test list

# Should show all discovered tests
```

### 2. Run Basic Tests

```bash
# Run a simple module test
nf-test test modules/nf-core/fastqc/tests/main.nf.test --profile test,docker

# Run pipeline test
nf-test test tests/default.nf.test --profile test,docker
```

## Best Practices for nf-core Pipeline Testing

1. **Use remote test data**: Links to GitHub test-datasets for consistency across environments
2. **Resource limits**: Always set sensible limits for CI/CD environments
3. **Anonymous AWS access**: Essential for public repositories accessing S3 data
4. **Minimal datasets**: Keep test data small but representative of real use cases
5. **Mixed test scenarios**: Include appropriate test cases for your pipeline's input types
6. **Plugin management**: Use specific plugin versions for reproducibility
7. **Ignore patterns**: Properly configure .nftignore for stable snapshots
8. **Profile separation**: Use different profiles for different testing scenarios

## Next Steps

With your project properly configured using standard nf-core patterns, proceed to [Testing Modules](./04_testing_modules.md) to start writing comprehensive module tests. 