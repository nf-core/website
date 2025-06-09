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

## Project Initialization Steps

### 1. Initialize nf-test (if not present. Most nf-core pipelines already have nf-test initialized.)

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

## Best Practices for nf-core Pipeline Testing

1. **Use remote test data**: Links to GitHub test-datasets for consistency across environments
2. **Resource limits**: Always set sensible limits for CI/CD environments
3. **Anonymous AWS access**: Essential for public repositories accessing S3 data
4. **Minimal datasets**: Keep test data small but representative of real use cases
5. **Mixed test scenarios**: Include appropriate test cases for your pipeline's input types
6. **Plugin management**: Use specific plugin versions for reproducibility
7. **Ignore patterns**: Properly configure .nftignore for stable snapshots
8. **Profile separation**: Use different profiles for different testing scenarios

For comprehensive information about nf-core test datasets and how to use them in your tests, see the [nf-core test-datasets documentation](https://nf-co.re/docs/tutorials/tests_and_test_data/test_data).


## Next Steps

With your project properly configured using standard nf-core patterns, proceed to [Testing Modules](./04_testing_modules.md) to start writing comprehensive module tests.