---
title: "8. Configuration Management"
subtitle: Managing test configurations and parameters
weight: 80
---

## Configuration Hierarchy

nf-test follows a hierarchical configuration system:

1. Global `.nf-test.config`
2. Test-specific `nextflow.config` in test directories
3. Inline configuration in test files
4. Process-specific configurations

## Global Configuration

### .nf-test.config

```groovy
config {
    // Test execution settings
    testsDir = "tests"
    workDir = ".nf-test"

    // List of filenames or patterns that should trigger a full test run
    triggers 'nextflow.config', 'nf-test.config', 'conf/test.config', 'conf/test_full.config'

    // Default parameters for all tests
    params {
        modules_testdata_base_path = 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/'
        pipelines_testdata_base_path = 'https://raw.githubusercontent.com/nf-core/test-datasets/nf-core-pipeline/testdata/'

        // Resource defaults
        max_memory = '6.GB'
        max_cpus = 2
        max_time = '6.h'
    }

    // Profile configurations
    profiles {
        docker {
            docker.enabled = true
            singularity.enabled = false
        }

        singularity {
            singularity.enabled = true
            docker.enabled = false
        }
    }

    // Plugin configurations
    plugins {
        load "nft-utils@0.0.3"
    }
}
```

## Trigger Configuration

The `triggers` directive specifies files that should cause a full test run when changed:

```groovy
config {
    // Standard nf-core triggers
    triggers 'nextflow.config', 'nf-test.config', 'conf/test.config', 'conf/test_full.config'

    // Additional custom triggers
    triggers 'assets/samplesheet_check.py', 'lib/WorkflowMain.groovy', 'bin/*'
}
```

This is essential for CI/CD systems to determine when comprehensive testing is needed versus targeted module/subworkflow testing.

## Test-Specific Configuration

### nextflow.config in Test Directories

```groovy
// modules/nf-core/fastqc/tests/nextflow.config
params {
    // Test-specific parameters
    publish_dir_mode = 'copy'

    // Override defaults for this module
    max_memory = '2.GB'
    max_cpus = 1
}

process {
    // Module-specific settings
    withName: "FASTQC" {
        ext.args = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
}

// Container settings
docker.runOptions = '-u $(id -u):$(id -g)'
```

## Inline Configuration

### Process Configuration

```groovy
nextflow_process {
    name "Test Process FASTQC"
    script "../main.nf"
    process "FASTQC"
    config "./nextflow.config"

    test("Custom configuration") {
        when {
            params {
                // Test-specific parameters
                fastqc_args = '--extract --threads 2'
                outdir = './custom_results'
            }
            process {
                """
                input[0] = [
                    [ id:'test', single_end:true ],
                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)
                ]
                """
            }
        }
        then {
            assertAll(
                { assert process.success },
                { assert snapshot(process.out).match() }
            )
        }
    }
}
```

### Dynamic Configuration

```groovy
test("Resource scaling") {
    when {
        params {
            // Dynamically set based on input size
            max_cpus = file('large_input.fastq.gz').size() > 1000000 ? 4 : 2
            max_memory = file('large_input.fastq.gz').size() > 1000000 ? '8.GB' : '4.GB'
        }
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('large_input.fastq.gz', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            { assert snapshot(process.out).match() }
        )
    }
}
```

## Parameter Management

### Parameter Validation

```groovy
nextflow_process {
    name "Test Process with Validation"
    script "../main.nf"
    process "TOOL"

    test("Valid parameters") {
        when {
            params {
                min_quality = 20
                max_length = 500
            }
            process {
                """
                // Validate parameters before use
                assert params.min_quality >= 0 && params.min_quality <= 40
                assert params.max_length > 0

                input[0] = [
                    [ id:'test', min_qual: params.min_quality, max_len: params.max_length ],
                    file('input.fastq.gz', checkIfExists: true)
                ]
                """
            }
        }
        then {
            assertAll(
                { assert process.success },
                { assert snapshot(process.out).match() }
            )
        }
    }
}
```

### Conditional Parameters

```groovy
test("Conditional configuration") {
    when {
        params {
            mode = 'advanced'
            // Set parameters based on mode
            tool_args = params.mode == 'advanced' ? '--comprehensive --detailed' : '--basic'
            memory_multiplier = params.mode == 'advanced' ? 2 : 1
        }
        process {
            """
            input[0] = [
                [ id:'test', mode: params.mode ],
                file('input.txt', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            { assert snapshot(process.out).match() }
        )
    }
}
```

## Environment-Specific Configuration

### Container Configuration

```groovy
// nextflow.config
profiles {
    docker {
        docker.enabled = true
        docker.userEmulation = true
        process.container = 'biocontainers/fastqc:0.11.9--0'
    }

    singularity {
        singularity.enabled = true
        singularity.autoMounts = true
        process.container = 'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0'
    }

    conda {
        conda.enabled = true
        process.conda = 'bioconda::fastqc=0.11.9'
    }
}
```

### Resource Profiles

```groovy
profiles {
    test {
        params.max_memory = '6.GB'
        params.max_cpus = 2
        params.max_time = '6.h'
    }

    test_full {
        params.max_memory = '128.GB'
        params.max_cpus = 16
        params.max_time = '72.h'
    }

    github_actions {
        params.max_memory = '6.GB'
        params.max_cpus = 2
        params.max_time = '6.h'

        // CI-specific settings
        docker.enabled = true
        singularity.enabled = false
    }
}
```

## Advanced Configuration Patterns

### Shared Configuration

```groovy
// shared/test.config
params {
    // Common test parameters
    outdir = './results'
    publish_dir_mode = 'copy'

    // Common test data paths
    test_data_base = 'https://raw.githubusercontent.com/nf-core/test-datasets/'
}

process {
    // Common process settings
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

Include in test configs:

```groovy
// modules/nf-core/tool/tests/nextflow.config
includeConfig '../../../shared/test.config'

// Module-specific overrides
process {
    withName: "TOOL" {
        ext.args = '--specific-args'
    }
}
```

### Configuration Functions

```groovy
// Helper functions for configuration
def getResourceConfig(inputSize) {
    if (inputSize > 1000000000) { // 1GB
        return [
            cpus: 4,
            memory: '8.GB',
            time: '4.h'
        ]
    } else {
        return [
            cpus: 2,
            memory: '4.GB',
            time: '2.h'
        ]
    }
}

test("Dynamic resource allocation") {
    when {
        params {
            def inputFile = file('input.bam', checkIfExists: true)
            def resources = getResourceConfig(inputFile.size())

            max_cpus = resources.cpus
            max_memory = resources.memory
            max_time = resources.time
        }
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('input.bam', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            { assert snapshot(process.out).match() }
        )
    }
}
```

## Configuration Best Practices

### 1. Use Hierarchical Configuration

- Global defaults in `.nf-test.config`
- Component-specific settings in local `nextflow.config`
- Test-specific overrides inline

### 2. Environment Separation

```groovy
profiles {
    test { /* minimal resources */ }
    test_full { /* full resources */ }
    ci { /* CI-specific settings */ }
    local { /* local development */ }
}
```

### 3. Parameter Documentation

```groovy
// Document parameters in config files
params {
    // Input/Output
    input = null        // Path to input samplesheet
    outdir = './results' // Output directory

    // Analysis options
    skip_fastqc = false  // Skip FastQC steps
    skip_multiqc = false // Skip MultiQC report

    // Resource limits
    max_memory = '128.GB' // Maximum memory allocation
    max_cpus = 16        // Maximum CPU cores
    max_time = '240.h'   // Maximum execution time
}
```

## Next Steps

Continue to [Snapshot Management](./09_snapshot_management.md) to learn about working with test snapshots.
