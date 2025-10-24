---
title: "Testing Subworkflows"
subtitle: Comprehensive testing strategies for nf-core subworkflows with nf-test
weight: 40
---

## Overview

The nf-test framework enables comprehensive testing of subworkflows, which combine multiple modules into integrated analysis steps.

This chapter builds upon the concepts introduced in the module testing chapter, covering testing strategies for subworkflows, from basic syntax to complex multi-module integration scenarios.

If you have not already, we highly recommend reading the [module testing](03_testing_modules.md) chapter first.

## Key differences from module testing

Before diving into subworkflow test syntax, it's important to understand how subworkflow testing differs from module testing:

**Workflow vs Process outputs:**

- **Process outputs**: Direct file outputs from a single tool (e.g., `process.out.fastq`, `process.out.versions`)
- **Workflow outputs**: Aggregated outputs from multiple processes with named channels (e.g., `workflow.out.reports`, `workflow.out.multiqc_html`)

**Success checking:**

- **`process.success`**: Checks if a single process completed successfully
- **`workflow.success`**: Checks if the entire workflow (all processes) completed successfully

**Testing scope:**

- **Module tests**: Focus on individual tool functionality, input/output validation, and parameter handling
- **Subworkflow tests**: Focus on integration between modules, data flow, and combined functionality

**Example comparison:**

```groovy
// Module test - single process
assert process.success
assert process.out.fastq.size() == 1

// Subworkflow test - multiple processes
assert workflow.success
assert workflow.trace.succeeded().size() == 3  // 3 processes succeeded
```

## Creating a new subworkflow with tests

When creating a new subworkflow using `nf-core/tools`, a test file is automatically generated based on the template.

```bash
# Create a new subworkflow using nf-core/tools
cd path/to/subworkflows
nf-core subworkflows create fastq_align_qc
```

This creates the following subworkflow structure:

```
subworkflows/nf-core/fastq_align_qc/
├── main.nf
├── meta.yml
└── tests/
    ├── main.nf.test
    └── nextflow.config
```

Once completed, the final test file (`tests/main.nf.test`) will include comprehensive tagging for all modules in the subworkflow:

```groovy
nextflow_workflow {
    name "Test Subworkflow FASTQ_ALIGN_QC"
    script "../main.nf"
    workflow "FASTQ_ALIGN_QC"
    config "./nextflow.config"

    tag "subworkflows"
    tag "subworkflows_nfcore"
    tag "subworkflows/fastq_align_qc"
    tag "fastqc"
    tag "trimgalore"
    tag "bwa/mem"
    tag "samtools/sort"
    tag "samtools/index"
    tag "samtools/stats"
    tag "samtools/flagstat"
    tag "picard/markduplicates"

    test("BWA alignment single-end | default") {
        when {
            workflow {
                """
                // Single-end fastq reads
                input[0] = Channel.of([
                            [ id:'test', single_end:true ],
                            file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)
                ])
                // Reference genome fasta file
                input[1] = Channel.of([
                            [ id:'test' ],
                            file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/genome.fasta', checkIfExists: true)
                ])
                // BWA index
                input[2] = Channel.of([
                            [ id:'test' ],
                            file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/bwa/genome.fasta.{amb,ann,bwt,pac,sa}', checkIfExists: true)
                ])
                """
            }
        }

        then {
            assertAll(
                { assert workflow.success },
                { assert snapshot(workflow.out).match() }
            )
        }
    }
}
```

Run the tests:

```bash
nf-core subworkflows test fastq_align_qc --profile docker
```

## Testing subworkflows with setup dependencies

Like module tests, subworkflows can use `setup` blocks to handle prerequisite steps such as generating a reference genome index. Here's an example for a BWA alignment subworkflow:

```groovy
nextflow_workflow {
    name "Test Subworkflow FASTQ_ALIGN_QC"
    script "../main.nf"
    workflow "FASTQ_ALIGN_QC"
    config "./nextflow.config"

    tag "subworkflows"
    tag "subworkflows_nfcore"
    tag "subworkflows/fastq_align_qc"
    tag "bwa/mem"
    tag "samtools/sort"
    tag "samtools/index"

    setup {
        run("BWA_INDEX") {
            // Path to the module's main.nf script, relative to the subworkflow's tests/ directory
            script "../../../../modules/nf-core/bwa/index/main.nf"
            process {
                """
                input[0] = [
                    [ id:'test' ],
                    file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/genome.fasta', checkIfExists: true)
                ]
                """
            }
        }
    }

    test("BWA alignment single-end | default") {
        when {
            workflow {
                """
                input[0] = Channel.of([
                            [ id:'test', single_end:true ],
                            file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)
                ])
                input[1] = Channel.of([
                            [ id:'test' ],
                            file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/genome.fasta', checkIfExists: true)
                ])
                // Use the output from the setup block
                input[2] = BWA_INDEX.out.index
                """
            }
        }

        then {
            assertAll(
                { assert workflow.success },
                { assert snapshot(
                    workflow.out.bam.collect { meta, bamfile -> bam(bamfile).getReadsMD5() },
                    workflow.out.bai.collect { meta, bai -> file(bai).name },
                    workflow.out.stats.collect { meta, stats -> file(stats).name },
                    workflow.out.flagstat.collect { meta, flagstat -> file(flagstat).name },
                    workflow.out.versions
                    ).match() }
            )
        }
    }
}
```

> `bam(bamfile).getReadsMD5()` is a function from the [`nft-bam`](https://nvnieuwk.github.io/nft-bam/latest/) plugin to get the stable MD5 checksum of the BAM file.

## Testing parameter variations

Test different parameter combinations that could affect subworkflow behavior.

You can do this with Nextflow config files that sit alongside the `main.nf.test` file, or with `params` blocks in the test file itself.

### Creating parameter-specific configuration

Create `tests/nextflow.config` for subworkflow-specific process configuration:

```groovy
params {
    aligner            = "bismark"
    cytosine_report    = false
    skip_deduplication = false
}

process {
    withName: 'BISMARK_ALIGN' {
        ext.args = { params.aligner == 'bismark_hisat' ? ' --hisat2' : ' --bowtie2' }
    }

    withName: 'SAMTOOLS_SORT' {
        ext.prefix = { "${meta.id}.sorted" }
    }
}
```

### Overriding parameters with the `params` block

While a `nextflow.config` file is useful for setting global parameters for all tests of a subworkflow, `nf-test` provides a `params` block within the `when` block to override Nextflow's input `params` for a specific test.

This is useful for testing different parameter combinations without creating multiple config files.

For example, you can override the `aligner` parameter from the `nextflow.config` example to test a different alignment tool:

```groovy
when {
    params {
        aligner = "bismark_hisat"
    }
    workflow {
        """
        // workflow inputs
        """
    }
}
```

## Subworkflow testing principles

Subworkflow testing follows the same core principles as module testing, but adapted for the broader scope of a subworkflow. Each nf-core subworkflow should include comprehensive tests that:

- **Prefer automated MD5 checksums using Snapshots** for output verification when possible, then file content checks, then existence checks as fallbacks
- **Test both regular process and stub modes** to verify functionality and stub outputs
- **Use appropriate test data** from the **nf-core test-datasets** repository
- **Minimal viable tests** that cover the core functionality without excessive complexity
- **Test all different parameter combinations**

:::note
For more nf-test assertion patterns, see the [nf-test assertions examples documentation](./07_assertions.md).
:::

## Next steps

Continue to [Testing Pipelines](./05_testing_pipelines.md) to learn about end-to-end pipeline testing.
