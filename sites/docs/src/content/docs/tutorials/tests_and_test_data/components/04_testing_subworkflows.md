---
title: "4. Testing Subworkflows"
subtitle: Testing nf-core subworkflows
weight: 40
---

## Workflow Testing with nf-test

nf-test allows you to test specific workflows defined in a workflow file. Subworkflows combine multiple modules and require comprehensive testing strategies to ensure proper integration:

```mermaid
flowchart TD
    A[Subworkflow Testing] --> B[Subworkflow Components]

    B --> B1[Multiple Modules]
    B --> B2[Shared Parameters]
    B --> B3[Output Integration]

    B1 --> C[Test Strategy]
    B2 --> C
    B3 --> C

    C --> D{Dependencies Required?}

    D --> |Yes| E[Setup Block]
    D --> |No| F[Direct Testing]

    E --> E1[Run dependency processes]
    E1 --> E2[Capture outputs]
    E2 --> E3[Pass to main workflow]
    E3 --> G[Test Execution]

    F --> G

    G --> H[Multi-Module Validation]
    H --> H1[All modules succeed]
    H2[Parameter passing]
    H3[Output channel integrity]
    H4[Resource usage]
    H --> H2
    H --> H3
    H --> H4

    H1 --> I[Comprehensive Assertions]
    H2 --> I
    H3 --> I
    H4 --> I

    I --> I1[File-based snapshots]
    I --> I2[Content verification]
    I --> I3[Channel structure]
    I --> I4[Version tracking]

    I1 --> J[Tag-based Organization]
    I2 --> J
    I3 --> J
    I4 --> J

    J --> J1[subworkflows tag]
    J --> J2[Individual module tags]
    J --> J3[Algorithm tags]

    J1 --> K["nf-core subworkflows test"]
    J2 --> K
    J3 --> K

    K --> L{Results}
    L --> |Pass| M[✅ Subworkflow Ready]
    L --> |Fail| N[Debug Integration Issues]
    N --> O[Check module compatibility]
    O --> P[Verify parameter flow]
    P --> Q[Review output channels]
    Q --> G
```

The basic syntax for a workflow test follows this structure:

```groovy
nextflow_workflow {
    name "<NAME>"
    script "<PATH/TO/NEXTFLOW_SCRIPT.nf>"
    workflow "<WORKFLOW_NAME>"

    test("<TEST_NAME>") {
        // Test implementation
    }
}
```

**Key Points:**

- Script paths starting with `./` or `../` are relative to the test script location
- Use relative paths to reference files within the same directory or parent directories

### Essential Assertions

Workflow tests commonly use these assertions:

```groovy
// Workflow status
assert workflow.success
assert workflow.failed
assert workflow.exitStatus == 0

// Error handling
assert workflow.errorReport.contains("....")

// Trace analysis
assert workflow.trace.succeeded().size() == 3  // succeeded tasks
assert workflow.trace.failed().size() == 0     // failed tasks
assert workflow.trace.tasks().size() == 3      // all tasks

// Output validation
assert workflow.stdout.contains("Hello World") == 3
```

## Philosophy of nf-test for nf-core Subworkflows

Following the [nf-core testing guidelines](https://nf-co.re/docs/tutorials/tests_and_test_data/nf-test_writing_tests), each nf-core subworkflow should include comprehensive tests that:

- Each subworkflow contains a `tests/` folder beside the `main.nf` of the subworkflow itself
- Test files come with snapshots of subworkflow output channels
- Tests verify both functionality and expected outputs of all included modules
- Support testing with different parameter combinations
- Include proper setup blocks for complex dependencies

## 1. Creating a New Subworkflow with Tests

Creating a new subworkflow automatically creates a test file based on the template.

```bash
# Create a new subworkflow using nf-core tools
cd path/to/subworkflows
nf-core subworkflows create fastq_align_qc

# This creates the subworkflow structure:
# subworkflows/nf-core/fastq_align_qc/
# ├── main.nf
# ├── meta.yml
# └── tests/
#     ├── main.nf.test
#     ├── nextflow.config
#     └── tags.yml
```

The generated test file (`tests/main.nf.test`) will include comprehensive tagging for all modules in the subworkflow:

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
                input[0] = Channel.of([
                            [ id:'test', single_end:true ],
                            file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)
                ])
                input[1] = Channel.of([
                            [ id:'test' ],
                            file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/genome.fasta', checkIfExists: true)
                ])
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

## 2. Testing Subworkflows with Setup Dependencies

For subworkflows that require setup (like index generation), use setup blocks. Here's an example for a BWA alignment subworkflow:

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

## 4. Testing Parameter Variations

Test different parameter combinations that affect subworkflow behavior:

### Creating Parameter-Specific Configuration

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

## 5. Testing Output Channels Comprehensively

### BAM File Testing with MD5 Checksums

Always use MD5 checksums for BAM files to ensure content consistency:

```groovy
{ assert snapshot(
    workflow.out.bam.collect { meta, bamfile -> bam(bamfile).getReadsMD5() },
    workflow.out.bai.collect { meta, bai -> file(bai).name },
    workflow.out.stats.collect { meta, stats -> file(stats).name },
    workflow.out.flagstat.collect { meta, flagstat -> file(flagstat).name },
    workflow.out.versions
    ).match() }
```

### File Name Testing for Stable Names

For files with stable names but variable content:

```groovy
workflow.out.bai.collect { meta, bai -> file(bai).name },
workflow.out.picard_metrics.collect { meta, metrics -> file(metrics).name },
workflow.out.multiqc.flatten().collect { path -> file(path).name }
```

## 7. Updating Subworkflow Snapshots

When subworkflow outputs change (e.g., due to module version bumps), update snapshots:

```bash
nf-core subworkflows test fastq_align_qc --profile docker --update
```

---

Read more nf-test assertion patterns in the [nf-test assertions examples doc](07_assertions.md)

---

## Next Steps

Continue to [Testing Pipelines](./06_testing_pipelines.md) to learn about end-to-end pipeline testing.
