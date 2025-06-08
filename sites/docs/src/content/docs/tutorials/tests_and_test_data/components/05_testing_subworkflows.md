---
title: "5. Testing Subworkflows"
subtitle: Testing subworkflow components with real examples from nf-core/methylseq
weight: 50
---

## Subworkflow Test Structure

Subworkflows combine multiple processes and require comprehensive testing. Here's the structure used in nf-core/methylseq for the Bismark alignment and deduplication subworkflow:

```groovy
nextflow_workflow {
    name "Test Subworkflow FASTQ_ALIGN_DEDUP_BISMARK"
    script "../main.nf"
    workflow "FASTQ_ALIGN_DEDUP_BISMARK"
    config "./nextflow.config"

    tag "subworkflows"
    tag "subworkflows_nfcore"
    tag "subworkflows/fastq_align_dedup_bismark"
    tag "bismark/align"
    tag "samtools/sort"
    tag "samtools/index"
    tag "bismark/deduplicate"
    tag "bismark/methylationextractor"
    tag "bismark/coverage2cytosine"
    tag "bismark/report"
    tag "bismark/summary"
    tag "untar"

    setup {
        run("UNTAR", alias: "BOWTIE2") {
            script "../../../../modules/nf-core/untar/main.nf"
            process {
                """
                input[0] = [
                    [:],
                    file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/Bowtie2_Index.tar.gz', checkIfExists: true)
                ]
                """
            }
        }

        run("UNTAR", alias: "HISAT2") {
            script "../../../../modules/nf-core/untar/main.nf"
            process {
                """
                input[0] = [
                    [:],
                    file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/Hisat2_Index.tar.gz', checkIfExists: true)
                ]
                """
            }
        }
    }

    test("Params: bismark single-end | default") {
        when {
            params {
                aligner            = "bismark"
                cytosine_report    = false
                skip_deduplication = false
            }

            workflow {
                """
                input[0] = Channel.of([
                            [ id:'test', single_end:true ],
                            file('https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub1.fastq.gz', checkIfExists: true)
                ])
                input[1] = Channel.of([
                            [:],
                            file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/genome.fa', checkIfExists: true)
                ])
                input[2] = BOWTIE2.out.untar
                input[3] = params.skip_deduplication
                input[4] = params.cytosine_report
                """
            }
        }

        then {
            assertAll(
                { assert workflow.success },
                { assert snapshot(
                    workflow.out.bam.collect { meta, bamfile -> bam(bamfile).getReadsMD5() },
                    workflow.out.bai.collect { meta, bai -> file(bai).name },
                    workflow.out.coverage2cytosine_coverage,
                    workflow.out.coverage2cytosine_report,
                    workflow.out.coverage2cytosine_summary,
                    workflow.out.methylation_bedgraph,
                    workflow.out.methylation_calls,
                    workflow.out.methylation_coverage,
                    workflow.out.methylation_report,
                    workflow.out.methylation_mbias,
                    workflow.out.bismark_report.collect { meta, report -> file(report).name },
                    workflow.out.bismark_summary[0][1],
                    workflow.out.multiqc.flatten().collect { path -> file(path).name },
                    workflow.out.versions
                    ).match() }
            )
        }
    }
}
```

## Essential Subworkflow Testing Components

### 1. Comprehensive Tagging

Tag all modules that are part of the subworkflow:

```groovy
tag "subworkflows"
tag "subworkflows_nfcore"
tag "subworkflows/fastq_align_dedup_bismark"
tag "bismark/align"
tag "samtools/sort"
tag "samtools/index"
tag "bismark/deduplicate"
tag "bismark/methylationextractor"
tag "bismark/coverage2cytosine"
tag "bismark/report"
tag "bismark/summary"
tag "untar"
```

### 2. Complex Setup Dependencies

Use setup blocks with aliases for multiple dependency modules:

```groovy
setup {
    run("UNTAR", alias: "BOWTIE2") {
        script "../../../../modules/nf-core/untar/main.nf"
        process {
            """
            input[0] = [
                [:],
                file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/Bowtie2_Index.tar.gz', checkIfExists: true)
            ]
            """
        }
    }

    run("UNTAR", alias: "HISAT2") {
        script "../../../../modules/nf-core/untar/main.nf"
        process {
            """
            input[0] = [
                [:],
                file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/Hisat2_Index.tar.gz', checkIfExists: true)
            ]
            """
        }
    }
}
```

### 3. Multiple Output Channel Testing

Test all output channels with appropriate validation:

```groovy
{ assert snapshot(
    workflow.out.bam.collect { meta, bamfile -> bam(bamfile).getReadsMD5() },
    workflow.out.bai.collect { meta, bai -> file(bai).name },
    workflow.out.coverage2cytosine_coverage,
    workflow.out.coverage2cytosine_report,
    workflow.out.coverage2cytosine_summary,
    workflow.out.methylation_bedgraph,
    workflow.out.methylation_calls,
    workflow.out.methylation_coverage,
    workflow.out.methylation_report,
    workflow.out.methylation_mbias,
    workflow.out.bismark_report.collect { meta, report -> file(report).name },
    workflow.out.bismark_summary[0][1],
    workflow.out.multiqc.flatten().collect { path -> file(path).name },
    workflow.out.versions
    ).match() }
```

## Real Subworkflow Testing Examples

### BWA-meth Alignment Subworkflow

```groovy
nextflow_workflow {
    name "Test Subworkflow FASTQ_ALIGN_DEDUP_BWAMETH"
    script "../main.nf"
    workflow "FASTQ_ALIGN_DEDUP_BWAMETH"
    config "./nextflow.config"

    tag "subworkflows"
    tag "subworkflows_nfcore"
    tag "subworkflows/fastq_align_dedup_bwameth"
    tag "bwameth/align"
    tag "parabricks/fq2bammeth"
    tag "samtools/sort"
    tag "samtools/index"
    tag "samtools/flagstat"
    tag "samtools/stats"
    tag "picard/markduplicates"
    tag "samtools/index"
    tag "methyldackel/extract"
    tag "methyldackel/mbias"
    tag "untar"

    setup {
        run("UNTAR") {
            script "../../../../modules/nf-core/untar/main.nf"
            process {
                """
                input[0] = [
                    [:],
                    file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/Bwameth_Index.tar.gz', checkIfExists: true)
                ]
                """
            }
        }
    }

    test("Params: bwameth single-end | default") {
        when {
            workflow {
                """
                input[0] = Channel.of([
                            [ id:'test', single_end:true ],
                            file('https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub1.fastq.gz', checkIfExists: true)
                ])
                input[1] = Channel.of([
                            [:],
                            file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/genome.fa', checkIfExists: true)
                ])
                input[2] = Channel.of([
                            [:],
                            file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/genome.fa.fai', checkIfExists: true)
                ])
                input[3] = UNTAR.out.untar
                input[4] = false // skip_deduplication
                input[5] = false // use_gpu
                """
            }
        }

        then {
            assertAll(
                { assert workflow.success},
                { assert snapshot(
                    workflow.out.bam.collect { meta, bamfile -> bam(bamfile).getReadsMD5() },
                    workflow.out.bai.collect { meta, bai -> file(bai).name },
                    workflow.out.samtools_flagstat,
                    workflow.out.samtools_stats,
                    workflow.out.methydackel_extract_bedgraph,
                    workflow.out.methydackel_extract_methylkit,
                    workflow.out.methydackel_mbias,
                    workflow.out.picard_metrics.collect { meta, metrics -> file(metrics).name },
                    workflow.out.multiqc.flatten().collect { path -> file(path).name },
                    workflow.out.versions
                    ).match() }
            )
        }
    }

    test("Params: bwameth paired-end | default") {
        when {
            workflow {
                """
                input[0] = Channel.of([
                            [ id:'test', single_end:true ],
                            file('https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R1.fastq.gz', checkIfExists: true),
                            file('https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R2.fastq.gz', checkIfExists: true)
                ])
                input[1] = Channel.of([
                            [:],
                            file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/genome.fa', checkIfExists: true)
                ])
                input[2] = Channel.of([
                            [:],
                            file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/genome.fa.fai', checkIfExists: true)
                ])
                input[3] = UNTAR.out.untar
                input[4] = false // skip_deduplication
                input[5] = false // use_gpu
                """
            }
        }

        then {
            assertAll(
                { assert workflow.success},
                { assert snapshot(
                    workflow.out.bam.collect { meta, bamfile -> bam(bamfile).getReadsMD5() },
                    workflow.out.bai.collect { meta, bai -> file(bai).name },
                    workflow.out.samtools_flagstat,
                    workflow.out.samtools_stats,
                    workflow.out.methydackel_extract_bedgraph,
                    workflow.out.methydackel_extract_methylkit,
                    workflow.out.methydackel_mbias,
                    workflow.out.picard_metrics.collect { meta, metrics -> file(metrics).name },
                    workflow.out.multiqc.flatten().collect { path -> file(path).name },
                    workflow.out.versions
                    ).match() }
            )
        }
    }
}
```

## Advanced Subworkflow Testing Patterns

### Parameter-Based Testing

Test different parameter combinations that affect subworkflow behavior:

```groovy
test("Params: bismark paired-end | skip_deduplication") {
    when {
        params {
            aligner            = "bismark"
            skip_deduplication = true
        }

        workflow {
            """
            input[0] = Channel.of([
                        [ id:'test', single_end:true ],
                        file('https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R1.fastq.gz', checkIfExists: true),
                        file('https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R2.fastq.gz', checkIfExists: true)
            ])
            input[1] = Channel.of([
                        [:],
                        file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/genome.fa', checkIfExists: true)
            ])
            input[2] = BOWTIE2.out.untar
            input[3] = params.skip_deduplication
            input[4] = params.cytosine_report
            """
        }
    }

    then {
        assertAll(
            { assert workflow.success },
            { assert snapshot(
                workflow.out.bam.collect { meta, bamfile -> bam(bamfile).getReadsMD5() },
                workflow.out.bai.collect { meta, bai -> file(bai).name },
                workflow.out.coverage2cytosine_coverage,
                workflow.out.coverage2cytosine_report,
                workflow.out.coverage2cytosine_summary,
                workflow.out.methylation_bedgraph,
                workflow.out.methylation_calls,
                workflow.out.methylation_coverage,
                workflow.out.methylation_report,
                workflow.out.methylation_mbias,
                workflow.out.bismark_report.collect { meta, report -> file(report).name },
                workflow.out.bismark_summary[0][1],
                workflow.out.multiqc.flatten().collect { path -> file(path).name },
                workflow.out.versions
                ).match() }
        )
    }
}

test("Params: bismark paired-end | cytosine_report") {
    when {
        params {
            aligner         = "bismark"
            cytosine_report = true
        }

        workflow {
            """
            input[0] = Channel.of([
                        [ id:'test', single_end:true ],
                        file('https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R1.fastq.gz', checkIfExists: true),
                        file('https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R2.fastq.gz', checkIfExists: true)
            ])
            input[1] = Channel.of([
                        [:],
                        file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/genome.fa', checkIfExists: true)
            ])
            input[2] = BOWTIE2.out.untar
            input[3] = params.skip_deduplication
            input[4] = params.cytosine_report
            """
        }
    }

    then {
        assertAll(
            { assert workflow.success },
            { assert snapshot(
                workflow.out.bam.collect { meta, bamfile -> bam(bamfile).getReadsMD5() },
                workflow.out.bai.collect { meta, bai -> file(bai).name },
                workflow.out.coverage2cytosine_coverage,
                workflow.out.coverage2cytosine_report,
                workflow.out.coverage2cytosine_summary,
                workflow.out.methylation_bedgraph,
                workflow.out.methylation_calls,
                workflow.out.methylation_coverage,
                workflow.out.methylation_report,
                workflow.out.methylation_mbias,
                workflow.out.bismark_report.collect { meta, report -> file(report).name },
                workflow.out.bismark_summary[0][1],
                workflow.out.multiqc.flatten().collect { path -> file(path).name },
                workflow.out.versions
                ).match() }
        )
    }
}
```

## Testing with Configuration

### Subworkflow Configuration

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

### Multi-Input Channel Testing

Test complex input channel combinations:

```groovy
workflow {
    """
    input[0] = Channel.of([
                [ id:'test', single_end:true ],
                file('https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub1.fastq.gz', checkIfExists: true)
    ])
    input[1] = Channel.of([
                [:],
                file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/genome.fa', checkIfExists: true)
    ])
    input[2] = Channel.of([
                [:],
                file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/genome.fa.fai', checkIfExists: true)
    ])
    input[3] = UNTAR.out.untar  // From setup block
    input[4] = false            // skip_deduplication
    input[5] = false            // use_gpu
    """
}
```

## Best Practices for Subworkflow Testing

### 1. Comprehensive Output Testing

Test all output channels that the subworkflow produces:

```groovy
{ assert snapshot(
    workflow.out.bam.collect { meta, bamfile -> bam(bamfile).getReadsMD5() },
    workflow.out.bai.collect { meta, bai -> file(bai).name },
    workflow.out.samtools_flagstat,
    workflow.out.samtools_stats,
    workflow.out.methydackel_extract_bedgraph,
    workflow.out.methydackel_extract_methylkit,
    workflow.out.methydackel_mbias,
    workflow.out.picard_metrics.collect { meta, metrics -> file(metrics).name },
    workflow.out.multiqc.flatten().collect { path -> file(path).name },
    workflow.out.versions
    ).match() }
```

### 2. Use Real Test Data

Use domain-specific test data (methylated FASTQ files for methylation workflows):

```groovy
file('https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub1.fastq.gz', checkIfExists: true)
file('https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R1.fastq.gz', checkIfExists: true)
```

### 3. Test Single-End and Paired-End Scenarios

Always test both sequencing modes:

```groovy
test("Params: bwameth single-end | default") {
    when {
        workflow {
            """
            input[0] = Channel.of([
                        [ id:'test', single_end:true ],
                        file('https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub1.fastq.gz', checkIfExists: true)
            ])
            // ... other inputs
            """
        }
    }
    // ... assertions
}

test("Params: bwameth paired-end | default") {
    when {
        workflow {
            """
            input[0] = Channel.of([
                        [ id:'test', single_end:false ],
                        file('https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R1.fastq.gz', checkIfExists: true),
                        file('https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R2.fastq.gz', checkIfExists: true)
            ])
            // ... other inputs
            """
        }
    }
    // ... assertions
}
```

### 4. File Name Stability

For files with stable names but variable content, collect only the filename:

```groovy
workflow.out.bai.collect { meta, bai -> file(bai).name },
workflow.out.picard_metrics.collect { meta, metrics -> file(metrics).name },
workflow.out.multiqc.flatten().collect { path -> file(path).name }
```

### 5. BAM File Testing

Always use MD5 checksums for BAM files:

```groovy
workflow.out.bam.collect { meta, bamfile -> bam(bamfile).getReadsMD5() }
```

### 6. Parameter-Driven Testing

Test critical parameter combinations:

```groovy
params {
    aligner            = "bismark"  // or "bwameth"
    cytosine_report    = false      // or true
    skip_deduplication = false      // or true
}
```

### 7. Flatten MultiQC Outputs

MultiQC outputs may be nested, so flatten them:

```groovy
workflow.out.multiqc.flatten().collect { path -> file(path).name }
```

## Advanced Subworkflow Patterns

### Multiple Reference Indexes

Test with different aligner indexes:

```groovy
setup {
    run("UNTAR", alias: "BOWTIE2") {
        script "../../../../modules/nf-core/untar/main.nf"
        process {
            """
            input[0] = [
                [:],
                file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/Bowtie2_Index.tar.gz', checkIfExists: true)
            ]
            """
        }
    }

    run("UNTAR", alias: "HISAT2") {
        script "../../../../modules/nf-core/untar/main.nf"
        process {
            """
            input[0] = [
                [:],
                file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/Hisat2_Index.tar.gz', checkIfExists: true)
            ]
            """
        }
    }
}
```

### Testing Conditional Outputs

Some outputs may be conditional based on parameters:

```groovy
then {
    assertAll(
        { assert workflow.success },
        {
            if (params.skip_deduplication) {
                // Test outputs when deduplication is skipped
                assert workflow.out.picard_metrics.isEmpty()
            } else {
                // Test outputs when deduplication is performed
                assert !workflow.out.picard_metrics.isEmpty()
            }
        },
        { assert snapshot(workflow.out).match() }
    )
}
```

### Complex Meta Map Propagation

Verify that meta maps are properly propagated through the subworkflow:

```groovy
then {
    assertAll(
        { assert workflow.success },
        {
            workflow.out.bam.each { meta, bam ->
                assert meta.id == 'test'
                assert meta.single_end == true  // or false for paired-end
            }
        },
        { assert snapshot(workflow.out).match() }
    )
}
```

## Common Subworkflow Testing Pitfalls

### 1. Incorrect Channel Structure

Ensure channel structure matches subworkflow expectations:

```groovy
// Correct: Channel with meta and files
Channel.of([
    [ id:'test', single_end:true ],
    file('test.fastq.gz', checkIfExists: true)
])

// Incorrect: Missing meta map
Channel.of(file('test.fastq.gz', checkIfExists: true))
```

### 2. Missing Setup Dependencies

Always include required setup blocks:

```groovy
setup {
    run("INDEX_GENERATION") {
        script "path/to/index/main.nf"
        process {
            """
            input[0] = reference_channel
            """
        }
    }
}
```

### 3. Incomplete Output Testing

Test all relevant outputs, not just success:

```groovy
// Good: Test multiple outputs
{ assert snapshot(
    workflow.out.bam,
    workflow.out.reports,
    workflow.out.metrics,
    workflow.out.versions
).match() }

// Poor: Only test success
{ assert workflow.success }
```

## Next Steps

Continue to [Testing Pipelines](./06_testing_pipelines.md) to learn about end-to-end pipeline testing. 