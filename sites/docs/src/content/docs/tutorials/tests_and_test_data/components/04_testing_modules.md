---
title: "4. Testing Modules"
subtitle: Testing individual nf-core modules with comprehensive examples
weight: 40
---

## Module Test Structure

A typical nf-core module test follows established patterns. Here's the structure commonly used for FastQC:

```groovy
nextflow_process {
    name "Test Process FASTQC"
    script "../main.nf"
    process "FASTQC"

    tag "modules"
    tag "modules_nfcore"
    tag "fastqc"

    test("sarscov2 single-end [fastq]") {
        when {
            process {
                """
                input[0] = Channel.of([
                    [ id: 'test', single_end:true ],
                    [ file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true) ]
                ])
                """
            }
        }

        then {
            assertAll (
                { assert process.success },
                // NOTE The report contains the date inside it, which means that the md5sum is stable per day, but not longer than that. So you can't md5sum it.
                // looks like this: <div id="header_filename">Mon 2 Oct 2023<br/>test.gz</div>
                // https://github.com/nf-core/modules/pull/3903#issuecomment-1743620039
                { assert process.out.html[0][1] ==~ ".*/test_fastqc.html" },
                { assert process.out.zip[0][1] ==~ ".*/test_fastqc.zip" },
                { assert path(process.out.html[0][1]).text.contains("<tr><td>File type</td><td>Conventional base calls</td></tr>") },
                { assert snapshot(process.out.versions).match() }
            )
        }
    }
}
```

## Essential Module Testing Components

### 1. Standard Test Tags

All nf-core modules use consistent tagging:

```groovy
tag "modules"
tag "modules_nfcore"
tag "fastqc"  // Tool-specific tag
```

### 2. Meta Map Structure

Meta maps follow nf-core conventions:

```groovy
// Single-end reads
[ id: 'test', single_end: true ]

// Paired-end reads  
[ id: 'test', single_end: false ]
```

### 3. Test Data Usage

Use nf-core test datasets with `checkIfExists: true`:

```groovy
file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)
```

## Real Module Testing Examples

### FastQC Module - Multiple Input Types

```groovy
nextflow_process {
    name "Test Process FASTQC"
    script "../main.nf"
    process "FASTQC"
    
    tag "modules"
    tag "modules_nfcore"
    tag "fastqc"

    test("sarscov2 paired-end [fastq]") {
        when {
            process {
                """
                input[0] = Channel.of([
                    [id: 'test', single_end: false], // meta map
                    [ file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true),
                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_2.fastq.gz', checkIfExists: true) ]
                ])
                """
            }
        }

        then {
            assertAll (
                { assert process.success },
                { assert process.out.html[0][1][0] ==~ ".*/test_1_fastqc.html" },
                { assert process.out.html[0][1][1] ==~ ".*/test_2_fastqc.html" },
                { assert process.out.zip[0][1][0] ==~ ".*/test_1_fastqc.zip" },
                { assert process.out.zip[0][1][1] ==~ ".*/test_2_fastqc.zip" },
                { assert path(process.out.html[0][1][0]).text.contains("<tr><td>File type</td><td>Conventional base calls</td></tr>") },
                { assert path(process.out.html[0][1][1]).text.contains("<tr><td>File type</td><td>Conventional base calls</td></tr>") },
                { assert snapshot(process.out.versions).match() }
            )
        }
    }

    test("sarscov2 paired-end [bam]") {
        when {
            process {
                """
                input[0] = Channel.of([
                    [id: 'test', single_end: false], // meta map
                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam', checkIfExists: true)
                ])
                """
            }
        }

        then {
            assertAll (
                { assert process.success },
                { assert process.out.html[0][1] ==~ ".*/test_fastqc.html" },
                { assert process.out.zip[0][1] ==~ ".*/test_fastqc.zip" },
                { assert path(process.out.html[0][1]).text.contains("<tr><td>File type</td><td>Conventional base calls</td></tr>") },
                { assert snapshot(process.out.versions).match() }
            )
        }
    }
}
```

### Bismark Align - Complex Module with Setup

```groovy
nextflow_process {
    name "Test Process BISMARK_ALIGN"
    script "../main.nf"
    process "BISMARK_ALIGN"
    config './nextflow.config'

    tag "modules"
    tag "modules_nfcore"
    tag "bismark"
    tag "bismark/align"
    tag "bismark/genomepreparation"

    setup {
        run("BISMARK_GENOMEPREPARATION") {
            script "../../genomepreparation/main.nf"
            process {
                """
                input[0] = Channel.of([
                            [ id:'test' ], // meta map
                            file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/genome.fasta', checkIfExists: true)
                ])
                """
            }
        }
    }

    test("bowtie2 | single-end | sarscov2 genome [fasta]") {
        when {
            params {
                bismark_args = '--bowtie2'
            }
            process {
                """
                input[0] = Channel.of([
                            [ id:'test', single_end:true ],
                            file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test.methylated_1.fastq.gz', checkIfExists: true)
                ])
                input[1] = Channel.of([
                            [ id:'sarscov2'],
                            file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/genome.fasta', checkIfExists: true)
                ])
                input[2] = BISMARK_GENOMEPREPARATION.out.index
                """
            }
        }

        then {
            assertAll(
                { assert process.success },
                { assert snapshot(
                    process.out.bam.collect { meta, bamfile ->
                                                bam(bamfile).getReadsMD5()
                                            },
                    process.out.report.collect { meta, report ->
                                                    file(report).readLines().contains("Number of alignments with a unique best hit from the different alignments:\t5009")
                                                },
                    process.out.unmapped,
                    process.out.versions
                    ).match()
                }
            )
        }
    }
}
```

### BWA-meth Align - Methylation-Specific Testing

```groovy
nextflow_process {
    name "Test Process BWAMETH_ALIGN"
    script "../main.nf"
    process "BWAMETH_ALIGN"

    tag "modules"
    tag "modules_nfcore"
    tag "bwameth"
    tag "bwameth/align"
    tag "bwameth/index"

    setup {
        run("BWAMETH_INDEX") {
            script "../../../bwameth/index/main.nf"
            process {
                """
                input[0] = Channel.of([
                    [ id:'test' ], // meta map
                    file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/genome.fasta', checkIfExists: true)
                ])
                """
            }
        }
    }

    test("sarscov2 methylated single_end [fastq]") {
        when {
            process {
                """
                input[0] = Channel.of([
                                [ id:'test', single_end:true ], // meta map
                                file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test.methylated_1.fastq.gz', checkIfExists: true)
				])
                input[1] = Channel.of([
                                [ id:'test' ], // meta map
                                file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/genome.fasta', checkIfExists: true)
                ])
				input[2] = BWAMETH_INDEX.out.index
                """
            }
        }

        then {
            assertAll(
                { assert process.success },
                { assert snapshot(
					bam(process.out.bam[0][1]).getReadsMD5(),
					process.out.versions
					).match()
				}
            )
        }
    }

    test("sarscov2 methylated paired_end [fastq]") {
        when {
            process {
                """
                input[0] = Channel.of([
                                [ id:'test', single_end:false ], // meta map
                                [
                                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test.methylated_1.fastq.gz', checkIfExists: true),
                                    file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test.methylated_2.fastq.gz', checkIfExists: true)
                                ]
				])
				input[1] = Channel.of([
                                [ id:'test' ], // meta map
                                file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/genome.fasta', checkIfExists: true)
                ])
				input[2] = BWAMETH_INDEX.out.index
                """
            }
        }

        then {
            assertAll(
                { assert process.success },
                { assert snapshot(
					bam(process.out.bam[0][1]).getReadsMD5(),
					process.out.versions
					).match()
				}
            )
        }
    }
}
```

## Advanced Module Testing Patterns

### Multiple Input Scenarios - FastQC Multiple Files

```groovy
test("sarscov2 multiple [fastq]") {
    when {
        process {
            """
            input[0] = Channel.of([
                [id: 'test', single_end: false], // meta map
                [ file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true),
                file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_2.fastq.gz', checkIfExists: true),
                file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test2_1.fastq.gz', checkIfExists: true),
                file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test2_2.fastq.gz', checkIfExists: true) ]
            ])
            """
        }
    }

    then {
        assertAll (
            { assert process.success },
            { assert process.out.html[0][1][0] ==~ ".*/test_1_fastqc.html" },
            { assert process.out.html[0][1][1] ==~ ".*/test_2_fastqc.html" },
            { assert process.out.html[0][1][2] ==~ ".*/test_3_fastqc.html" },
            { assert process.out.html[0][1][3] ==~ ".*/test_4_fastqc.html" },
            { assert process.out.zip[0][1][0] ==~ ".*/test_1_fastqc.zip" },
            { assert process.out.zip[0][1][1] ==~ ".*/test_2_fastqc.zip" },
            { assert process.out.zip[0][1][2] ==~ ".*/test_3_fastqc.zip" },
            { assert process.out.zip[0][1][3] ==~ ".*/test_4_fastqc.zip" },
            { assert path(process.out.html[0][1][0]).text.contains("<tr><td>File type</td><td>Conventional base calls</td></tr>") },
            { assert path(process.out.html[0][1][1]).text.contains("<tr><td>File type</td><td>Conventional base calls</td></tr>") },
            { assert path(process.out.html[0][1][2]).text.contains("<tr><td>File type</td><td>Conventional base calls</td></tr>") },
            { assert path(process.out.html[0][1][3]).text.contains("<tr><td>File type</td><td>Conventional base calls</td></tr>") },
            { assert snapshot(process.out.versions).match() }
        )
    }
}
```

### Custom Prefix Testing

```groovy
test("sarscov2 custom_prefix") {
    when {
        process {
            """
            input[0] = Channel.of([
                [ id:'mysample', single_end:true ], // meta map
                file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)
            ])
            """
        }
    }

    then {
        assertAll (
            { assert process.success },
            { assert process.out.html[0][1] ==~ ".*/mysample_fastqc.html" },
            { assert process.out.zip[0][1] ==~ ".*/mysample_fastqc.zip" },
            { assert path(process.out.html[0][1]).text.contains("<tr><td>File type</td><td>Conventional base calls</td></tr>") },
            { assert snapshot(process.out.versions).match() }
        )
    }
}
```

## Testing with Configuration

### Module-Specific Configuration

Create `tests/nextflow.config` for module-specific parameters:

```groovy
process {
    withName: 'METHYLDACKEL_EXTRACT' {
        ext.args = '--methylKit'
    }
}
```

### Configuration-Based Test

```groovy
test("test-methyldackel-extract-methyl-kit") {
    config "./nextflow.config"

    when {
        process {
            """
            input[0] = [
                [ id:'test', single_end:false ], // meta map
                file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/bam/test.paired_end.methylated.sorted.bam', checkIfExists: true),
                file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/bam/test.paired_end.methylated.sorted.bam.bai', checkIfExists: true)
            ]
            input[1] = file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/genome.fasta', checkIfExists: true)
            input[2] = file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/genome.fasta.fai', checkIfExists: true)
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

## Stub Testing for Modules

Stub tests verify module structure without running actual tools:

```groovy
test("sarscov2 single-end [fastq] - stub") {
    options "-stub"
    
    when {
        process {
            """
            input[0] = Channel.of([
                [ id: 'test', single_end:true ],
                [ file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true) ]
            ])
            """
        }
    }

    then {
        assertAll (
            { assert process.success },
            { assert snapshot(process.out).match() }
        )
    }
}
```

## nf-test Best Practices for Modules

### 1. Use BAM MD5 Checksums

For BAM files, use `bam().getReadsMD5()` for stable testing:

```groovy
{ assert snapshot(
    bam(process.out.bam[0][1]).getReadsMD5(),
    process.out.versions
    ).match()
}
```

### 2. Content Validation

Check specific content in output files:

```groovy
{ assert path(process.out.html[0][1]).text.contains("<tr><td>File type</td><td>Conventional base calls</td></tr>") }
```

### 3. Report Validation

For tool reports, validate specific metrics:

```groovy
{ assert file(report).readLines().contains("Number of alignments with a unique best hit from the different alignments:\t5009") }
```

### 4. File Pattern Matching

Use regex patterns to verify output file names:

```groovy
{ assert process.out.html[0][1] ==~ ".*/test_fastqc.html" }
{ assert process.out.zip[0][1] ==~ ".*/test_fastqc.zip" }
```

### 5. Comprehensive Testing

Test multiple scenarios systematically:

- Single-end vs paired-end reads
- Different input formats (FASTQ vs BAM)
- Multiple files
- Custom prefixes
- Edge cases
- Stub mode

### 6. Parameter Testing

Test different tool parameters:

```groovy
when {
    params {
        bismark_args = '--bowtie2'  // vs '--hisat2'
    }
    // ... test implementation
}
```

## Common Module Testing Patterns

### Multi-Input Modules

```groovy
// Input with meta, files, and reference
input[0] = Channel.of([meta, files])
input[1] = Channel.of([meta, reference])  
input[2] = index_channel
```

### Optional Inputs

```groovy
// With optional BAM index
input[0] = [
    [ id:'test', single_end:false ],
    file('test.bam', checkIfExists: true),
    file('test.bam.bai', checkIfExists: true)  // Optional index
]
```

### Module Chains

Use `setup` blocks for dependent modules:

```groovy
setup {
    run("MODULE_INDEX") {
        script "../../index/main.nf"
        process {
            """
            input[0] = reference_channel
            """
        }
    }
}
```

## Next Steps

Continue to [Testing Subworkflows](./05_testing_subworkflows.md) to learn about testing more complex components. 