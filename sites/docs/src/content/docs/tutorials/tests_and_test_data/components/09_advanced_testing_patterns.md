---
title: "10. Advanced Testing Patterns"
subtitle: Complex testing scenarios and patterns
weight: 100
---

## Parameterized Testing

### Matrix Testing

```groovy
nextflow_process {
    name "Test Process TOOL"
    script "../main.nf"
    process "TOOL"
    
    ['single', 'paired'].each { read_type ->
        ['30', '20'].each { quality ->
            test("${read_type}_reads_quality_${quality}") {
                when {
                    params {
                        min_quality = quality as Integer
                        read_mode = read_type
                    }
                    process {
                        """
                        input[0] = [
                            [ id: 'test', single_end: read_type == 'single' ],
                            read_type == 'single' ? 
                                file('single.fastq.gz', checkIfExists: true) :
                                [file('read1.fastq.gz', checkIfExists: true), file('read2.fastq.gz', checkIfExists: true)]
                        ]
                        """
                    }
                }
                then {
                    assertAll(
                        { assert process.success },
                        { assert snapshot(process.out).match("${read_type}_quality_${quality}") }
                    )
                }
            }
        }
    }
}
```

### Data-Driven Testing

```groovy
nextflow_process {
    name "Test Process ALIGNER"
    script "../main.nf"
    process "ALIGNER"
    
    // Test data configurations
    def testConfigs = [
        [
            name: "human_genome",
            genome: "GRCh38",
            expected_alignments: 1000,
            reference: "human_ref.fa"
        ],
        [
            name: "mouse_genome", 
            genome: "GRCm39",
            expected_alignments: 800,
            reference: "mouse_ref.fa"
        ]
    ]
    
    testConfigs.each { config ->
        test("Alignment with ${config.name}") {
            when {
                params {
                    genome = config.genome
                }
                process {
                    """
                    input[0] = [
                        [ id: 'test', genome: '${config.genome}' ],
                        file('reads.fastq.gz', checkIfExists: true)
                    ]
                    input[1] = file('${config.reference}', checkIfExists: true)
                    """
                }
            }
            then {
                assertAll(
                    { assert process.success },
                    { 
                        process.out.bam.each { meta, bam ->
                            def lines = bam.readLines()
                            assert lines.size() >= config.expected_alignments
                        }
                    },
                    { assert snapshot(process.out).match(config.name) }
                )
            }
        }
    }
}
```

## Complex Setup Patterns

### Multi-Step Setup

```groovy
nextflow_workflow {
    name "Test Complex Workflow"
    script "../main.nf"
    workflow "COMPLEX_ANALYSIS"
    
    test("Full analysis pipeline") {
        setup {
            // Step 1: Prepare reference
            run("PREPARE_REFERENCE") {
                script "../../prepare_reference/main.nf"
                workflow {
                    """
                    input[0] = file('reference.fa', checkIfExists: true)
                    """
                }
            }
            
            // Step 2: Index reference
            run("INDEX_REFERENCE") {
                script "../../index_reference/main.nf"
                workflow {
                    """
                    input[0] = PREPARE_REFERENCE.out.processed_ref
                    """
                }
            }
            
            // Step 3: Prepare samples
            run("PREPARE_SAMPLES") {
                script "../../prepare_samples/main.nf"
                workflow {
                    """
                    input[0] = Channel.fromList([
                        [[ id: 'sample1' ], file('sample1.fastq.gz', checkIfExists: true)],
                        [[ id: 'sample2' ], file('sample2.fastq.gz', checkIfExists: true)]
                    ])
                    """
                }
            }
        }
        
        when {
            workflow {
                """
                input[0] = PREPARE_SAMPLES.out.processed_samples
                input[1] = INDEX_REFERENCE.out.indexed_ref
                input[2] = PREPARE_REFERENCE.out.annotation
                """
            }
        }
        
        then {
            assertAll(
                { assert workflow.success },
                { assert workflow.out.results.size() == 2 },
                { assert snapshot(workflow.out).match() }
            )
        }
    }
}
```

### Conditional Setup

```groovy
test("Analysis with optional preprocessing") {
    setup {
        if (params.skip_preprocessing) {
            // Use pre-processed data
            run("LOAD_PREPROCESSED") {
                script "../../load_data/main.nf"
                workflow {
                    """
                    input[0] = file('preprocessed_data.txt', checkIfExists: true)
                    """
                }
            }
        } else {
            // Run preprocessing pipeline
            run("PREPROCESS_DATA") {
                script "../../preprocess/main.nf"
                workflow {
                    """
                    input[0] = file('raw_data.txt', checkIfExists: true)
                    """
                }
            }
        }
    }
    
    when {
        params {
            skip_preprocessing = false
        }
        workflow {
            """
            input[0] = params.skip_preprocessing ? 
                LOAD_PREPROCESSED.out.data : 
                PREPROCESS_DATA.out.processed_data
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
```

## Error and Edge Case Testing

### Expected Failures

```groovy
test("Invalid input format") {
    when {
        process {
            """
            input[0] = [
                [ id: 'test' ],
                file('invalid_format.txt', checkIfExists: true)  // Wrong format
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.failed },
            { assert process.stderr.contains("Invalid file format") }
        )
    }
}
```

### Edge Cases

```groovy
test("Empty input file") {
    when {
        process {
            """
            input[0] = [
                [ id: 'test' ],
                file('empty_file.txt', checkIfExists: true)  // Empty file
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            { assert process.out.size() == 1 },
            // Should handle empty input gracefully
            { assert snapshot(process.out).match("empty_input") }
        )
    }
}

test("Maximum input size") {
    when {
        process {
            """
            input[0] = [
                [ id: 'test' ],
                file('large_file.txt', checkIfExists: true)  // Very large file
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            // Verify resource usage within limits
            { 
                assert process.trace.memory <= params.max_memory.replace('GB', '').toLong() * 1000000000
            },
            { assert snapshot(process.out).match("large_input") }
        )
    }
}
```

## Performance Testing

### Resource Monitoring

```groovy
test("Resource usage verification") {
    when {
        process {
            """
            input[0] = [
                [ id: 'test' ],
                file('benchmark_data.txt', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            // Check execution time
            { assert process.trace.duration.toMillis() < 300000 }, // 5 minutes
            // Check memory usage
            { assert process.trace.peakRss < 4000000000 }, // 4GB
            // Check CPU efficiency
            { assert process.trace.cpuTime / process.trace.duration.toMillis() > 0.8 },
            { assert snapshot(process.out).match() }
        )
    }
}
```

### Scalability Testing

```groovy
[1, 2, 4, 8].each { cpu_count ->
    test("Scaling with ${cpu_count} CPUs") {
        when {
            params {
                max_cpus = cpu_count
            }
            process {
                """
                input[0] = [
                    [ id: 'test', cpus: ${cpu_count} ],
                    file('scalability_test.txt', checkIfExists: true)
                ]
                """
            }
        }
        then {
            assertAll(
                { assert process.success },
                { assert process.trace.cpus <= cpu_count },
                // Performance should improve with more CPUs
                { 
                    if (cpu_count > 1) {
                        assert process.trace.duration.toMillis() < baseline_duration
                    }
                },
                { assert snapshot(process.out).match("cpu_${cpu_count}") }
            )
        }
    }
}
```

## Integration Testing

### Cross-Module Integration

```groovy
nextflow_workflow {
    name "Test Module Integration"
    script "../integration_test.nf"
    workflow "MODULE_INTEGRATION"
    
    test("Module chain integration") {
        when {
            workflow {
                """
                input[0] = Channel.fromList([
                    [[ id: 'sample1' ], file('sample1.fastq.gz', checkIfExists: true)],
                    [[ id: 'sample2' ], file('sample2.fastq.gz', checkIfExists: true)]
                ])
                """
            }
        }
        then {
            assertAll(
                { assert workflow.success },
                // Verify each module in chain executed
                { assert workflow.trace.findAll { it.name == 'FASTQC' }.size() == 2 },
                { assert workflow.trace.findAll { it.name == 'TRIMMOMATIC' }.size() == 2 },
                { assert workflow.trace.findAll { it.name == 'BWA_MEM' }.size() == 2 },
                // Verify data flow between modules
                {
                    workflow.out.bam.each { meta, bam ->
                        assert meta.id in ['sample1', 'sample2']
                        assert bam.exists()
                    }
                },
                { assert snapshot(workflow.out).match() }
            )
        }
    }
}
```

### End-to-End Testing

```groovy
nextflow_pipeline {
    name "Test Full Pipeline"
    script "main.nf"
    
    test("Complete pipeline execution") {
        when {
            params {
                input = 'tests/full_samplesheet.csv'
                genome = 'GRCh38'
                outdir = './e2e_results'
            }
        }
        then {
            assertAll(
                { assert workflow.success },
                // Verify all expected outputs
                { assert new File("${params.outdir}/multiqc/multiqc_report.html").exists() },
                { assert new File("${params.outdir}/pipeline_info/execution_report.html").exists() },
                // Check sample-specific outputs
                {
                    ['sample1', 'sample2', 'sample3'].each { sample ->
                        assert new File("${params.outdir}/bwa/${sample}.sorted.bam").exists()
                        assert new File("${params.outdir}/variants/${sample}.vcf.gz").exists()
                    }
                },
                // Verify pipeline metrics
                { assert workflow.trace.succeeded().size() >= 15 },
                { assert workflow.duration.toHours() < 2 },
                { assert snapshot(workflow.out).match() }
            )
        }
    }
}
```

## Custom Assertion Patterns

### Custom Validation Functions

```groovy
// Define custom validation functions
def validateFastqOutput(fastqFile) {
    def lines = fastqFile.readLines()
    assert lines.size() % 4 == 0  // FASTQ format validation
    assert lines.every { it.trim().length() > 0 }  // No empty lines
    return true
}

def validateBAMHeader(bamFile) {
    def header = "samtools view -H ${bamFile}".execute().text
    assert header.contains("@HD")  // Header present
    assert header.contains("@SQ")  // Sequence dictionary present
    return true
}

test("Custom validation") {
    when {
        process {
            """
            input[0] = [
                [ id: 'test' ],
                file('input.fastq.gz', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            // Use custom validation
            {
                process.out.fastq.each { meta, fastq ->
                    validateFastqOutput(fastq)
                }
            },
            {
                process.out.bam.each { meta, bam ->
                    validateBAMHeader(bam)
                }
            },
            { assert snapshot(process.out).match() }
        )
    }
}
```

### Conditional Assertions

```groovy
then {
    assertAll(
        { assert process.success },
        // Conditional assertions based on input
        {
            if (meta.single_end) {
                assert process.out.bam.size() == 1
            } else {
                assert process.out.bam.size() == 1
                assert process.out.metrics.size() == 1
            }
        },
        // Platform-specific assertions
        {
            def os = System.getProperty("os.name").toLowerCase()
            if (os.contains("linux")) {
                // Linux-specific checks
                assert process.trace.container != null
            }
        },
        { assert snapshot(process.out).match() }
    )
}
```

## Next Steps

Continue to [Test Data Management](./11_test_data_management.md) to learn about organizing and managing test datasets. 