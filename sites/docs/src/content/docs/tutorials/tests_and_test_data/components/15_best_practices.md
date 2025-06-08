---
title: "16. Best Practices"
subtitle: Testing best practices and conventions from nf-core/methylseq
weight: 160
---

## Test Organization Best Practices

### Directory Structure

Follow consistent directory structure as used in nf-core pipelines:

```
tests/
├── default.nf.test                    # Main pipeline test
├── bwameth.nf.test                    # Alternative pathway tests
├── bismark_targeted_sequencing.nf.test # Specialized feature tests
├── .nftignore                         # Files to ignore in snapshots
├── nextflow.config                    # Test-specific configuration
└── default.nf.test.snap              # Snapshot files
```

### Test Naming Conventions

Based on methylseq's current implementation:

```groovy
nextflow_pipeline {
    name "Test pipeline"
    script "../main.nf"
    config "./nextflow.config"
    tag "pipeline"
    tag "cpu"
    
    // Good: Descriptive parameter-based names
    test("-profile test") { /* Default pathway */ }
    test("Params: bwameth") { /* Alternative aligner */ }
    test("Params: bismark | run_targeted_sequencing") { /* Feature-specific */ }
    test("Params: bismark | skip_trimming") { /* Skip options */ }
    
    // Avoid: Generic test names
    test("test1") { /* ... */ }
    test("basic") { /* ... */ }
    test("default") { /* ... */ }
}
```

### Pipeline Test Pattern

Standard nf-core pipeline test structure from methylseq:

```groovy
test("-profile test") {
    when {
        params {
            outdir = "$outputDir"
        }
    }
    then {
        // Standard pattern: stable names, stable content, specific file types
        def stable_name = getAllFilesFromDir(params.outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}'])
        def stable_path = getAllFilesFromDir(params.outdir, ignoreFile: 'tests/.nftignore')
        def bam_files = getAllFilesFromDir(params.outdir, include: ['**/*.bam'])
        assertAll(
            { assert workflow.success },
            { assert snapshot(
                workflow.trace.succeeded().size(),
                removeNextflowVersion("$outputDir/pipeline_info/nf_core_methylseq_software_mqc_versions.yml"),
                stable_name,
                stable_path,
                bam_files.collect{ file -> [ file.getName(), bam(file.toString()).getReadsMD5() ] }
            ).match() }
        )
    }
}
```

## nf-core Testing Standards

### Essential Components

Every nf-core pipeline test **MUST** include:

1. **Task count verification**: `workflow.trace.succeeded().size()`
2. **Version file processing**: `removeNextflowVersion()` for reproducible tests
3. **File structure verification**: `getAllFilesFromDir()` with stable names
4. **Content verification**: `getAllFilesFromDir()` with `.nftignore` filtering
5. **Format-specific verification**: e.g., `bam().getReadsMD5()` for BAM files

### Required Plugins

All nf-core pipelines should use these plugins in `nf-test.config`:

```groovy
plugins {
    load "nft-bam@0.5.0"      // BAM file testing utilities
    load "nft-utils@0.0.3"    // File management and utility functions
}
```

### Standard Test Configuration

Based on methylseq's implementation:

```groovy
// nf-test.config
config {
    testsDir "."
    workDir System.getenv("NFT_WORKDIR") ?: ".nf-test"
    configFile "tests/nextflow.config"
    ignore 'modules/nf-core/**/*', 'subworkflows/nf-core/**/*'
    profile "test"
    triggers 'nextflow.config', 'nf-test.config', 'conf/test.config', 'tests/nextflow.config', 'tests/.nftignore'
    plugins {
        load "nft-bam@0.5.0"
        load "nft-utils@0.0.3"
    }
}
```

### Test Data Management

Use standardized test data paths:

```groovy
params {
    // For modules
    modules_testdata_base_path = 's3://ngi-igenomes/testdata/nf-core/modules/'
    // For pipelines
    pipelines_testdata_base_path = 'https://raw.githubusercontent.com/nf-core/test-datasets/methylseq/'
    
    // Pipeline-specific test data
    input = "${projectDir}/assets/samplesheet.csv"
    fasta = "${params.pipelines_testdata_base_path}/reference/genome.fa.gz"
}
```

## Test Data Best Practices

### Use Appropriate Data Sizes

```groovy
// CI tests: minimal data
test("CI test") {
    when {
        process {
            """
            input[0] = [
                [ id:'test' ],
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

// Full tests: realistic data
test("Full test", tags: ["full"]) {
    when {
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('tests/data/large_sample.fastq.gz', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            { assert snapshot(process.out).match("full_dataset") }
        )
    }
}
```

### Data Documentation

```groovy
// Document test data in test files
test("Reference genome alignment") {
    // Test data: Human chromosome 22 (50Mb subset)
    // Source: 1000 Genomes Project
    // Format: BWA-indexed FASTA
    // Size: ~15MB compressed
    when {
        process {
            """
            input[0] = [
                [ id:'test', genome:'GRCh38' ],
                file('tests/data/chr22.fasta.gz', checkIfExists: true)
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

## Assertion Best Practices

### Comprehensive Assertions

```groovy
test("Comprehensive output validation") {
    when {
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('input.fastq.gz', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            // Basic success check
            { assert process.success },
            
            // Output existence
            { assert process.out.html.size() == 1 },
            { assert process.out.zip.size() == 1 },
            { assert process.out.versions.size() == 1 },
            
            // File content validation
            {
                process.out.html.each { meta, html ->
                    assert html.exists()
                    assert html.size() > 0
                    assert html.text.contains("FastQC Report")
                }
            },
            
            // Metadata validation
            {
                process.out.html.each { meta, html ->
                    assert meta.id == 'test'
                    assert meta.containsKey('single_end')
                }
            },
            
            // Snapshot comparison
            { assert snapshot(process.out).match() }
        )
    }
}
```

### Error-Specific Assertions

```groovy
test("Invalid input handling") {
    when {
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('invalid.txt', checkIfExists: true)  // Not FASTQ
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.failed },
            { assert process.exitStatus == 1 },
            // Check for specific error messages
            { assert process.stderr.contains("not in FASTQ format") || 
                     process.stderr.contains("Invalid file format") },
            // Ensure no partial outputs
            { assert process.out.html.size() == 0 },
            { assert process.out.zip.size() == 0 }
        )
    }
}
```

## Snapshot Management Best Practices

### Selective Snapshots

```groovy
test("Selective snapshot example") {
    when {
        process {
            """
            input[0] = [
                [ id:'test' ],
                file('input.fastq.gz', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            
            // Snapshot only stable outputs
            { assert snapshot(
                process.out.html,
                process.out.zip
            ).match("main_outputs") },
            
            // Separate snapshot for versions (may change frequently)
            { assert snapshot(process.out.versions).match("versions") }
        )
    }
}
```

### Content Filtering

```groovy
test("Filtered snapshot content") {
    when {
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
            
            // Filter variable content before snapshotting
            { assert snapshot(
                process.out.stats.collect { meta, stats ->
                    // Remove timestamps and process IDs
                    def content = stats.text
                        .replaceAll(/\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}/, 'TIMESTAMP')
                        .replaceAll(/PID:\d+/, 'PID:FILTERED')
                    [meta: [id: meta.id], content: content]
                }
            ).match("filtered_stats") }
        )
    }
}
```

## Configuration Best Practices

### Environment-Specific Configs

```groovy
// tests/nextflow.config
profiles {
    test {
        params {
            max_memory = '6.GB'
            max_cpus = 2
            max_time = '6.h'
        }
        process {
            withName: 'FASTQC' {
                ext.args = '--quiet'
            }
        }
    }
    
    test_full {
        params {
            max_memory = '64.GB'
            max_cpus = 16
            max_time = '48.h'
        }
    }
    
    ci {
        params {
            max_memory = '6.GB'
            max_cpus = 2
            max_time = '1.h'
        }
        docker.enabled = true
        singularity.enabled = false
    }
}
```

### Parameterized Testing

```groovy
// Use parameters for flexible testing
def testCases = [
    [quality: 20, description: "low quality threshold"],
    [quality: 30, description: "standard quality threshold"],
    [quality: 40, description: "high quality threshold"]
]

testCases.each { testCase ->
    test("Quality filtering - ${testCase.description}") {
        when {
            params {
                quality_threshold = testCase.quality
            }
            process {
                """
                input[0] = [
                    [ id:'test', quality: ${testCase.quality} ],
                    file('input.fastq.gz', checkIfExists: true)
                ]
                """
            }
        }
        then {
            assertAll(
                { assert process.success },
                { assert snapshot(process.out).match("quality_${testCase.quality}") }
            )
        }
    }
}
```

## Performance Best Practices

### Resource-Aware Testing

```groovy
test("Resource usage validation") {
    when {
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
            
            // Validate resource usage
            { assert process.trace.duration.toMinutes() < 10 },
            { assert process.trace.peakRss < 2_000_000_000 }, // 2GB
            { assert process.trace.cpus <= 4 },
            
            // Performance regression detection
            {
                def expectedDuration = 5 * 60 * 1000 // 5 minutes in ms
                def actualDuration = process.trace.duration.toMillis()
                def tolerance = 0.2 // 20% tolerance
                
                assert actualDuration < expectedDuration * (1 + tolerance)
            },
            
            { assert snapshot(process.out).match() }
        )
    }
}
```

### Scalability Testing

```groovy
[1, 2, 4].each { cpus ->
    test("CPU scaling - ${cpus} cores") {
        when {
            params {
                max_cpus = cpus
            }
            process {
                """
                input[0] = [
                    [ id:'test', cpus: ${cpus} ],
                    file('large_input.bam', checkIfExists: true)
                ]
                """
            }
        }
        then {
            assertAll(
                { assert process.success },
                { assert process.trace.cpus <= cpus },
                { assert snapshot(process.out).match("cpus_${cpus}") }
            )
        }
    }
}
```

## Error Handling Best Practices

### Graceful Failure Testing

```groovy
test("Graceful failure handling") {
    when {
        process {
            """
            # Test graceful handling of missing dependencies
            if ! command -v optional_tool >/dev/null 2>&1; then
                echo "Optional tool not found, using fallback method"
                # Implement fallback logic
                basic_tool input.txt > output.txt
            else
                optional_tool input.txt > output.txt
            fi
            
            input[0] = [
                [ id:'test' ],
                file('input.txt', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            // Should succeed regardless of optional tool availability
            { assert snapshot(process.out).match() }
        )
    }
}
```

### Input Validation

```groovy
test("Input validation") {
    when {
        process {
            """
            # Comprehensive input validation
            validate_input() {
                local file="\$1"
                
                # Check file exists and is readable
                if [[ ! -r "\$file" ]]; then
                    echo "ERROR: Cannot read file: \$file" >&2
                    return 1
                fi
                
                # Check file size
                if [[ ! -s "\$file" ]]; then
                    echo "ERROR: File is empty: \$file" >&2
                    return 1
                fi
                
                # Format-specific validation
                if [[ "\$file" == *.fastq.gz ]]; then
                    if ! zcat "\$file" | head -4 | awk 'NR==1{if(!\$0~/^@/)exit 1}NR==3{if(!\$0~/^\+/)exit 1}'; then
                        echo "ERROR: Invalid FASTQ format: \$file" >&2
                        return 1
                    fi
                fi
                
                return 0
            }
            
            # Validate all inputs
            for input_file in input.fastq.gz reference.fa; do
                if ! validate_input "\$input_file"; then
                    exit 1
                fi
            done
            
            echo "Input validation passed"
            
            input[0] = [
                [ id:'test' ],
                file('input.fastq.gz', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            { assert process.stdout.contains("Input validation passed") },
            { assert snapshot(process.out).match() }
        )
    }
}
```

## Documentation Best Practices

### Self-Documenting Tests

```groovy
nextflow_process {
    name "Test Process BWA_MEM"
    script "../main.nf"
    process "BWA_MEM"
    
    /**
     * Test BWA-MEM alignment with single-end reads
     * 
     * This test verifies:
     * - Basic alignment functionality
     * - SAM output format
     * - Proper handling of single-end data
     * - Resource usage within limits
     */
    test("Single-end alignment") {
        when {
            process {
                """
                // Test data: E. coli reads (100k reads, ~15MB)
                // Reference: E. coli K-12 MG1655 complete genome
                input[0] = [
                    [ id:'test', single_end:true ],
                    file(params.modules_testdata_base_path + 'genomics/prokaryotes/bacteroides_fragilis/illumina/fastq/test1_1.fastq.gz', checkIfExists: true)
                ]
                input[1] = file(params.modules_testdata_base_path + 'genomics/prokaryotes/bacteroides_fragilis/genome/genome.fasta', checkIfExists: true)
                input[2] = file(params.modules_testdata_base_path + 'genomics/prokaryotes/bacteroides_fragilis/genome/genome.fasta.fai', checkIfExists: true)
                """
            }
        }
        then {
            assertAll(
                { assert process.success },
                
                // Verify SAM output structure
                {
                    process.out.sam.each { meta, sam ->
                        def lines = sam.readLines()
                        assert lines.any { it.startsWith("@HD") } // Header present
                        assert lines.any { it.startsWith("@SQ") } // Reference sequences present
                        assert lines.count { !it.startsWith("@") } > 0 // Alignment records present
                    }
                },
                
                // Resource validation
                { assert process.trace.duration.toMinutes() < 5 },
                { assert process.trace.peakRss < 1_000_000_000 }, // 1GB
                
                { assert snapshot(process.out).match() }
            )
        }
    }
}
```

### Test Metadata

```groovy
// tests/tags.yml
BWA_MEM:
  - modules
  - modules_nfcore
  - alignment
  - genomics
  - bwa
  
test_descriptions:
  "Single-end alignment": "Tests basic BWA-MEM functionality with single-end reads"
  "Paired-end alignment": "Tests BWA-MEM with paired-end reads and proper mate handling"
  "Large genome alignment": "Performance test with human genome reference"
  
test_requirements:
  "Single-end alignment":
    memory: "1GB"
    time: "5min"
    data_size: "15MB"
  "Large genome alignment":
    memory: "8GB"
    time: "30min"
    data_size: "500MB"
```

## Maintenance Best Practices

### Version Compatibility

```groovy
test("Version compatibility check") {
    when {
        process {
            """
            # Check tool version and adjust behavior accordingly
            TOOL_VERSION=\$(my_tool --version | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+')
            MAJOR_VERSION=\$(echo \$TOOL_VERSION | cut -d. -f1)
            
            if [[ \$MAJOR_VERSION -ge 2 ]]; then
                # New version syntax
                my_tool --new-flag input.txt > output.txt
            else
                # Legacy version syntax
                my_tool --old-flag input.txt > output.txt
            fi
            
            input[0] = [
                [ id:'test', tool_version: "\$TOOL_VERSION" ],
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

### Deprecation Handling

```groovy
test("Feature deprecation handling") {
    when {
        process {
            """
            # Handle deprecated parameters gracefully
            if [[ "\${params.legacy_parameter:-}" != "" ]]; then
                echo "WARNING: legacy_parameter is deprecated, use new_parameter instead" >&2
                NEW_PARAM="\${params.legacy_parameter}"
            else
                NEW_PARAM="\${params.new_parameter}"
            fi
            
            my_tool --param "\$NEW_PARAM" input.txt > output.txt
            
            input[0] = [
                [ id:'test' ],
                file('input.txt', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            { 
                if (params.legacy_parameter) {
                    assert process.stderr.contains("WARNING: legacy_parameter is deprecated")
                }
            },
            { assert snapshot(process.out).match() }
        )
    }
}
```

## Continuous Improvement

### Test Coverage Analysis

```bash
#!/bin/bash
# scripts/analyze_test_coverage.sh

echo "=== nf-core Test Coverage Analysis ==="

# Find all modules
MODULES=$(find modules/nf-core -name "main.nf" | wc -l)
echo "Total modules: $MODULES"

# Find modules with tests
TESTED_MODULES=$(find modules/nf-core -name "main.nf.test" | wc -l)
echo "Modules with tests: $TESTED_MODULES"

# Calculate coverage
COVERAGE=$(echo "scale=2; $TESTED_MODULES * 100 / $MODULES" | bc)
echo "Test coverage: ${COVERAGE}%"

# Find modules without tests
echo -e "\nModules missing tests:"
find modules/nf-core -name "main.nf" -exec dirname {} \; | while read module_dir; do
    if [[ ! -f "$module_dir/tests/main.nf.test" ]]; then
        echo "  - $module_dir"
    fi
done
```

### Performance Monitoring

```groovy
// Add to all performance-critical tests
test("Performance regression check") {
    when {
        process {
            """
            START_TIME=\$(date +%s.%N)
            
            # Run the actual process
            my_tool input.txt > output.txt
            
            END_TIME=\$(date +%s.%N)
            DURATION=\$(echo "\$END_TIME - \$START_TIME" | bc)
            
            echo "Execution time: \${DURATION}s"
            
            input[0] = [
                [ id:'test', execution_time: "\$DURATION" ],
                file('input.txt', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            
            // Performance regression check (baseline: 30s, tolerance: 20%)
            {
                def executionTime = process.stdout.find(/Execution time: ([0-9.]+)s/) { match, time ->
                    time as Double
                }
                def baseline = 30.0
                def tolerance = 0.2
                
                assert executionTime < baseline * (1 + tolerance), 
                       "Performance regression detected: ${executionTime}s > ${baseline * (1 + tolerance)}s"
            },
            
            { assert snapshot(process.out).match() }
        )
    }
}
```

## Next Steps

Continue to [Troubleshooting Guide](./16_troubleshooting_guide.md) to learn how to debug common testing issues. 