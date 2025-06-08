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

## Using Custom Functions

### Custom Utility Classes

#### File Validation Utilities

```groovy
// tests/utils/FileValidators.groovy
class FileValidators {
    
    static class FastqValidator {
        static boolean isValid(file) {
            def lines = file.readLines()
            if (lines.size() % 4 != 0) return false
            
            for (int i = 0; i < lines.size(); i += 4) {
                if (!lines[i].startsWith("@")) return false      // Header
                if (lines[i+1].isEmpty()) return false           // Sequence
                if (!lines[i+2].startsWith("+")) return false    // Separator
                if (lines[i+3].isEmpty()) return false           // Quality
                if (lines[i+1].length() != lines[i+3].length()) return false
            }
            return true
        }
        
        static Map getStats(file) {
            def lines = file.readLines()
            def numReads = lines.size() / 4
            def avgLength = lines.findAll { !it.startsWith("@") && !it.startsWith("+") }
                                .collect { it.length() }
                                .sum() / numReads
            
            return [
                num_reads: numReads,
                avg_length: avgLength,
                total_bases: numReads * avgLength
            ]
        }
    }
    
    static class BamValidator {
        static boolean isValid(file) {
            try {
                def header = "samtools view -H ${file}".execute().text
                return header.contains("@HD") && header.contains("VN:")
            } catch (Exception e) {
                return false
            }
        }
        
        static Map getStats(file) {
            def stats = [:]
            try {
                def idxstats = "samtools idxstats ${file}".execute().text
                def lines = idxstats.split('\n').findAll { it.trim() }
                
                stats.chromosomes = lines.size() - 1  // Exclude unmapped
                stats.mapped_reads = lines.collect { 
                    it.split('\t')[2] as Integer 
                }.sum()
                stats.unmapped_reads = lines[-1].split('\t')[3] as Integer
                
            } catch (Exception e) {
                stats.error = e.message
            }
            return stats
        }
    }
}
```

#### Data Generation Utilities

```groovy
// tests/utils/DataGenerators.groovy
class DataGenerators {
    
    static class GenomicDataGenerator {
        private static Random random = new Random()
        private static String[] bases = ['A', 'T', 'C', 'G']
        
        static File generateFasta(String filename, Map options = [:]) {
            def numSequences = options.numSequences ?: 5
            def seqLength = options.seqLength ?: 1000
            def prefix = options.prefix ?: "seq"
            
            def file = new File(filename)
            file.withWriter { writer ->
                (1..numSequences).each { i ->
                    writer.println(">${prefix}_${i}")
                    def sequence = (1..seqLength).collect { 
                        bases[random.nextInt(4)] 
                    }.join('')
                    
                    // Write in 80-character lines
                    sequence.toList().collate(80).each { chunk ->
                        writer.println(chunk.join(''))
                    }
                }
            }
            return file
        }
        
        static File generateFastq(String filename, Map options = [:]) {
            def numReads = options.numReads ?: 1000
            def readLength = options.readLength ?: 150
            def qualityOffset = options.qualityOffset ?: 33
            def prefix = options.prefix ?: "read"
            
            def file = new File(filename)
            file.withWriter { writer ->
                (1..numReads).each { i ->
                    // Header
                    writer.println("@${prefix}_${i}")
                    
                    // Sequence
                    def sequence = (1..readLength).collect { 
                        bases[random.nextInt(4)] 
                    }.join('')
                    writer.println(sequence)
                    
                    // Separator
                    writer.println("+")
                    
                    // Quality scores
                    def quality = (1..readLength).collect {
                        (char)(random.nextInt(40) + qualityOffset)
                    }.join('')
                    writer.println(quality)
                }
            }
            return file
        }
    }
}
```

#### Assertion Utilities

```groovy
// tests/utils/AssertionUtils.groovy
class AssertionUtils {
    
    static void assertFileExists(file, String message = null) {
        def msg = message ?: "File should exist: ${file}"
        assert file instanceof File ? file.exists() : new File(file.toString()).exists(), msg
    }
    
    static void assertFileNotEmpty(file, String message = null) {
        def fileObj = file instanceof File ? file : new File(file.toString())
        def msg = message ?: "File should not be empty: ${fileObj}"
        assert fileObj.exists() && fileObj.size() > 0, msg
    }
    
    static void assertFileContentMatches(file, pattern, String message = null) {
        def fileObj = file instanceof File ? file : new File(file.toString())
        def content = fileObj.text
        def msg = message ?: "File content should match pattern: ${pattern}"
        assert content.matches(pattern), msg
    }
    
    static void assertProcessResourceUsage(trace, Map limits, String message = null) {
        if (limits.maxMemory) {
            def memoryBytes = limits.maxMemory.replace('GB', '').toLong() * 1_000_000_000
            assert trace.peakRss <= memoryBytes, 
                   message ?: "Memory usage ${trace.peakRss} exceeds limit ${memoryBytes}"
        }
        
        if (limits.maxCpus) {
            assert trace.cpus <= limits.maxCpus,
                   message ?: "CPU usage ${trace.cpus} exceeds limit ${limits.maxCpus}"
        }
        
        if (limits.maxTime) {
            def timeMillis = limits.maxTime.replace('h', '').toLong() * 3_600_000
            assert trace.duration.toMillis() <= timeMillis,
                   message ?: "Execution time ${trace.duration} exceeds limit ${limits.maxTime}"
        }
    }
}
```

### Container and Environment Testing

#### Multi-Container Testing

```groovy
test("Docker container") {
    when {
        process {
            """
            container 'biocontainers/samtools:1.15.1--h1170115_0'
            
            input[0] = [
                [ id:'test' ],
                file('test.bam', checkIfExists: true)
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

test("Container environment validation") {
    when {
        process {
            """
            // Validate container environment
            echo "Container info:"
            uname -a
            echo "Available tools:"
            which samtools && samtools --version
            which bcftools && bcftools --version
            
            input[0] = [
                [ id:'test' ],
                file('test.bam', checkIfExists: true)
            ]
            """
        }
    }
    then {
        assertAll(
            { assert process.success },
            { assert process.stdout.contains("samtools") },
            { assert snapshot(process.out).match() }
        )
    }
}
```

### Mock and Stub Testing

#### Mocking External Services

```groovy
test("Mock external service") {
    setup {
        // Create mock service response
        def mockResponse = new File("tests/mock_api_response.json")
        mockResponse.text = '''
        {
            "status": "success",
            "data": {
                "gene_id": "ENSG00000123456",
                "gene_name": "TEST_GENE",
                "chromosome": "1"
            }
        }
        '''
    }
    
    when {
        process {
            """
            # Mock external API call
            if [[ "\${params.use_mock_api}" == "true" ]]; then
                echo "Using mock API response"
                cp tests/mock_api_response.json api_response.json
            else
                echo "Calling real API"
                curl -s "https://api.example.com/gene/ENSG00000123456" > api_response.json
            fi
            
            input[0] = [
                [ id:'test' ],
                file('gene_list.txt', checkIfExists: true)
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

#### Tool Stubbing

```groovy
test("Stubbed external tool") {
    setup {
        // Create stub for external tool
        def stubScript = new File("tests/stubs/external_tool")
        stubScript.text = '''#!/bin/bash
echo "STUB: external_tool called with arguments: $@"
echo "Creating mock output..."
echo "mock_result" > $2  # Assuming $2 is output file
exit 0
'''
        stubScript.setExecutable(true)
    }
    
    when {
        process {
            """
            # Add stub to PATH
            export PATH="tests/stubs:\$PATH"
            
            # Verify stub is used
            which external_tool
            
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
            { assert process.stdout.contains("STUB:") },
            { assert snapshot(process.out).match() }
        )
    }
}
```

### Snapshot and Output Normalization

#### Advanced Snapshot Management

```groovy
// tests/utils/SnapshotUtils.groovy
class SnapshotUtils {
    
    static Map normalizeProcessOutput(processOut) {
        def normalized = [:]
        
        processOut.each { channelName, channelData ->
            normalized[channelName] = channelData.collect { item ->
                if (item instanceof List && item.size() >= 2) {
                    def meta = item[0]
                    def files = item[1..-1]
                    
                    // Normalize meta
                    def normalizedMeta = normalizeMeta(meta)
                    
                    // Normalize files
                    def normalizedFiles = files.collect { file ->
                        normalizeFile(file)
                    }
                    
                    return [normalizedMeta] + normalizedFiles
                } else {
                    return normalizeFile(item)
                }
            }
        }
        
        return normalized
    }
    
    private static Map normalizeMeta(meta) {
        def normalized = [:]
        meta.each { key, value ->
            if (key in ['id', 'single_end']) {
                normalized[key] = value  // Keep essential fields
            }
            // Skip timestamp and path fields
        }
        return normalized
    }
    
    private static Object normalizeFile(file) {
        if (file instanceof File || file.toString().contains('/')) {
            // For files, return just filename and checksum
            def filename = file.toString().split('/').last()
            try {
                def checksum = file.digest('MD5')
                return "${filename}:md5:${checksum}"
            } catch (Exception e) {
                return "${filename}:exists:${file.exists()}"
            }
        }
        return file
    }
    
    static Map excludeVariableContent(processOut, List excludePatterns = []) {
        def defaultPatterns = [
            /\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}/,  // Timestamps
            /\/tmp\/\w+/,                            // Temp directories
            /process_id:\d+/,                        // Process IDs
            /execution_id:\w+-\w+-\w+-\w+-\w+/       // Execution IDs
        ]
        
        def allPatterns = defaultPatterns + excludePatterns
        
        def filtered = [:]
        processOut.each { channelName, channelData ->
            filtered[channelName] = channelData.collect { item ->
                def itemStr = item.toString()
                allPatterns.each { pattern ->
                    itemStr = itemStr.replaceAll(pattern, 'FILTERED')
                }
                return itemStr
            }
        }
        
        return filtered
    }
}
```

### Error Handling and Graceful Degradation

#### Graceful Tool Degradation

```groovy
test("Graceful tool degradation") {
    when {
        process {
            """
            # Try to use preferred tool, fallback to alternative
            if command -v fast_tool >/dev/null 2>&1; then
                echo "Using fast_tool"
                fast_tool --input input.txt --output output.txt
            elif command -v slow_tool >/dev/null 2>&1; then
                echo "fast_tool not available, using slow_tool"
                slow_tool --input input.txt --output output.txt
            else
                echo "No suitable tool found, using built-in method"
                # Built-in fallback implementation
                cat input.txt | sed 's/pattern/replacement/g' > output.txt
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
            { assert snapshot(process.out).match() }
        )
    }
}
```

#### Timeout and Retry Logic

```groovy
test("External tool with timeout") {
    when {
        process {
            """
            # Function to run command with timeout
            run_with_timeout() {
                local cmd="\$1"
                local timeout="\$2"
                
                timeout "\$timeout" bash -c "\$cmd" || {
                    echo "Command timed out after \$timeout seconds"
                    return 1
                }
            }
            
            # Function to retry command
            retry_command() {
                local cmd="\$1"
                local max_attempts=3
                local attempt=1
                
                while [[ \$attempt -le \$max_attempts ]]; do
                    echo "Attempt \$attempt of \$max_attempts"
                    if run_with_timeout "\$cmd" "300s"; then
                        echo "Command succeeded on attempt \$attempt"
                        return 0
                    else
                        echo "Command failed on attempt \$attempt"
                        attempt=\$((attempt + 1))
                        sleep 10
                    fi
                done
                
                echo "Command failed after \$max_attempts attempts"
                return 1
            }
            
            # Use retry logic for unreliable external command
            retry_command "external_unreliable_tool --input input.txt --output output.txt"
            
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
            { assert snapshot(process.out).match() }
        )
    }
}
```

### Utility Integration Example

```groovy
// Example test using multiple utilities
@Grab('tests/utils/DataGenerators.groovy')
@Grab('tests/utils/AssertionUtils.groovy')
@Grab('tests/utils/SnapshotUtils.groovy')

import DataGenerators
import AssertionUtils
import SnapshotUtils

nextflow_process {
    name "Test Process with Full Utilities"
    script "../main.nf"
    process "COMPREHENSIVE_TOOL"
    
    test("Comprehensive test with utilities") {
        setup {
            // Generate test data
            DataGenerators.GenomicDataGenerator.generateFastq(
                "tests/generated_reads.fastq", 
                [numReads: 100, readLength: 150]
            )
            
            DataGenerators.GenomicDataGenerator.generateFasta(
                "tests/reference.fa",
                [numSequences: 1, seqLength: 5000]
            )
        }
        
        when {
            process {
                """
                input[0] = [
                    [ id:'test', type:'comprehensive' ],
                    file('tests/generated_reads.fastq', checkIfExists: true)
                ]
                input[1] = file('tests/reference.fa', checkIfExists: true)
                """
            }
        }
        
        then {
            assertAll(
                { assert process.success },
                
                // File existence checks
                {
                    AssertionUtils.assertFileExists(process.out.bam[0][1])
                    AssertionUtils.assertFileNotEmpty(process.out.bam[0][1])
                },
                
                // Resource usage validation
                {
                    AssertionUtils.assertProcessResourceUsage(
                        process.trace,
                        [maxMemory: '4GB', maxCpus: 2, maxTime: '1h']
                    )
                },
                
                // Normalized snapshot
                { 
                    assert snapshot(
                        SnapshotUtils.normalizeProcessOutput(process.out)
                    ).match() 
                }
            )
        }
    }
}
```

## Next Steps

Continue to [Test Data Management](./11_test_data_management.md) to learn about organizing and managing test datasets. 