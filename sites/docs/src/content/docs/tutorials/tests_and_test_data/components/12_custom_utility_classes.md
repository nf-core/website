---
title: "13. Custom Utility Classes"
subtitle: Creating reusable test utilities
weight: 130
---

## Test Utility Organization

### Utility Class Structure

```groovy
// tests/utils/TestUtils.groovy
class TestUtils {
    
    // File validation utilities
    static boolean validateFastq(file) {
        def lines = file.readLines()
        return lines.size() % 4 == 0 && 
               lines[0].startsWith("@") && 
               lines[2].startsWith("+")
    }
    
    static boolean validateFasta(file) {
        def content = file.text
        return content.contains(">") && 
               content.matches(/(?s)>[^\n]+\n[ATCGN\n]+/)
    }
    
    // Data generation utilities
    static File generateSyntheticFastq(String filename, int numReads = 100) {
        def file = new File(filename)
        def bases = ['A', 'T', 'C', 'G']
        def random = new Random()
        
        file.withWriter { writer ->
            (1..numReads).each { i ->
                writer.println("@read_${i}")
                writer.println((1..100).collect { bases[random.nextInt(4)] }.join(''))
                writer.println("+")
                writer.println("I" * 100)
            }
        }
        return file
    }
    
    // Snapshot utilities
    static Map sanitizeOutput(Map output) {
        // Remove timestamps and variable content
        output.collectEntries { key, value ->
            if (value instanceof String) {
                value = value.replaceAll(/\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}/, 'TIMESTAMP')
                value = value.replaceAll(/\/tmp\/\w+/, '/tmp/TEMP_DIR')
            }
            [key, value]
        }
    }
}
```

### Using Utility Classes

```groovy
// Include utility class
@Grab('tests/utils/TestUtils.groovy')
import TestUtils

nextflow_process {
    name "Test Process with Utilities"
    script "../main.nf"
    process "TOOL"
    
    test("Using test utilities") {
        setup {
            // Generate test data using utility
            TestUtils.generateSyntheticFastq("tests/synthetic.fastq", 50)
        }
        
        when {
            process {
                """
                input[0] = [
                    [ id:'test' ],
                    file('tests/synthetic.fastq', checkIfExists: true)
                ]
                """
            }
        }
        
        then {
            assertAll(
                { assert process.success },
                // Validate output using utility
                {
                    process.out.fastq.each { meta, fastq ->
                        assert TestUtils.validateFastq(fastq)
                    }
                },
                // Sanitize before snapshot
                { assert snapshot(TestUtils.sanitizeOutput(process.out)).match() }
            )
        }
    }
}
```

## Data Validation Utilities

### Comprehensive File Validators

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
    
    static class VcfValidator {
        static boolean isValid(file) {
            def lines = file.readLines()
            return lines.any { it.startsWith("##fileformat=VCF") } &&
                   lines.any { it.startsWith("#CHROM") }
        }
        
        static Map getStats(file) {
            def lines = file.readLines()
            def variants = lines.findAll { !it.startsWith("#") }
            
            return [
                num_variants: variants.size(),
                chromosomes: variants.collect { it.split('\t')[0] }.unique().size(),
                variant_types: variants.collect { 
                    def cols = it.split('\t')
                    cols[3].length() == 1 && cols[4].length() == 1 ? 'SNP' : 'INDEL'
                }.countBy { it }
            ]
        }
    }
}
```

## Test Data Generators

### Genomic Data Generator

```groovy
// tests/utils/DataGenerators.groovy
class DataGenerators {
    
    static class GenomicDataGenerator {
        private static Random random = new Random()
        private static String[] bases = ['A', 'T', 'C', 'G']
        private static String[] qualities = ['!', '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '-', '.', '/', 
                                           '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', 
                                           '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
        
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
        
        static File generateVcf(String filename, Map options = [:]) {
            def numVariants = options.numVariants ?: 100
            def chromosomes = options.chromosomes ?: ['chr1', 'chr2', 'chr3']
            def samples = options.samples ?: ['sample1']
            
            def file = new File(filename)
            file.withWriter { writer ->
                // Header
                writer.println("##fileformat=VCFv4.2")
                writer.println("##source=TestDataGenerator")
                writer.println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
                
                // Column headers
                def headerLine = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
                samples.each { sample ->
                    headerLine += "\t${sample}"
                }
                writer.println(headerLine)
                
                // Variants
                (1..numVariants).each { i ->
                    def chr = chromosomes[random.nextInt(chromosomes.size())]
                    def pos = random.nextInt(1000000) + 1
                    def ref = bases[random.nextInt(4)]
                    def alt = bases[random.nextInt(4)]
                    def qual = random.nextInt(100) + 20
                    
                    def line = "${chr}\t${pos}\tvar${i}\t${ref}\t${alt}\t${qual}\tPASS\t.\tGT"
                    samples.each { sample ->
                        def gt = random.nextBoolean() ? "0/1" : "1/1"
                        line += "\t${gt}"
                    }
                    writer.println(line)
                }
            }
            return file
        }
    }
}
```

## Assertion Utilities

### Custom Assertion Library

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
    
    static void assertFileLineCount(file, expectedCount, String message = null) {
        def fileObj = file instanceof File ? file : new File(file.toString())
        def actualCount = fileObj.readLines().size()
        def msg = message ?: "File should have ${expectedCount} lines, but has ${actualCount}"
        assert actualCount == expectedCount, msg
    }
    
    static void assertFilesHaveSameChecksum(file1, file2, String message = null) {
        def checksum1 = file1.digest('MD5')
        def checksum2 = file2.digest('MD5')
        def msg = message ?: "Files should have same checksum"
        assert checksum1 == checksum2, msg
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
    
    static void assertOutputChannelStructure(channel, Map expectedStructure) {
        channel.each { item ->
            if (item instanceof List && item.size() >= 2) {
                def meta = item[0]
                def files = item[1..-1]
                
                expectedStructure.metaFields?.each { field ->
                    assert meta.containsKey(field), "Meta missing required field: ${field}"
                }
                
                if (expectedStructure.fileCount) {
                    assert files.size() == expectedStructure.fileCount,
                           "Expected ${expectedStructure.fileCount} files, got ${files.size()}"
                }
                
                expectedStructure.fileExtensions?.each { ext ->
                    assert files.any { it.toString().endsWith(ext) },
                           "Missing file with extension: ${ext}"
                }
            }
        }
    }
}
```

## Snapshot Utilities

### Advanced Snapshot Management

```groovy
// tests/utils/SnapshotUtils.groovy
class SnapshotUtils {
    
    static Map normalizeProcessOutput(processOut) {
        // Normalize file paths and timestamps
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
    
    static Map createSelectiveSnapshot(processOut, List channels) {
        def selective = [:]
        channels.each { channel ->
            if (processOut.containsKey(channel)) {
                selective[channel] = processOut[channel]
            }
        }
        return normalizeProcessOutput(selective)
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

## Configuration Utilities

### Dynamic Configuration Generator

```groovy
// tests/utils/ConfigUtils.groovy
class ConfigUtils {
    
    static String generateTestConfig(Map options) {
        def config = []
        
        // Process configuration
        if (options.processes) {
            config << "process {"
            options.processes.each { processName, processConfig ->
                config << "    withName: '${processName}' {"
                processConfig.each { key, value ->
                    config << "        ${key} = ${formatConfigValue(value)}"
                }
                config << "    }"
            }
            config << "}"
        }
        
        // Parameter configuration
        if (options.params) {
            config << "params {"
            options.params.each { key, value ->
                config << "    ${key} = ${formatConfigValue(value)}"
            }
            config << "}"
        }
        
        // Profile configuration
        if (options.profiles) {
            config << "profiles {"
            options.profiles.each { profileName, profileConfig ->
                config << "    ${profileName} {"
                profileConfig.each { key, value ->
                    config << "        ${key} = ${formatConfigValue(value)}"
                }
                config << "    }"
            }
            config << "}"
        }
        
        return config.join('\n')
    }
    
    private static String formatConfigValue(value) {
        if (value instanceof String) {
            return "'${value}'"
        } else if (value instanceof Map) {
            def entries = value.collect { k, v -> "${k}: ${formatConfigValue(v)}" }
            return "[${entries.join(', ')}]"
        } else if (value instanceof List) {
            def items = value.collect { formatConfigValue(it) }
            return "[${items.join(', ')}]"
        } else {
            return value.toString()
        }
    }
    
    static File createConfigFile(String filename, Map options) {
        def configContent = generateTestConfig(options)
        def file = new File(filename)
        file.text = configContent
        return file
    }
    
    static Map loadTestDataConfig(String configFile) {
        def config = new ConfigSlurper().parse(new File(configFile).toURI().toURL())
        return config.test_data ?: [:]
    }
    
    static String resolveTestDataPath(String path, Map testDataConfig) {
        if (path.startsWith('http')) {
            return path  // Already a URL
        }
        
        def basePath = testDataConfig.base_path ?: 'tests/data'
        return "${basePath}/${path}"
    }
}
```

## Integration with Test Framework

### Utility Integration Example

```groovy
// Example test using multiple utilities
@Grab('tests/utils/TestUtils.groovy')
@Grab('tests/utils/DataGenerators.groovy')
@Grab('tests/utils/AssertionUtils.groovy')
@Grab('tests/utils/SnapshotUtils.groovy')

import TestUtils
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
                
                // Output structure validation
                {
                    AssertionUtils.assertOutputChannelStructure(
                        process.out.bam,
                        [
                            metaFields: ['id', 'type'],
                            fileCount: 1,
                            fileExtensions: ['.bam']
                        ]
                    )
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

Continue to [Error Testing & Edge Cases](./13_error_testing_edge_cases.md) to learn about testing failure scenarios. 