---
# Copied from docs/utilities/NfcoreVersionUtils.md of nf-core/nf-core-utils 0.5.0 by .github/workflows/add-nf-core-utils-docs.yml: edit it there
title: "NfcoreVersionUtils"
shortTitle: "NfcoreVersionUtils"
weight: 9
md_github_url: https://github.com/nf-core/nf-core-utils/blob/main/docs/utilities/NfcoreVersionUtils.md
---

# NfcoreVersionUtils

Version management and topic channel utilities for nf-core pipelines.

## Overview

The `NfcoreVersionUtils` utility provides comprehensive version management functionality for nf-core pipelines, supporting both traditional file-based approaches and modern topic channel patterns. This utility enables progressive migration from legacy `versions.yml` files to the more flexible and efficient topic channel system.

### Topic Channel Migration Strategy

The utility supports a three-phase migration approach:

1. **Legacy Stage**: Traditional `versions.yml` file processing
2. **Transition Stage**: Mixed approach supporting both files and topic channels
3. **Modern Stage**: Pure topic channel approach with enhanced functionality

This migration strategy ensures backward compatibility while enabling pipelines to benefit from modern Nextflow features.

## Available Functions

### `getWorkflowVersion()`

**Description:**  
Generates a formatted version string for the workflow, combining version information with git commit details when available. This is the primary function for getting workflow version information.

**Function Signature:**

```nextflow
String getWorkflowVersion(String version = null, String commitId = null)
```

**Parameters:**

- `version` (String, optional): Explicit version string to use (defaults to workflow.manifest.version)
- `commitId` (String, optional): Explicit commit ID to use (defaults to workflow.commitId)

**Returns:**

- `String`: Formatted version string (e.g., `v1.2.3-gabcdef1`)

**Usage Example:**

```nextflow
include { getWorkflowVersion } from 'plugin/nf-core-utils'

// Basic usage - uses workflow metadata automatically
def version = getWorkflowVersion()
log.info "Pipeline version: ${version}"

// With explicit parameters
def customVersion = getWorkflowVersion('2.0.0', 'abc123')
log.info "Custom version: ${customVersion}"
```

**Integration in Pipeline:**

```nextflow
workflow {
    // Use in logging
    log.info "Starting ${workflow.manifest.name} ${getWorkflowVersion()}"

    // Include in MultiQC reports
    def version_info = [
        'Workflow': getWorkflowVersion(),
        'Nextflow': workflow.nextflow.version
    ]
}
```

---

### `softwareVersionsToYAML()`

**Description:**  
Combines software versions from multiple sources into a single YAML string. This is the recommended function for collecting and formatting version information in pipelines. It accepts mixed input types including YAML strings, file paths, topic tuples, and maps, making it ideal for both legacy and modern pipelines.

**Function Signature:**

```nextflow
// Positional syntax
String softwareVersionsToYAML(Object versions)
String softwareVersionsToYAML(Object versions, Map options)

// Named parameter syntax
String softwareVersionsToYAML(softwareVersions: channel, nextflowVersion: version)
```

**Parameters:**

- `versions` (Object): Can be a Nextflow channel or list containing:
  - YAML strings (legacy inline format)
  - File paths (String), File, or Path objects pointing to `versions.yml`
  - Topic tuples: `[process, tool, version]`
  - Maps of `tool -> version`
  - Any combination of the above
- `options` (Map, optional):
  - `nextflowVersion` (String/VersionNumber): Override the detected Nextflow version

**Returns:**

- `String` (if input is a List): Combined YAML string directly
- `Channel<String>` (if input is a Channel): Channel emitting the combined YAML string

**Usage Examples:**

**Basic usage with a channel:**

```nextflow
include { softwareVersionsToYAML } from 'plugin/nf-core-utils'

workflow {
    // Collect versions from multiple sources
    ch_versions = Channel.empty()

    ch_versions = ch_versions.mix(FASTQC.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS.out.versions)

    // Combine into YAML
    ch_software_versions = softwareVersionsToYAML(ch_versions)

    // Save to file
    ch_software_versions.collectFile(
        name: 'software_versions.yml',
        storeDir: "${params.outdir}/pipeline_info"
    )
}
```

**With custom Nextflow version:**

```nextflow
// Override the detected Nextflow version
ch_software_versions = softwareVersionsToYAML(
    ch_versions,
    nextflowVersion: '24.10.0'
)

// Or using named parameters
ch_software_versions = softwareVersionsToYAML(
    softwareVersions: ch_versions,
    nextflowVersion: workflow.nextflow.version
)
```

**Mixed input types:**

```nextflow
workflow {
    // Mix different version formats
    ch_mixed = Channel.of(
        'fastqc: 0.12.1',                              // YAML string
        ['NFCORE_SAMTOOLS', 'samtools', '1.17'],       // Topic tuple
        [multiqc: '1.15', python: '3.9.0'],            // Map
        file('legacy_versions.yml')                     // File path
    )

    // Function handles all types seamlessly
    ch_versions_yaml = softwareVersionsToYAML(ch_mixed)
}
```

**Direct list processing (no channel):**

```nextflow
workflow.onComplete {
    // Collect versions as a list
    def versionsList = [
        ['FASTQC', 'fastqc', '0.12.1'],
        ['SAMTOOLS', 'samtools', '1.17']
    ]

    // Get YAML string directly
    def versionsYaml = softwareVersionsToYAML(
        versionsList,
        [nextflowVersion: workflow.nextflow.version]
    )

    file("${params.outdir}/pipeline_info/versions.yml").text = versionsYaml
}
```

**Output Format:**

The function generates a structured YAML with processes and tools sorted alphabetically:

```yaml
FASTQC:
  fastqc: 0.12.1
SAMTOOLS:
  samtools: 1.17
Software:
  multiqc: 1.15
  python: 3.9.0
Workflow:
  nf-core/mypipeline: v1.0.0
  Nextflow: 24.10.0
```

**Key Features:**

- **Flexible Input**: Accepts any combination of YAML strings, files, tuples, and maps
- **Process Extraction**: Automatically extracts process names from full paths (e.g., `NFCORE_RNASEQ:TRIMGALORE:FASTQC` → `FASTQC`)
- **Tool Name Cleaning**: Removes prefixes like `tool:` from tool names
- **Alphabetical Sorting**: Sorts both processes and tools for consistent output
- **Graceful Error Handling**: Continues processing even if some entries are invalid
- **Automatic Merging**: Combines versions from the same tool across different sources
- **Workflow Metadata**: Automatically includes workflow name, version, and Nextflow version

**Migration Patterns:**

```nextflow
// Legacy pattern (old versions.yml files)
ch_versions = Channel.fromPath('**/versions.yml')
ch_versions_yaml = softwareVersionsToYAML(ch_versions)

// Modern pattern (topic channels)
ch_versions = PROCESS.out.versions  // [process, tool, version]
ch_versions_yaml = softwareVersionsToYAML(ch_versions)

// Hybrid pattern (transitioning)
ch_legacy = Channel.fromPath('legacy/**/versions.yml')
ch_modern = MODERN_PROCESS.out.versions
ch_all = ch_legacy.mix(ch_modern)
ch_versions_yaml = softwareVersionsToYAML(ch_all)
```

---

### `processMixedVersionSources()`

**Description:**  
Processes version information from both topic channels (modern approach) and traditional files (legacy approach). This is the key function for migration scenarios where both formats coexist.

**Function Signature:**

```nextflow
String processMixedVersionSources(List<List> topicVersions, List<String> versionsFiles)
```

**Parameters:**

- `topicVersions` (List<List>): Topic channel data in `[process, name, version]` format
- `versionsFiles` (List<String>): List of paths to legacy `versions.yml` files

**Returns:**

- `String`: Combined YAML string with all version information

**Usage Example:**

```nextflow
include { processMixedVersionSources } from 'plugin/nf-core-utils'

workflow {
    // Collect from both sources during migration
    def topicVersions = PROCESS_A.out.versions.collect()
    def legacyFiles = ['process_b/versions.yml', 'process_c/versions.yml']

    def allVersions = processMixedVersionSources(topicVersions, legacyFiles)

    // Use in MultiQC
    file("${params.outdir}/pipeline_info/software_versions.yml").text = allVersions
}
```

**Migration Scenario:**

```nextflow
// Stage 2: Mixed approach during pipeline migration
workflow {
    // Modern processes using topic channels
    FASTQC(samples)
    MULTIQC(reports.collect())

    // Legacy processes still using files (being gradually updated)
    LEGACY_PROCESS(data)

    // Combine both approaches seamlessly
    def modernVersions = FASTQC.out.versions
        .mix(MULTIQC.out.versions)
        .collect()

    def legacyFiles = LEGACY_PROCESS.out.versions.collect()

    def combinedVersions = processMixedVersionSources(modernVersions, legacyFiles)
}
```

---

### `convertLegacyYamlToEvalSyntax()`

**Description:**  
Converts traditional YAML version content to the new eval syntax format for topic channels. Useful for migrating existing version data or integrating legacy systems.

**Function Signature:**

```nextflow
List<List> convertLegacyYamlToEvalSyntax(String yamlContent, String processName = 'LEGACY')
```

**Parameters:**

- `yamlContent` (String): YAML content from legacy versions.yml files
- `processName` (String, optional): Process name to assign (defaults to 'LEGACY')

**Returns:**

- `List<List>`: List of `[process, name, version]` tuples ready for topic channel processing

**Usage Example:**

```nextflow
include { convertLegacyYamlToEvalSyntax } from 'plugin/nf-core-utils'

// Convert existing YAML data to topic format
def legacyYaml = '''
FASTQC:
    fastqc: 0.11.9
MULTIQC:
    multiqc: 1.12
    python: 3.9.0
'''

def topicData = convertLegacyYamlToEvalSyntax(legacyYaml, 'CONVERTED_LEGACY')
// Result: [['CONVERTED_LEGACY', 'fastqc', '0.11.9'], ['CONVERTED_LEGACY', 'multiqc', '1.12'], ...]

// Use in topic channels
channel.fromList(topicData).set { converted_versions_ch }
```

**File Migration Example:**

```nextflow
// Migrate existing versions.yml files to topic format
def legacyVersionsFile = file("legacy_versions.yml")
def yamlContent = legacyVersionsFile.text

def migratedVersions = convertLegacyYamlToEvalSyntax(yamlContent, 'MIGRATED_PROCESS')

// Emit to topic channel
channel.fromList(migratedVersions)
    .set { migrated_versions_ch }
```

---

### `generateYamlFromEvalSyntax()`

**Description:**  
Converts topic channel eval syntax data back to YAML format for reporting and MultiQC integration. Essential for generating traditional version reports from modern topic channels.

**Function Signature:**

```nextflow
String generateYamlFromEvalSyntax(List<List> evalData, boolean includeWorkflow = true)
```

**Parameters:**

- `evalData` (List<List>): Topic channel data in `[process, name, version]` format
- `includeWorkflow` (boolean, optional): Include workflow version information (defaults to true)

**Returns:**

- `String`: YAML string suitable for MultiQC and reporting

**Usage Example:**

```nextflow
include { generateYamlFromEvalSyntax } from 'plugin/nf-core-utils'

workflow {
    // Collect topic channel versions
    def topicVersions = [
        ['FASTQC', 'fastqc', '0.11.9'],
        ['MULTIQC', 'multiqc', '1.12'],
        ['MULTIQC', 'python', '3.9.0']
    ]

    // Generate YAML for reporting
    def versionsYaml = generateYamlFromEvalSyntax(topicVersions, true)

    // Write to MultiQC-compatible file
    file("${params.outdir}/pipeline_info/software_versions.yml").text = versionsYaml
}
```

**MultiQC Integration:**

```nextflow
// Generate versions file for MultiQC
workflow {
    // Collect all topic versions
    def allVersions = ch_versions.collect()

    allVersions.map { versions ->
        def yamlContent = generateYamlFromEvalSyntax(versions)
        return tuple('software_versions.yml', yamlContent)
    }.set { ch_versions_yaml }

    // Use in MultiQC process
    MULTIQC(
        multiqc_files.mix(ch_versions_yaml).collect()
    )
}
```

---

### `processVersionsFromTopic()`

**Description:**  
Processes version information exclusively from topic channel format. This is the pure modern approach for pipelines that have fully migrated to topic channels.

**Function Signature:**

```nextflow
String processVersionsFromTopic(List<List> topicData)
```

**Parameters:**

- `topicData` (List<List>): Topic channel data in `[process, name, version]` format

**Returns:**

- `String`: YAML string with processed versions including workflow metadata

**Usage Example:**

```nextflow
include { processVersionsFromTopic } from 'plugin/nf-core-utils'

workflow {
    // Pure topic channel approach
    def versions_ch = channel.empty()

    FASTQC(samples)
    MULTIQC(reports)

    versions_ch = versions_ch.mix(FASTQC.out.versions)
    versions_ch = versions_ch.mix(MULTIQC.out.versions)

    // Process topic versions
    def versionsYaml = processVersionsFromTopic(versions_ch.collect())

    // Output for MultiQC
    file("${params.outdir}/pipeline_info/software_versions.yml").text = versionsYaml
}
```

**Modern Pipeline Pattern:**

```nextflow
// Stage 3: Fully modern topic channel approach
workflow {
    ch_versions = channel.empty()

    // All processes emit to versions topic
    PROCESS_A(input_a)
    PROCESS_B(input_b)
    PROCESS_C(input_c)

    // Collect all versions
    ch_versions = ch_versions
        .mix(PROCESS_A.out.versions)
        .mix(PROCESS_B.out.versions)
        .mix(PROCESS_C.out.versions)

    // Generate final versions report
    ch_versions.collect().map { versions ->
        processVersionsFromTopic(versions)
    }.set { ch_versions_yaml }
}
```

---

### `processVersionsFromFile()`

**Description:**  
Processes version information exclusively from traditional YAML files. This maintains support for legacy pipelines and processes that haven't yet migrated to topic channels.

**Function Signature:**

```nextflow
String processVersionsFromFile(List<String> versionsFiles)
```

**Parameters:**

- `versionsFiles` (List<String>): List of paths to `versions.yml` files

**Returns:**

- `String`: YAML string with processed versions including workflow metadata

**Usage Example:**

```nextflow
include { processVersionsFromFile } from 'plugin/nf-core-utils'

// Legacy approach - still supported
workflow {
    LEGACY_PROCESS_A(input_a)
    LEGACY_PROCESS_B(input_b)

    // Collect file paths
    def versionFiles = [
        'work/process_a/versions.yml',
        'work/process_b/versions.yml'
    ]

    // Process legacy versions
    def versionsYaml = processVersionsFromFile(versionFiles)

    file("${params.outdir}/pipeline_info/software_versions.yml").text = versionsYaml
}
```

**File Collection Pattern:**

```nextflow
// Collect versions from legacy processes
workflow {
    LEGACY_PROCESS(samples)

    // Collect all versions.yml files
    def version_files = LEGACY_PROCESS.out.versions
        .collect()
        .map { files -> files.collect { it.toString() } }

    version_files.map { files ->
        processVersionsFromFile(files)
    }.set { ch_versions_yaml }
}
```

---

### `workflowVersionToChannel()`

**Description:**  
Converts workflow metadata to topic channel format, enabling consistent version handling across all pipeline components. Useful for including workflow version in topic channel processing.

**Function Signature:**

```nextflow
List<List> workflowVersionToChannel()
```

**Parameters:**

- None (uses session context automatically)

**Returns:**

- `List<List>`: Workflow information in `[process, name, version]` format

**Usage Example:**

```nextflow
include { workflowVersionToChannel } from 'plugin/nf-core-utils'

workflow {
    // Get workflow version in topic format
    def workflowVersions = workflowVersionToChannel()

    // Combine with process versions
    def allVersions = ch_process_versions
        .collect()
        .map { processVersions ->
            processVersions + workflowVersions
        }

    // Use in reporting
    allVersions.map { versions ->
        generateYamlFromEvalSyntax(versions)
    }.set { ch_versions_yaml }
}
```

**Complete Topic Channel Integration:**

```nextflow
// Include workflow metadata in topic channel processing
workflow {
    ch_versions = channel.empty()

    // Add workflow version to versions channel
    ch_versions = ch_versions.mix(
        channel.fromList(workflowVersionToChannel())
    )

    // Add process versions
    ch_versions = ch_versions
        .mix(PROCESS_A.out.versions)
        .mix(PROCESS_B.out.versions)

    // Process all versions together
    ch_versions.collect().map { versions ->
        processVersionsFromTopic(versions)
    }.set { ch_comprehensive_versions }
}
```

---

## Topic Channel Migration Guide

### Stage 1: Legacy File-Based Approach

Traditional approach using `versions.yml` files:

```nextflow
// Legacy pattern - still supported
process FASTQC {
    output:
    path "versions.yml", emit: versions

    script:
    """
    echo 'FASTQC:' > versions.yml
    echo '    fastqc: 0.11.9' >> versions.yml
    """
}

workflow {
    FASTQC(samples)

    // Process legacy versions
    def versionsYaml = processVersionsFromFile(FASTQC.out.versions.collect())
}
```

### Stage 2: Mixed Migration Approach

Gradual migration supporting both formats:

```nextflow
// Mixed approach during migration
workflow {
    // Modern processes using topic channels
    MODERN_PROCESS(input_modern)

    // Legacy processes still using files
    LEGACY_PROCESS(input_legacy)

    // Combine both approaches
    def topicVersions = MODERN_PROCESS.out.versions.collect()
    def legacyFiles = LEGACY_PROCESS.out.versions.collect()

    def allVersions = processMixedVersionSources(topicVersions, legacyFiles)

    file("${params.outdir}/pipeline_info/software_versions.yml").text = allVersions
}
```

### Stage 3: Modern Topic Channel Approach

Fully migrated to topic channels:

```nextflow
// Modern pattern - recommended
process FASTQC {
    output:
    tuple val('FASTQC'), val('fastqc'), val('0.11.9'), emit: versions

    script:
    """
    fastqc --version
    """
}

workflow {
    ch_versions = channel.empty()

    FASTQC(samples)
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // Process modern versions
    def versionsYaml = processVersionsFromTopic(ch_versions.collect())
}
```

---

## Best Practices

### 1. Progressive Migration

Migrate modules gradually to avoid disruption:

```nextflow
// Step-by-step migration approach
workflow {
    // Phase 1: Identify high-priority modules
    CRITICAL_PROCESS(input)  // Migrate first

    // Phase 2: Convert supporting modules
    SECONDARY_PROCESS(input)  // Migrate second

    // Phase 3: Handle legacy modules
    LEGACY_PROCESS(input)  // Migrate last

    // Use mixed processing during transition
    def versions = processMixedVersionSources(
        modern_versions.collect(),
        legacy_files.collect()
    )
}
```

### 2. Version Information Consistency

Ensure consistent version extraction across modules:

```nextflow
// Standardized version extraction
process EXAMPLE {
    output:
    tuple val('EXAMPLE'), val('tool_name'), val('${task.container.split(':')[-1]}'), emit: versions

    script:
    """
    tool_name --version | head -1
    """
}
```

### 3. Error Handling

Handle missing or malformed version data gracefully:

```nextflow
// Robust version processing
workflow {
    try {
        def versionsYaml = processVersionsFromTopic(ch_versions.collect())
        file("${params.outdir}/pipeline_info/software_versions.yml").text = versionsYaml
    } catch (Exception e) {
        log.warn "Version processing failed: ${e.message}"
        // Fallback to basic workflow version
        file("${params.outdir}/pipeline_info/software_versions.yml").text = """
Workflow:
    ${workflow.manifest.name}: ${getWorkflowVersion()}
"""
    }
}
```

### 4. Testing Strategies

Test version functions across migration stages:

```groovy
// Test different migration scenarios
test("Legacy version processing") {
    when {
        def legacyFiles = ['test_versions_1.yml', 'test_versions_2.yml']
        def result = processVersionsFromFile(legacyFiles)
    }
    then {
        assert result.contains('FASTQC')
        assert result.contains('fastqc: 0.11.9')
    }
}

test("Topic channel processing") {
    when {
        def topicVersions = [
            ['PROCESS_A', 'tool_a', '1.0.0'],
            ['PROCESS_B', 'tool_b', '2.0.0']
        ]
        def result = processVersionsFromTopic(topicVersions)
    }
    then {
        assert result.contains('PROCESS_A')
        assert result.contains('tool_a: 1.0.0')
    }
}

test("Mixed version processing") {
    when {
        def topicVersions = [['MODERN', 'tool_modern', '3.0.0']]
        def legacyFiles = ['legacy_versions.yml']
        def result = processMixedVersionSources(topicVersions, legacyFiles)
    }
    then {
        assert result.contains('MODERN')
        assert result.contains('tool_modern: 3.0.0')
        // Also contains legacy content
    }
}
```

---

## Integration Examples

### Complete Pipeline Integration

```nextflow
#!/usr/bin/env nextflow

include { getWorkflowVersion; processMixedVersionSources; generateYamlFromEvalSyntax } from 'plugin/nf-core-utils'

workflow {
    // Initialize version tracking
    ch_versions = channel.empty()
    log.info "Starting ${workflow.manifest.name} ${getWorkflowVersion()}"

    // Modern processes with topic channel versions
    FASTQC(samples)
    MULTIQC(reports)

    ch_versions = ch_versions
        .mix(FASTQC.out.versions)
        .mix(MULTIQC.out.versions)

    // Legacy processes (during migration)
    LEGACY_TOOL(data)

    // Combine modern and legacy versions
    def topicVersions = ch_versions.collect()
    def legacyFiles = LEGACY_TOOL.out.versions.collect()

    // Generate comprehensive versions report
    def versionsYaml = processMixedVersionSources(topicVersions, legacyFiles)

    // Write versions file for MultiQC
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    file("software_versions.yml").text = versionsYaml
}
```

### MultiQC Integration

```nextflow
// Full MultiQC integration with version reporting
include { processVersionsFromTopic; workflowVersionToChannel } from 'plugin/nf-core-utils'

workflow {
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // Add workflow metadata to versions
    ch_versions = ch_versions.mix(
        channel.fromList(workflowVersionToChannel())
    )

    // Process modules
    PROCESS_MODULE(samples)
    ch_versions = ch_versions.mix(PROCESS_MODULE.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PROCESS_MODULE.out.reports)

    // Generate versions YAML for MultiQC
    ch_versions.collect().map { versions ->
        def yamlContent = processVersionsFromTopic(versions)
        return tuple('software_versions_mqc.yaml', yamlContent)
    }.set { ch_versions_mqc }

    // Run MultiQC with versions
    MULTIQC(
        ch_multiqc_files.mix(ch_versions_mqc).collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
}
```

---

## Common Issues and Solutions

### Issue: Empty Version Output

**Problem:**

```nextflow
def versionsYaml = processVersionsFromTopic([])
// Result: Only contains workflow metadata
```

**Solution:**

```nextflow
// Verify version channels are populated
ch_versions.view { "Collected version: $it" }

// Ensure proper collection
ch_versions.collect().view { "All versions: $it" }
```

### Issue: Mixed Format Conversion Errors

**Problem:**
Converting between formats produces malformed YAML

**Solution:**

```nextflow
// Validate input formats
def topicVersions = ch_versions.collect().map { versions ->
    versions.findAll { version ->
        version instanceof List && version.size() == 3
    }
}

def versionsYaml = processVersionsFromTopic(topicVersions)
```

### Issue: Legacy File Path Issues

**Problem:**
Version files not found during processing

**Solution:**

```nextflow
// Verify file existence before processing
def validFiles = versionsFiles.findAll { file ->
    def versionFile = file instanceof File ? file : new File(file.toString())
    if (!versionFile.exists()) {
        log.warn "Version file not found: ${versionFile}"
        return false
    }
    return true
}

def versionsYaml = processVersionsFromFile(validFiles)
```

---

## Performance Considerations

### Topic Channel Benefits

- **Reduced I/O**: No file system operations for version collection
- **Better Parallelization**: Versions flow through channels naturally
- **Memory Efficiency**: No intermediate file storage required
- **Simplified Debugging**: Clear data flow through channels

### Migration Performance

```nextflow
// Optimize mixed processing during migration
workflow {
    // Process in parallel where possible
    def topicVersionsFuture = ch_topic_versions.collect()
    def legacyFilesFuture = ch_legacy_files.collect()

    // Combine results efficiently
    tuple(topicVersionsFuture, legacyFilesFuture).map { topic, files ->
        processMixedVersionSources(topic, files)
    }.set { ch_final_versions }
}
```

For more information on version management and topic channels, see the [Nextflow documentation](https://www.nextflow.io/docs/latest/channel.html#topic) and [nf-core best practices](https://nf-co.re/docs/contributing/modules#software-versions).
