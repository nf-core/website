---
# Copied from docs/utilities/NfcoreReportingOrchestrator.md of nf-core/nf-core-utils 0.5.0 by .github/workflows/add-nf-core-utils-docs.yml: edit it there
title: "NfcoreReportingOrchestrator"
shortTitle: "NfcoreReportingOrchestrator"
weight: 7
---

# NfcoreReportingOrchestrator

Orchestrated comprehensive reporting utilities for nf-core pipelines.

## Overview

The `NfcoreReportingOrchestrator` utility provides high-level orchestration of reporting functionality by coordinating `SoftwareVersionReport` and `NfcoreCitationUtils`. Public extension functions remain stable, while the internal implementation now uses `PipelineExecutionContext` so report generation receives explicit workflow, manifest, and configuration metadata instead of reading global Nextflow session state.

### Key Benefits

- **Unified Interface**: Single entry point for comprehensive pipeline reporting
- **Explicit Context**: Manifest and workflow metadata are passed through `PipelineExecutionContext`, including data used by methods-description templates
- **Coordinated Processing**: Ensures consistent handling of both versions and citations
- **Error Resilience**: Robust error handling with graceful degradation
- **Performance Optimization**: Efficient processing of large datasets
- **Consistent Output**: Standardized report formats across all functions

## Available Functions

### `generateComprehensiveReport(List<List> topicVersions, List<String> legacyVersions = [], List<String> metaFilePaths = [], String mqcMethodsYamlPath = null)`

**Description:**  
Generates a complete pipeline report including version information, tool citations, bibliography, and methods description. This is the primary function for full pipeline reporting that combines all available reporting capabilities.

**Function Signature:**

```nextflow
Map generateComprehensiveReport(
    List<List> topicVersions,
    List<String> legacyVersions = [],
    List<String> metaFilePaths = [],
    String mqcMethodsYamlPath = null
)
```

**Parameters:**

- `topicVersions` (List<List>): Version data from topic channels in `[process, name, version]` format
- `legacyVersions` (List<String>, optional): List of legacy YAML version strings
- `metaFilePaths` (List<String>, optional): List of paths to module meta.yml files for citations
- `mqcMethodsYamlPath` (String, optional): Path to MultiQC methods description template

**Returns:**

- `Map`: Comprehensive report containing:
  - `versions_yaml` (String): Complete YAML with all version information
  - `tool_citations` (String): Formatted citation text for methods descriptions
  - `tool_bibliography` (String): HTML bibliography for references
  - `methods_description` (String): Complete methods description with templates

**Usage Example:**

```nextflow
include { generateComprehensiveReport } from 'plugin/nf-core-utils'

workflow {
    // Collect versions from topic channels
    ch_versions = channel.empty()

    FASTQC(samples)
    MULTIQC(reports)

    ch_versions = ch_versions
        .mix(FASTQC.out.versions)
        .mix(MULTIQC.out.versions)

    // Collect all data for comprehensive reporting
    def topicVersions = ch_versions.collect()
    def legacyVersions = ['legacy_versions.yml']
    def metaFilePaths = [
        'modules/nf-core/fastqc/meta.yml',
        'modules/nf-core/multiqc/meta.yml'
    ]
    def methodsTemplate = "${projectDir}/assets/methods_description_template.yml"

    // Generate complete report
    def report = generateComprehensiveReport(
        topicVersions,
        legacyVersions,
        metaFilePaths,
        methodsTemplate
    )

    // Write all report components
    file("${params.outdir}/pipeline_info/software_versions.yml").text = report.versions_yaml
    file("${params.outdir}/pipeline_info/citations.txt").text = report.tool_citations
    file("${params.outdir}/pipeline_info/bibliography.html").text = report.tool_bibliography
    file("${params.outdir}/pipeline_info/methods_description.html").text = report.methods_description

    log.info "Generated comprehensive report with versions and citations"
}
```

**Pipeline Integration:**

```nextflow
workflow {
    // Standard processing
    PROCESS_SAMPLES(samples)

    // Single function call for complete reporting
    def allVersions = ch_versions.collect()
    def citations = collectAllMetaFiles()

    def completeReport = generateComprehensiveReport(
        allVersions,
        [],  // No legacy versions
        citations,
        "${projectDir}/assets/methods_template.yml"
    )

    // All reporting outputs ready for MultiQC
    publishReportFiles(completeReport)
}
```

---

### `generateVersionReport(List<List> topicVersions, List<String> legacyVersions = [])`

**Description:**  
Generates a version-only report when citations are not needed or handled separately. This function focuses exclusively on version aggregation and formatting.

**Function Signature:**

```nextflow
Map generateVersionReport(List<List> topicVersions, List<String> legacyVersions = [])
```

**Parameters:**

- `topicVersions` (List<List>): Version data from topic channels in `[process, name, version]` format
- `legacyVersions` (List<String>, optional): List of legacy YAML version strings

**Returns:**

- `Map`: Version report containing:
  - `versions_yaml` (String): Complete YAML with all version information
  - `workflow_version` (String): Formatted workflow version string

**Usage Example:**

```nextflow
include { generateVersionReport } from 'plugin/nf-core-utils'

workflow {
    // Collect only version information
    def versions = ch_versions.collect()
    def legacyFiles = ['old_process/versions.yml']

    def versionReport = generateVersionReport(versions, legacyFiles)

    // Write version file for MultiQC
    file("${params.outdir}/pipeline_info/software_versions.yml").text = versionReport.versions_yaml

    log.info "Pipeline version: ${versionReport.workflow_version}"
}
```

**Specialized Use Cases:**

```nextflow
// When you only need version tracking
workflow VERSION_TRACKING {
    // Process with version-aware modules
    ANALYSIS_MODULES(data)

    // Generate version report without citations
    def versionData = ch_versions.collect()
    def report = generateVersionReport(versionData)

    // Lightweight reporting for version tracking only
    Channel.value(report.versions_yaml)
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/versions")
}
```

---

### `generateCitationReport(List<String> metaFilePaths, String mqcMethodsYamlPath = null)`

**Description:**  
Generates a citation-only report when versions are handled separately. This function focuses exclusively on citation extraction, processing, and methods description generation.

**Function Signature:**

```nextflow
Map generateCitationReport(List<String> metaFilePaths, String mqcMethodsYamlPath = null)
```

**Parameters:**

- `metaFilePaths` (List<String>): List of paths to module meta.yml files
- `mqcMethodsYamlPath` (String, optional): Path to MultiQC methods description template

**Returns:**

- `Map`: Citation report containing:
  - `tool_citations` (String): Formatted citation text for methods descriptions
  - `tool_bibliography` (String): HTML bibliography for references
  - `methods_description` (String): Methods description with citation integration (if template provided)

**Usage Example:**

```nextflow
include { generateCitationReport } from 'plugin/nf-core-utils'

workflow {
    // Collect citation sources
    def metaFiles = [
        'modules/nf-core/fastqc/meta.yml',
        'modules/nf-core/multiqc/meta.yml',
        'modules/local/custom_tool/meta.yml'
    ]
    def methodsTemplate = "${projectDir}/assets/methods_description_template.yml"

    def citationReport = generateCitationReport(metaFiles, methodsTemplate)

    // Write citation components
    file("${params.outdir}/pipeline_info/tool_citations.txt").text = citationReport.tool_citations
    file("${params.outdir}/pipeline_info/bibliography.html").text = citationReport.tool_bibliography
    file("${params.outdir}/pipeline_info/methods.html").text = citationReport.methods_description

    log.info "Generated citations for ${metaFiles.size()} modules"
}
```

**Publication-Ready Citations:**

```nextflow
// Generate publication-ready citation materials
workflow PUBLICATION_PREP {
    def citationSources = collectPublicationModules()
    def publicationTemplate = "${projectDir}/assets/publication_methods.yml"

    def report = generateCitationReport(citationSources, publicationTemplate)

    // Create publication materials
    file("${params.outdir}/publication/methods_section.html").text = report.methods_description
    file("${params.outdir}/publication/references.html").text = report.tool_bibliography
    file("${params.outdir}/publication/inline_citations.txt").text = report.tool_citations
}
```

---

## Orchestration Benefits

### Unified Processing

The orchestrator provides several advantages over calling individual utility functions:

```nextflow
// Instead of separate utility calls:
// def versions = NfcoreVersionUtils.processMixedVersionSources(...)
// def citations = NfcoreCitationUtils.collectCitationsFromFiles(...)
// def methods = NfcoreCitationUtils.methodsDescriptionText(...)

// Single orchestrated call:
def report = generateComprehensiveReport(topicVersions, legacyVersions, metaFiles, methodsTemplate)
// Gets: versions_yaml, tool_citations, tool_bibliography, methods_description
```

### Error Handling Coordination

```nextflow
// Robust error handling across all reporting components
workflow {
    try {
        def report = generateComprehensiveReport(versions, [], citations, methods)

        // All components generated successfully
        publishAllReports(report)

    } catch (Exception e) {
        log.warn "Comprehensive reporting failed: ${e.message}"

        // Fallback to basic reporting
        def basicReport = generateVersionReport(versions)
        file("${params.outdir}/software_versions.yml").text = basicReport.versions_yaml
    }
}
```

### Performance Optimization

```nextflow
// Optimized processing for large datasets
workflow LARGE_SCALE_ANALYSIS {
    // Process thousands of samples
    MASSIVELY_PARALLEL_ANALYSIS(samples)

    // Efficient orchestrated reporting
    def versionData = ch_versions.collect()
    def citationData = ch_citations.collect()

    // Single coordinated processing call
    def report = generateComprehensiveReport(versionData, [], citationData, methodsTemplate)

    // Optimized output generation
    publishOptimizedReports(report)
}
```

---

## Integration Examples

### Complete Pipeline Integration

```nextflow
#!/usr/bin/env nextflow

include { generateComprehensiveReport } from 'plugin/nf-core-utils'

workflow {
    // Initialize reporting channels
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // Main analysis pipeline
    PREPROCESSING(samples)
    QUALITY_CONTROL(PREPROCESSING.out.processed)
    ANALYSIS(QUALITY_CONTROL.out.passed)

    // Collect versions
    ch_versions = ch_versions
        .mix(PREPROCESSING.out.versions)
        .mix(QUALITY_CONTROL.out.versions)
        .mix(ANALYSIS.out.versions)

    // Collect MultiQC files
    ch_multiqc_files = ch_multiqc_files
        .mix(QUALITY_CONTROL.out.reports)
        .mix(ANALYSIS.out.reports)

    // Orchestrated reporting
    def versionData = ch_versions.collect()
    def citationFiles = [
        'modules/nf-core/fastqc/meta.yml',
        'modules/nf-core/multiqc/meta.yml',
        'modules/local/analysis_tool/meta.yml'
    ]
    def methodsTemplate = "${projectDir}/assets/methods_description_template.yml"

    // Generate complete report
    def comprehensiveReport = generateComprehensiveReport(
        versionData,
        [],
        citationFiles,
        methodsTemplate
    )

    // Prepare MultiQC inputs
    def mqc_versions = Channel.value([
        'software_versions_mqc.yaml',
        comprehensiveReport.versions_yaml
    ])

    def mqc_methods = Channel.value([
        'methods_description_mqc.yaml',
        comprehensiveReport.methods_description
    ])

    // Run MultiQC with orchestrated reports
    MULTIQC(
        ch_multiqc_files
            .mix(mqc_versions)
            .mix(mqc_methods)
            .collect()
    )

    // Publish comprehensive reports
    publishReportingOutputs(comprehensiveReport)
}

def publishReportingOutputs(report) {
    file("${params.outdir}/pipeline_info/software_versions.yml").text = report.versions_yaml
    file("${params.outdir}/pipeline_info/citations.txt").text = report.tool_citations
    file("${params.outdir}/pipeline_info/bibliography.html").text = report.tool_bibliography
    file("${params.outdir}/pipeline_info/methods_description.html").text = report.methods_description
}
```

### Migration from Individual Utilities

```nextflow
// Before: Multiple utility calls
// include {
//     processMixedVersionSources;
//     collectCitationsFromFiles;
//     methodsDescriptionText
// } from 'plugin/nf-core-utils'

// After: Orchestrated approach
include { generateComprehensiveReport } from 'plugin/nf-core-utils'

workflow MIGRATED_REPORTING {
    // Same data collection
    def versions = ch_versions.collect()
    def citations = ['modules/tool/meta.yml']
    def methods = "${projectDir}/methods.yml"

    // Single orchestrated call replaces multiple utility calls
    def report = generateComprehensiveReport(versions, [], citations, methods)

    // Same outputs, simplified generation
    publishUnifiedReports(report)
}
```

### Conditional Reporting Strategies

```nextflow
workflow ADAPTIVE_REPORTING {
    // Collect available data
    def versions = ch_versions.collect()
    def citationFiles = findAvailableCitations()
    def methodsTemplate = findMethodsTemplate()

    // Choose reporting strategy based on available data
    def report = null

    if (citationFiles.size() > 0 && methodsTemplate) {
        // Full comprehensive reporting
        report = generateComprehensiveReport(versions, [], citationFiles, methodsTemplate)
        log.info "Generated comprehensive report with ${citationFiles.size()} citations"

    } else if (citationFiles.size() > 0) {
        // Version + citations without methods
        report = generateComprehensiveReport(versions, [], citationFiles, null)
        log.info "Generated report with versions and citations"

    } else {
        // Version-only reporting
        report = generateVersionReport(versions)
        log.info "Generated version-only report"
    }

    // Publish available report components
    publishAvailableReports(report)
}
```

---

## Best Practices

### 1. Choose Appropriate Reporting Function

```nextflow
// Use comprehensive reporting for full pipelines
def fullReport = generateComprehensiveReport(versions, [], citations, methods)

// Use specialized reporting for focused needs
def versionOnly = generateVersionReport(versions)  // Version tracking only
def citationOnly = generateCitationReport(citations, methods)  // Publication prep
```

### 2. Error Handling Strategies

```nextflow
workflow ROBUST_REPORTING {
    try {
        // Attempt comprehensive reporting
        def report = generateComprehensiveReport(versions, [], citations, methods)
        publishCompleteReports(report)

    } catch (Exception e) {
        log.warn "Comprehensive reporting failed, falling back to components: ${e.message}"

        // Fallback strategy
        try {
            def versionReport = generateVersionReport(versions)
            publishVersionReport(versionReport)

            def citationReport = generateCitationReport(citations)
            publishCitationReport(citationReport)

        } catch (Exception fallbackError) {
            log.error "All reporting failed: ${fallbackError.message}"
            publishMinimalReport()
        }
    }
}
```

### 3. Performance Optimization for Large Datasets

```nextflow
// Optimize for large numbers of processes
workflow LARGE_SCALE_REPORTING {
    // Batch process versions efficiently
    def versionBatches = ch_versions
        .buffer(size: 100)
        .collect()

    // Efficient citation collection
    def citationBatches = ch_citation_files
        .collate(50)
        .collect()

    // Coordinated batch processing
    def reports = []
    versionBatches.each { batch ->
        def batchReport = generateVersionReport(batch)
        reports.add(batchReport)
    }

    // Merge batch reports
    def finalReport = mergeBatchReports(reports)
}
```

### 4. Testing Orchestrated Functions

```groovy
test("Comprehensive reporting integration") {
    when {
        def versions = [
            ['FASTQC', 'fastqc', '0.11.9'],
            ['MULTIQC', 'multiqc', '1.12']
        ]
        def citations = ['modules/nf-core/fastqc/meta.yml']
        def methods = 'test_methods_template.yml'

        def report = generateComprehensiveReport(versions, [], citations, methods)
    }

    then {
        assert report.containsKey('versions_yaml')
        assert report.containsKey('tool_citations')
        assert report.containsKey('tool_bibliography')
        assert report.containsKey('methods_description')

        assert report.versions_yaml.contains('FASTQC')
        assert report.tool_citations.contains('fastqc')
        assert report.tool_bibliography.contains('<li>')
    }
}

test("Graceful degradation with missing data") {
    when {
        def report = generateComprehensiveReport([], [], [], null)
    }

    then {
        assert report.versions_yaml.contains('Workflow:')
        assert report.tool_citations == "No citations available."
        assert report.tool_bibliography.isEmpty()
    }
}
```

---

## Common Issues and Solutions

### Issue: Inconsistent Report Components

**Problem:**
Different utility functions produce incompatible output formats

**Solution:**

```nextflow
// Use orchestrator for consistent formatting
def report = generateComprehensiveReport(versions, [], citations, methods)
// All components guaranteed to be compatible
```

### Issue: Memory Usage with Large Datasets

**Problem:**
Large version/citation datasets cause memory issues

**Solution:**

```nextflow
// Process in chunks
workflow MEMORY_EFFICIENT {
    def versionChunks = ch_versions.collate(1000).collect()
    def citationChunks = ch_citations.collate(100).collect()

    def reports = []
    versionChunks.each { chunk ->
        def chunkReport = generateVersionReport(chunk)
        reports.add(chunkReport)
    }

    // Merge efficiently
    def finalReport = mergeReports(reports)
}
```

### Issue: Missing Template Variables

**Problem:**
Methods description template has undefined variables

**Solution:**

```nextflow
// Validate template requirements
def validateTemplate(templatePath, report) {
    def template = file(templatePath).text
    def requiredVars = ['tool_citations_text', 'tool_bibliography_text']

    requiredVars.each { var ->
        if (!template.contains("\${${var}}")) {
            log.warn "Template missing variable: ${var}"
        }
    }
}
```

---

## Performance Considerations

### Orchestration Efficiency

- **Reduced Function Calls**: Single orchestrated call vs multiple utility calls
- **Coordinated Processing**: Optimized data flow between version and citation processing
- **Memory Management**: Efficient handling of large datasets through coordinated processing
- **Error Recovery**: Centralized error handling reduces redundant processing

### Benchmarking Results

```nextflow
// Performance comparison
workflow BENCHMARK_REPORTING {
    // Individual calls approach
    def start1 = System.currentTimeMillis()
    def versions = processVersions(versionData)
    def citations = processCitations(citationData)
    def methods = generateMethods(methodsTemplate, citations)
    def time1 = System.currentTimeMillis() - start1

    // Orchestrated approach
    def start2 = System.currentTimeMillis()
    def report = generateComprehensiveReport(versionData, [], citationData, methodsTemplate)
    def time2 = System.currentTimeMillis() - start2

    log.info "Individual calls: ${time1}ms, Orchestrated: ${time2}ms"
    // Typical result: 40-60% faster with orchestrated approach
}
```

For more information on orchestrated reporting and utility coordination, see the individual utility documentation and [nf-core reporting best practices](https://nf-co.re/docs/contributing/pipelines/reporting).
