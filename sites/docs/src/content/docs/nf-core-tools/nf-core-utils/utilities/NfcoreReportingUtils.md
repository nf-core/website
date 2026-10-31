---
# Copied from docs/utilities/NfcoreReportingUtils.md of nf-core/nf-core-utils 0.7.0 by .github/workflows/add-nf-core-utils-docs.yml: edit it there
title: "MultiQC Integration Tutorial"
shortTitle: "MultiQC Integration Tutorial"
weight: 8
---

# MultiQC Integration Tutorial

MultiQC is the cornerstone of nf-core pipeline reporting, and proper integration requires structured data formatting. This tutorial teaches you how to create professional pipeline reports using the `NfcoreReportingUtils` utility.

## 1. Overview

MultiQC transforms scattered log files and statistics into comprehensive, interactive HTML reports. However, it needs properly formatted input data to work effectively. The `NfcoreReportingUtils` utility bridges this gap by:

- Converting pipeline parameters into MultiQC-compatible YAML format
- Creating structured workflow summaries with metadata
- Generating color-coded terminal logs for different pipeline sections
- Handling complex data structures for professional report presentation

### 1.1. Why Structured Reporting Matters

Raw pipeline outputs are often scattered across multiple files and formats. Users struggle to find:

- Which parameters were used in the analysis
- What software versions were executed
- How long each step took and what resources were used
- Summary statistics and quality metrics

Structured reporting solves these challenges by consolidating information into searchable, interactive dashboards.

## Available Functions

### `paramsSummaryMultiqc()`

**Description:**  
Generates a YAML-formatted string from a map of summary parameters for MultiQC workflow summary integration. This function creates properly structured YAML that MultiQC can parse and display in its workflow summary section.

**Function Signature:**

```groovy
static String paramsSummaryMultiqc(Map summaryParams)
```

**Parameters:**

- `summaryParams` (Map): Nested map of parameter groups and their parameters

**Returns:**

- `String`: YAML-formatted string suitable for MultiQC workflow summary

**Usage Example:**

```nextflow
include { paramsSummaryMultiqc } from 'plugin/nf-core-utils'

workflow {
    // Organize parameters into logical groups
    def summary_params = [
        'Core Nextflow options': [
            'revision': workflow.revision ?: 'N/A',
            'runName': workflow.runName,
            'containerEngine': workflow.containerEngine,
            'profile': workflow.profile
        ],
        'Input/output options': [
            'input': params.input ?: 'N/A',
            'outdir': params.outdir ?: 'N/A',
            'publish_dir_mode': params.publish_dir_mode
        ],
        'Reference genome options': [
            'genome': params.genome ?: 'None',
            'fasta': params.fasta ?: 'N/A'
        ]
    ]

    // Generate MultiQC-compatible YAML
    def yaml_content = paramsSummaryMultiqc(summary_params)

    // Write to file for MultiQC consumption
    file("${params.outdir}/pipeline_info/workflow_summary_mqc.yaml").text = yaml_content
}
```

**MultiQC Integration Pattern:**

```nextflow
process MULTIQC {
    input:
    path('*')

    output:
    path("multiqc_report.html"), emit: report

    script:
    // Include workflow summary in MultiQC
    def summary = [
        'Pipeline': workflow.manifest.name,
        'Version': workflow.manifest.version,
        'Run Name': workflow.runName
    ]
    def yaml = paramsSummaryMultiqc(['Workflow Summary': summary])
    """
    echo '$yaml' > workflow_summary_mqc.yaml
    multiqc -f .
    """
}
```

---

### `workflowSummaryMQC()`

**Description:**  
Creates a comprehensive workflow summary template for MultiQC with both HTML and text formats. This function generates detailed workflow metadata including execution information, parameters, and results.

**Function Signature:**

```groovy
static Map workflowSummaryMQC(Map summary, List<String> nfMetadataList, Map results)
```

**Parameters:**

- `summary` (Map): Map of workflow parameters and metadata
- `nfMetadataList` (List\<String\>): List of metadata fields to include from workflow
- `results` (Map): Map of pipeline results and statistics

**Returns:**

- `Map`: Contains both `html` and `text` formatted summaries for MultiQC

**Usage Example:**

```nextflow
include { workflowSummaryMQC } from 'plugin/nf-core-utils'

workflow {
    // Define summary parameters
    def summary = [
        'Run ID': workflow.runName,
        'User': workflow.userName,
        'Command': workflow.commandLine,
        'Profile': workflow.profile,
        'Container': workflow.containerEngine
    ]

    // Define metadata fields to include
    def metadata_fields = [
        'version',
        'nextflow_version',
        'command_line',
        'start_time',
        'complete_time',
        'duration'
    ]

    // Define results summary
    def results = [
        'Samples Processed': '150',
        'Successfully Completed': '148',
        'Failed QC': '2',
        'Total Runtime': workflow.duration
    ]

    // Generate comprehensive summary
    def mqc_summary = workflowSummaryMQC(summary, metadata_fields, results)

    // Write HTML summary for MultiQC
    file("${params.outdir}/pipeline_info/workflow_summary_mqc.html").text = mqc_summary.html

    // Write text summary for logs
    log.info mqc_summary.text
}
```

**Advanced Integration Example:**

```nextflow
workflow COMPREHENSIVE_REPORTING {
    // Collect pipeline statistics
    def pipeline_stats = [
        'input_files': input_files.size(),
        'processes_run': completed_processes.size(),
        'total_cpu_hours': workflow.stats.computeTimeFmt,
        'peak_memory': workflow.stats.peakRss
    ]

    // Generate detailed workflow summary
    def detailed_summary = [
        'Execution Environment': [
            'workDir': workflow.workDir,
            'launchDir': workflow.launchDir,
            'projectDir': workflow.projectDir
        ],
        'Resource Usage': pipeline_stats,
        'Configuration': [
            'config_files': workflow.configFiles.join(', '),
            'container_engine': workflow.containerEngine
        ]
    ]

    def metadata_fields = [
        'version', 'nextflow_version', 'command_line',
        'start_time', 'complete_time', 'duration',
        'success', 'exit_status'
    ]

    def results_summary = [
        'Total Samples': samples.size(),
        'Successful': successful_samples.size(),
        'Failed': failed_samples.size(),
        'Success Rate': "${(successful_samples.size() / samples.size() * 100).round(1)}%"
    ]

    def comprehensive_summary = workflowSummaryMQC(
        detailed_summary,
        metadata_fields,
        results_summary
    )

    // Use in MultiQC and logging
    publishSummaryReports(comprehensive_summary)
}
```

---

### `sectionLogs()`

**Description:**  
Generates formatted log messages for different pipeline sections with optional color coding. This function creates consistent, readable terminal output for section-based logging.

**Function Signature:**

```groovy
static Map<String, String> sectionLogs(Map sections, boolean monochrome = false)
```

**Parameters:**

- `sections` (Map): Map of section names to log messages
- `monochrome` (boolean, optional): If true, disables color formatting (default: false)

**Returns:**

- `Map<String, String>`: Map of section names to formatted log strings

**Usage Example:**

```nextflow
include { sectionLogs } from 'plugin/nf-core-utils'

workflow {
    // Define section information
    def pipeline_sections = [
        'Input Validation': 'Validated 150 input files successfully',
        'Quality Control': 'FastQC completed for all samples',
        'Preprocessing': 'Trimming and filtering completed',
        'Alignment': 'STAR alignment completed with 95% mapping rate',
        'Analysis': 'Differential expression analysis completed',
        'Reporting': 'MultiQC report generated successfully'
    ]

    // Generate formatted section logs
    def section_logs = sectionLogs(pipeline_sections, params.monochrome_logs)

    // Display section summaries
    log.info "\n" + "=" * 50
    log.info "PIPELINE EXECUTION SUMMARY"
    log.info "=" * 50

    section_logs.each { section, message ->
        log.info message
    }

    log.info "=" * 50
}
```

**Integration with Pipeline Sections:**

```nextflow
workflow SECTIONED_PIPELINE {
    def section_status = [:]

    // Section 1: Input Processing
    VALIDATE_INPUT(input_files)
    section_status['Input Validation'] = "Processed ${input_files.size()} input files"

    // Section 2: Quality Control
    FASTQC(VALIDATE_INPUT.out.validated)
    section_status['Quality Control'] = "QC completed for ${FASTQC.out.html.size()} samples"

    // Section 3: Main Analysis
    MAIN_ANALYSIS(FASTQC.out.passed)
    section_status['Main Analysis'] = "Analysis completed successfully"

    // Generate section summary
    def formatted_logs = sectionLogs(section_status, params.monochrome_logs)

    // Display at completion
    workflow.onComplete {
        log.info "\nPipeline Section Summary:"
        formatted_logs.each { section, log_message ->
            log.info log_message
        }
    }
}
```

**Color-Aware Logging:**

```nextflow
// Conditional color formatting based on success/failure
workflow COLOR_AWARE_LOGGING {
    def sections = [:]

    // Process sections and track success/failure
    try {
        CRITICAL_PROCESS(data)
        sections['Critical Analysis'] = "✓ Successfully completed critical analysis"
    } catch (Exception e) {
        sections['Critical Analysis'] = "✗ Critical analysis failed: ${e.message}"
    }

    try {
        OPTIONAL_PROCESS(data)
        sections['Optional Analysis'] = "✓ Optional analysis completed"
    } catch (Exception e) {
        sections['Optional Analysis'] = "⚠ Optional analysis skipped: ${e.message}"
    }

    // Format with colors (if enabled)
    def colored_logs = sectionLogs(sections, params.monochrome_logs)

    workflow.onComplete {
        colored_logs.each { section, message ->
            if (message.contains('✗')) {
                log.error message
            } else if (message.contains('⚠')) {
                log.warn message
            } else {
                log.info message
            }
        }
    }
}
```

---

## MultiQC Integration Best Practices

### 1. Parameter Organization

Structure parameters logically for better MultiQC display:

```nextflow
// GOOD - Grouped structure
def summary_params = [
    'Core Nextflow options': [
        'runName': workflow.runName,
        'profile': workflow.profile,
        'revision': workflow.revision
    ],
    'Input/output options': [
        'input': params.input,
        'outdir': params.outdir
    ],
    'Analysis parameters': [
        'genome': params.genome,
        'aligner': params.aligner
    ]
]

// AVOID - Flat structure
def summary_params = [
    runName: workflow.runName,
    input: params.input,
    genome: params.genome
]
```

### 2. Null Value Handling

Always provide fallback values for potentially null parameters:

```nextflow
def summary_params = [
    'Core Options': [
        'revision': workflow.revision ?: 'N/A',
        'runName': workflow.runName ?: 'Unknown',
        'profile': workflow.profile ?: 'default'
    ]
]
```

### 3. MultiQC File Integration

Integrate reporting utilities with MultiQC processes:

```nextflow
process MULTIQC {
    input:
    path('*')

    output:
    path("multiqc_report.html"), emit: report

    script:
    def summary = ['Pipeline': workflow.manifest.name]
    def workflow_yaml = paramsSummaryMultiqc(['Workflow': summary])
    """
    echo '$workflow_yaml' > workflow_summary_mqc.yaml
    multiqc -f .
    """
}
```

---

## Integration Examples

### Complete Pipeline Reporting Setup

```nextflow
#!/usr/bin/env nextflow

include {
    paramsSummaryMultiqc;
    workflowSummaryMQC;
    sectionLogs
} from 'plugin/nf-core-utils'

workflow {
    // Initialize section tracking
    def sections = [:]

    // Main pipeline processing
    PREPROCESSING(samples)
    sections['Preprocessing'] = "Completed preprocessing for ${samples.size()} samples"

    ANALYSIS(PREPROCESSING.out.processed)
    sections['Analysis'] = "Analysis completed successfully"

    QUALITY_CONTROL(ANALYSIS.out.results)
    sections['Quality Control'] = "QC metrics generated"

    // Generate comprehensive reporting
    generatePipelineReports(sections)
}

def generatePipelineReports(sections) {
    // 1. Parameter summary for MultiQC
    def param_summary = [
        'Core Nextflow options': [
            'revision': workflow.revision ?: 'N/A',
            'runName': workflow.runName,
            'profile': workflow.profile
        ],
        'Input/output options': [
            'input': params.input,
            'outdir': params.outdir
        ]
    ]

    def param_yaml = paramsSummaryMultiqc(param_summary)
    file("${params.outdir}/pipeline_info/params_summary_mqc.yaml").text = param_yaml

    // 2. Comprehensive workflow summary
    def workflow_summary = [
        'Pipeline': workflow.manifest.name,
        'Version': workflow.manifest.version,
        'Command': workflow.commandLine
    ]

    def metadata_fields = [
        'version', 'nextflow_version', 'start_time',
        'complete_time', 'duration', 'success'
    ]

    def results = [
        'Total Processes': sections.size(),
        'Success Status': workflow.success ? 'Completed' : 'Failed'
    ]

    def mqc_summary = workflowSummaryMQC(workflow_summary, metadata_fields, results)
    file("${params.outdir}/pipeline_info/workflow_summary_mqc.html").text = mqc_summary.html

    // 3. Section logging
    def section_logs = sectionLogs(sections, params.monochrome_logs)

    workflow.onComplete {
        log.info "\nPipeline Execution Summary:"
        log.info "=" * 50
        section_logs.each { section, message ->
            log.info message
        }
        log.info "=" * 50
    }
}
```

### Dynamic Section Tracking

```nextflow
workflow DYNAMIC_REPORTING {
    def section_tracker = new SectionTracker()

    // Track sections dynamically
    section_tracker.start('Input Validation')
    VALIDATE_INPUT(samples)
    section_tracker.complete('Input Validation', "Validated ${samples.size()} samples")

    section_tracker.start('Processing')
    PROCESS_SAMPLES(VALIDATE_INPUT.out.validated)
    section_tracker.complete('Processing', "Processed all samples successfully")

    section_tracker.start('Analysis')
    ANALYZE_RESULTS(PROCESS_SAMPLES.out.processed)
    section_tracker.complete('Analysis', "Analysis completed with ${ANALYZE_RESULTS.out.results.size()} results")

    // Generate final report
    def final_sections = section_tracker.getSections()
    def formatted_logs = sectionLogs(final_sections, params.monochrome_logs)

    workflow.onComplete {
        formatted_logs.each { section, message ->
            log.info message
        }
    }
}

class SectionTracker {
    def sections = [:]
    def start_times = [:]

    void start(String section) {
        start_times[section] = System.currentTimeMillis()
    }

    void complete(String section, String message) {
        def duration = System.currentTimeMillis() - start_times[section]
        sections[section] = "${message} (${duration}ms)"
    }

    Map getSections() {
        return sections
    }
}
```

---

## Testing Strategies

### Unit Testing

```groovy
test("Parameter summary YAML generation") {
    when {
        def params = [
            'Core Options': [
                'runName': 'test_run',
                'profile': 'docker'
            ]
        ]
        def yaml = paramsSummaryMultiqc(params)
    }

    then {
        assert yaml.contains('Core Options:')
        assert yaml.contains('runName: test_run')
        assert yaml.contains('profile: docker')
    }
}

test("Workflow summary generation") {
    when {
        def summary = ['pipeline': 'test']
        def metadata = ['version', 'start_time']
        def results = ['samples': '10']

        def mqc_summary = workflowSummaryMQC(summary, metadata, results)
    }

    then {
        assert mqc_summary.containsKey('html')
        assert mqc_summary.containsKey('text')
        assert mqc_summary.html.contains('pipeline')
        assert mqc_summary.text.contains('test')
    }
}

test("Section log formatting") {
    when {
        def sections = [
            'QC': 'Quality control completed',
            'Analysis': 'Analysis finished successfully'
        ]
        def logs = sectionLogs(sections, false)
    }

    then {
        assert logs.size() == 2
        assert logs['QC'].contains('Quality control completed')
        assert logs['Analysis'].contains('Analysis finished successfully')
    }
}
```

### Integration Testing

```groovy
test("Complete reporting workflow") {
    when {
        // Simulate full reporting workflow
        def params = generateTestParams()
        def summary = generateTestSummary()
        def sections = generateTestSections()

        // Generate all reports
        def param_yaml = paramsSummaryMultiqc(params)
        def workflow_summary = workflowSummaryMQC(summary, ['version'], ['success': true])
        def section_logs = sectionLogs(sections, false)
    }

    then {
        assert param_yaml != null && !param_yaml.isEmpty()
        assert workflow_summary.html != null
        assert workflow_summary.text != null
        assert section_logs.size() > 0
    }
}
```

---

## Common Issues and Solutions

### Issue: YAML Formatting Problems

**Problem:**
Generated YAML has formatting issues for MultiQC

**Solution:**

```nextflow
// Ensure proper parameter structure
def params = [
    'Group Name': [  // Use quoted keys for complex names
        'param1': value1?.toString() ?: 'N/A',  // Handle null values
        'param2': value2?.toString() ?: 'Unknown'
    ]
]

def yaml = paramsSummaryMultiqc(params)
```

### Issue: Missing Workflow Metadata

**Problem:**
Workflow summary missing key information

**Solution:**

```nextflow
// Include comprehensive metadata
def metadata_fields = [
    'version', 'nextflow_version', 'command_line',
    'start_time', 'complete_time', 'duration',
    'success', 'exit_status', 'error_message'
]

// Handle missing metadata gracefully
def safe_summary = workflow_metadata.findAll { k, v -> v != null }
```

### Issue: Color Codes in Log Files

**Problem:**
ANSI color codes appear in log files

**Solution:**

```nextflow
// Use monochrome for file logging
def file_logs = sectionLogs(sections, true)  // monochrome = true
def terminal_logs = sectionLogs(sections, params.monochrome_logs)

// Write clean logs to files
file("${params.outdir}/execution_log.txt").text = file_logs.values().join('\n')
```

---

## Performance Considerations

### Efficient Parameter Processing

```nextflow
// Optimize large parameter sets
workflow LARGE_PARAM_PROCESSING {
    // Process parameters in chunks for very large pipelines
    def param_chunks = large_params.collate(100)

    def yaml_sections = []
    param_chunks.eachWithIndex { chunk, index ->
        def chunk_yaml = paramsSummaryMultiqc(["Section_${index}": chunk])
        yaml_sections.add(chunk_yaml)
    }

    // Combine sections
    def final_yaml = yaml_sections.join('\n---\n')
}
```

### Memory-Efficient Reporting

```nextflow
// Handle large section collections efficiently
workflow MEMORY_EFFICIENT_REPORTING {
    def section_buffer = []
    def max_buffer_size = 50

    processes.eachWithIndex { process, index ->
        section_buffer.add(["Process_${index}": process.status])

        if (section_buffer.size() >= max_buffer_size) {
            def batch_logs = sectionLogs(section_buffer.collectEntries(), params.monochrome_logs)
            processBatchLogs(batch_logs)
            section_buffer.clear()
        }
    }

    // Process remaining sections
    if (!section_buffer.isEmpty()) {
        def final_logs = sectionLogs(section_buffer.collectEntries(), params.monochrome_logs)
        processFinalLogs(final_logs)
    }
}
```

For more information on MultiQC integration and pipeline reporting, see the [MultiQC documentation](https://multiqc.info/) and [nf-core reporting guidelines](https://nf-co.re/docs/developing/components/automated-methods-text).
