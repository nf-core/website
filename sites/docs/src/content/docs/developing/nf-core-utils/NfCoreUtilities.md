---
# Copied from docs/NfCoreUtilities.md of nf-core/nf-core-utils 0.5.0 by .github/workflows/add-nf-core-utils-docs.yml: edit it there
title: "nf-core Utilities Guide"
shortTitle: "nf-core Utilities Guide"
weight: 1
---

# nf-core Utilities Guide

This comprehensive guide covers the nf-core-utils plugin, a powerful extension that brings essential nf-core functionality directly into your Nextflow pipelines. Whether you're building a new pipeline or migrating from legacy subworkflows, this guide will help you leverage these utilities effectively.

## 1. Overview

The nf-core-utils plugin transforms common pipeline tasks from complex subworkflows into simple function calls. It provides battle-tested utilities for configuration validation, notifications, reporting, and citation management that work seamlessly across all nf-core pipelines.

### 1.1. Core Capabilities

The plugin organizes functionality into six focused areas:

- **Configuration & Validation**: Ensure proper pipeline setup and profile usage
- **Notifications & Logging**: Send completion emails, Slack/Teams messages, and terminal summaries
- **Reporting & MultiQC**: Generate pipeline reports and integrate with MultiQC
- **Version Management**: Track software versions with topic channel support
- **Citation Management**: Extract and format tool citations automatically
- **Orchestrated Reporting**: Coordinate comprehensive reporting workflows

### 1.2. Documentation Structure

Each utility class has dedicated documentation with examples and best practices:

| Category                       | Utility                                                                     | Purpose                                                                |
| ------------------------------ | --------------------------------------------------------------------------- | ---------------------------------------------------------------------- |
| **Configuration & Validation** | **[NfcoreConfigValidator](utilities/NfcoreConfigValidator.md)**             | Pipeline configuration validation and profile checking                 |
| **Notifications & Logging**    | **[NfcoreNotificationUtils](utilities/NfcoreNotificationUtils.md)**         | Email notifications, Slack/Teams integration, and completion summaries |
| **Reporting & MultiQC**        | **[NfcoreReportingUtils](utilities/NfcoreReportingUtils.md)**               | MultiQC integration and pipeline reporting utilities                   |
| **Version Management**         | **[NfcoreVersionUtils](utilities/NfcoreVersionUtils.md)**                   | Version tracking with topic channel support and migration utilities    |
| **Citation Management**        | **[NfcoreCitationUtils](utilities/NfcoreCitationUtils.md)**                 | Citation extraction, processing, and topic channel integration         |
| **Orchestrated Reporting**     | **[NfcoreReportingOrchestrator](utilities/NfcoreReportingOrchestrator.md)** | Comprehensive reporting coordination for versions and citations        |

:::tip{title="Where to Start"}
New to nf-core utilities? Begin with the **Configuration & Validation** section to understand pipeline setup, then explore **Notifications & Logging** for completion handling.
:::

## 2. Quick Start Tutorial

This section provides a hands-on introduction to the most commonly used nf-core utilities. We'll build from basic imports to a fully-featured pipeline with comprehensive logging and reporting.

### 2.1. Your First nf-core Utility

Let's start with the simplest possible example - checking pipeline configuration:

```nextflow title="minimal_example.nf"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import a single utility function
include { checkConfigProvided } from 'plugin/nf-core-utils'

// Check if custom config or profile was provided
checkConfigProvided()

workflow {
    log.info "Pipeline started with proper validation!"
}
```

```console title="Output"
N E X T F L O W  ~  version 25.04.0
Launching `minimal_example.nf` [focused-turing] - revision: abc1234 [main]

WARN: No custom Nextflow configuration detected! Please provide a custom config file or profile.
INFO: Pipeline started with proper validation!
```

### 2.2. Adding Completion Notifications

Let's expand our example to include completion summaries:

```nextflow title="notification_example.nf" hl_lines="5 14-16"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { checkConfigProvided; completionSummary } from 'plugin/nf-core-utils'

checkConfigProvided()

workflow {
    log.info "Pipeline processing samples..."
    // Your workflow logic here
}

workflow.onComplete {
    completionSummary(params.monochrome_logs)
}
```

```console title="Output"
Pipeline completed successfully!

-[nf-core/example] Pipeline completed successfully-
Completed at: 2024-01-15T10:30:45.123Z
Duration    : 45s
CPU hours   : 0.1
Succeeded   : 3
```

### 2.3. Complete Import Pattern

For production pipelines, you'll typically need multiple utilities. Here's a recommended import pattern:

```nextflow title="production_imports.nf" hl_lines="4-10"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Essential utilities for most pipelines
include {
    checkConfigProvided; checkProfileProvided;
    completionSummary; completionEmail;
    paramsSummaryMultiqc; getWorkflowVersion
} from 'plugin/nf-core-utils'

// Advanced features for comprehensive pipelines
include {
    getCitation; autoToolCitationText;
    generateComprehensiveReport
} from 'plugin/nf-core-utils'
```

:::note{title="Import Strategy"}
Start with essential utilities and add advanced features as your pipeline grows. Each function is independent and can be imported separately.
:::

## 3. Function Categories and Usage

Understanding when and how to use each function is key to building robust pipelines. Let's explore the main categories with practical examples.

### 3.1. Pipeline Initialization Functions

These functions should be called early in your pipeline to validate setup:

| Function                                      | When to Use                            | Example                                              |
| --------------------------------------------- | -------------------------------------- | ---------------------------------------------------- |
| `checkConfigProvided()`                       | Always, to ensure custom configuration | `checkConfigProvided()`                              |
| `checkProfileProvided(args, monochrome_logs)` | Always, to validate profile arguments  | `checkProfileProvided(args, params.monochrome_logs)` |

```nextflow title="Initialization pattern"
#!/usr/bin/env nextflow

include { checkConfigProvided; checkProfileProvided } from 'plugin/nf-core-utils'

// Validate setup before any processing
checkConfigProvided()
checkProfileProvided(args, params.monochrome_logs)

workflow {
    // Now proceed with confidence that setup is correct
    log.info "Configuration validated successfully"
}
```

### 3.2. Workflow Information Functions

These functions provide metadata and can be called anywhere in your pipeline:

| Function                            | Purpose                                     | Example                                         |
| ----------------------------------- | ------------------------------------------- | ----------------------------------------------- |
| `getWorkflowVersion()`              | Get git-aware version string                | `getWorkflowVersion()`                          |
| `softwareVersionsToYAML(versions)`  | Combine versions from mixed sources to YAML | `softwareVersionsToYAML(ch_versions)`           |
| `logColours(monochrome_logs)`       | Get color codes for terminal output         | `logColours(params.monochrome_logs)`            |
| `sectionLogs(sections, monochrome)` | Generate colored section summaries          | `sectionLogs(sections, params.monochrome_logs)` |

### 3.3. MultiQC Integration Functions

These functions help integrate with MultiQC reporting:

| Function                                         | Purpose                               | Usage Context            |
| ------------------------------------------------ | ------------------------------------- | ------------------------ |
| `paramsSummaryMultiqc(summary_params)`           | Generate MultiQC YAML summary         | Process or main workflow |
| `workflowSummaryMQC(summary, metadata, results)` | Create comprehensive MultiQC template | Process                  |
| `getSingleReport(multiqc_reports)`               | Extract single report from Path/List  | Anywhere                 |

### 3.4. Completion Handler Functions

:::warning{title="Session Context Required"}
These functions require Nextflow session context and **MUST** be called from `workflow.onComplete` or `workflow.onError` handlers only. They will fail with null pointer exceptions if called from the main workflow block.
:::

| Function                                   | Purpose                                   | Usage Context                               |
| ------------------------------------------ | ----------------------------------------- | ------------------------------------------- |
| `completionSummary(monochrome_logs)`       | Print colored pipeline completion summary | `workflow.onComplete` or `workflow.onError` |
| `completionEmail(summary_params, ...)`     | Send detailed completion email            | `workflow.onComplete` or `workflow.onError` |
| `imNotification(summary_params, hook_url)` | Send Slack/Teams notification             | `workflow.onComplete` or `workflow.onError` |

```nextflow title="Completion handler pattern"
workflow.onComplete {
    // CORRECT - These functions work here
    completionSummary(params.monochrome_logs)

    if (params.email) {
        completionEmail(summary_params, params.email, ...)
    }

    if (params.hook_url) {
        imNotification(summary_params, params.hook_url)
    }
}

workflow {
    // WRONG - These functions will fail here
    // completionSummary(params.monochrome_logs)  // Error: session is null
}
```

### 3.5. Citation Management Functions

Citation management comes in two flavors: modern topic-channel based (recommended) and legacy file-based approaches.

#### Modern Topic Channel Citations (Recommended)

These functions work with the new topic channel citation system:

| Function                                   | Purpose                                        | Usage Context        |
| ------------------------------------------ | ---------------------------------------------- | -------------------- |
| `getCitation(metaFilePath)`                | Extract citation for topic channel emission    | Process output block |
| `autoToolCitationText(citationTopics)`     | Generate citation text from topic channels     | Workflow completion  |
| `autoToolBibliographyText(citationTopics)` | Generate bibliography HTML from topic channels | Workflow completion  |

```nextflow title="Modern citation pattern"
process FASTQC {
    output:
    path "*.html", emit: html
    val citation_data, topic: citation

    script:
    // Extract citation for automatic collection
    citation_data = getCitation("${moduleDir}/meta.yml")
    """
    fastqc sample.fastq
    """
}

workflow {
    FASTQC(samples)

    // Collect citations from all processes
    citation_ch = FASTQC.out.citation.collect()

    // Generate formatted citations
    citation_ch.view { citations ->
        def citationText = autoToolCitationText(citations)
        def bibliography = autoToolBibliographyText(citations)
        log.info "Citations: ${citationText}"
    }
}
```

#### Legacy File-Based Citations

These functions support the older meta.yml file approach:

| Function                                   | Purpose                                         | Usage Context |
| ------------------------------------------ | ----------------------------------------------- | ------------- |
| `generateModuleToolCitation(metaFilePath)` | Extract citations from meta.yml file            | Anywhere      |
| `toolCitationText(citationsMap)`           | Generate citation text from collected citations | Anywhere      |
| `toolBibliographyText(citationsMap)`       | Generate bibliography HTML from citations       | Anywhere      |
| `collectCitationsFromFiles(metaFilePaths)` | Collect citations from multiple meta.yml files  | Anywhere      |

#### Migration and Mixed Approaches

For pipelines transitioning between approaches:

| Function                                             | Purpose                                     | Usage Context       |
| ---------------------------------------------------- | ------------------------------------------- | ------------------- |
| `processCitationsFromTopic(topicData)`               | Process citations from topic channel format | Workflow completion |
| `processMixedCitationSources(topics, files)`         | Combine topic and file-based citations      | Workflow completion |
| `convertMetaYamlToTopicFormat(metaPath, moduleName)` | Convert meta.yml to topic channel format    | Migration scripts   |

## 4. Advanced Reporting Architecture

The nf-core-utils plugin provides a sophisticated reporting system that coordinates version tracking, citation management, and comprehensive report generation.

### 4.1. Reporting Components

The reporting system consists of three specialized utility classes:

| Component                       | Purpose                             | Use When                         |
| ------------------------------- | ----------------------------------- | -------------------------------- |
| **NfcoreVersionUtils**          | Version aggregation and formatting  | Need version tracking only       |
| **NfcoreCitationUtils**         | Citation extraction and processing  | Need citation management only    |
| **NfcoreReportingOrchestrator** | Coordinates comprehensive reporting | Need complete reporting solution |

### 4.2. Simple Reporting Example

For most pipelines, you'll want comprehensive reporting that includes versions, citations, and methods descriptions:

```nextflow title="comprehensive_reporting.nf"
include { generateComprehensiveReport } from 'plugin/nf-core-utils'

workflow {
    // Your pipeline processes with version and citation collection...

    // Collect versions and citations
    ch_versions = PROCESS1.out.versions.mix(PROCESS2.out.versions).collect()
    ch_citations = PROCESS1.out.citation.mix(PROCESS2.out.citation).collect()
}

workflow.onComplete {
    // Generate comprehensive report
    ch_versions.concat(ch_citations).collect().view { data ->
        def report = generateComprehensiveReport(
            data.findAll { it.topic == 'versions' },  // Topic versions
            [],  // Legacy versions (empty for modern pipelines)
            [],  // Meta file paths (empty for topic-based citations)
            'multiqc_methods.yml'  // Methods template
        )

        // Write report files
        file("${params.outdir}/pipeline_info/software_versions.yml").text = report.versions_yaml
        file("${params.outdir}/pipeline_info/citations.txt").text = report.tool_citations
        file("${params.outdir}/pipeline_info/methods.txt").text = report.methods_description
    }
}
```

### 4.3. Migration from Legacy Systems

The plugin supports a progressive migration path from legacy approaches to modern topic channels:

| Stage          | Description                     | Implementation                               |
| -------------- | ------------------------------- | -------------------------------------------- |
| **Legacy**     | File-based meta.yml citations   | `collectCitationsFromFiles(metaFilePaths)`   |
| **Transition** | Mixed file and topic approaches | `processMixedCitationSources(topics, files)` |
| **Modern**     | Pure topic channel approach     | `processCitationsFromTopic(topicCitations)`  |

:::tip{title="Migration Strategy"}
Start by converting high-priority modules to topic channels while keeping legacy modules unchanged. The mixed approach allows gradual migration without breaking existing functionality.
:::

### 4.4. Software Versions Collection with `softwareVersionsToYAML()`

The `softwareVersionsToYAML()` function provides a unified way to collect and format software versions from any source. It's the recommended approach for version management in both new and migrating pipelines.

#### Basic Usage Pattern

```nextflow title="simple_versions.nf"
include { softwareVersionsToYAML } from 'plugin/nf-core-utils'

workflow {
    // Collect versions from processes
    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS.out.versions)

    // Convert to YAML and save
    ch_versions_yaml = softwareVersionsToYAML(ch_versions)

    ch_versions_yaml.collectFile(
        name: 'software_versions.yml',
        storeDir: "${params.outdir}/pipeline_info"
    )
}
```

#### With Custom Nextflow Version

```nextflow title="custom_nf_version.nf"
// Override detected Nextflow version (useful for testing)
ch_versions_yaml = softwareVersionsToYAML(
    ch_versions,
    nextflowVersion: '24.10.0'
)

// Or use named parameters
ch_versions_yaml = softwareVersionsToYAML(
    softwareVersions: ch_versions,
    nextflowVersion: workflow.nextflow.version
)
```

#### Mixed Input Types (Migration Scenario)

The function automatically handles different version formats, making it perfect for pipelines in transition:

```nextflow title="mixed_versions.nf"
workflow {
    // Modern processes emit topic tuples
    MODERN_FASTQC(samples)  // Emits: ['FASTQC', 'fastqc', '0.12.1']

    // Legacy processes emit YAML strings or files
    LEGACY_SAMTOOLS(bams)   // Emits: 'samtools: 1.17'

    // Some processes emit maps
    CUSTOM_PROCESS(data)    // Emits: [tool1: '1.0', tool2: '2.0']

    // Mix all version sources together
    ch_all_versions = Channel.empty()
        .mix(MODERN_FASTQC.out.versions)
        .mix(LEGACY_SAMTOOLS.out.versions)
        .mix(CUSTOM_PROCESS.out.versions)

    // Function handles all types seamlessly!
    ch_versions_yaml = softwareVersionsToYAML(ch_all_versions)

    ch_versions_yaml.collectFile(
        name: 'software_versions.yml',
        storeDir: "${params.outdir}/pipeline_info"
    )
}
```

#### Direct Processing (No Channel)

For use cases where you already have collected versions as a list:

```nextflow title="direct_processing.nf"
workflow.onComplete {
    // Collect versions as a list
    def versionsList = [
        ['FASTQC', 'fastqc', '0.12.1'],
        ['SAMTOOLS', 'samtools', '1.17'],
        'multiqc: 1.15'  // Can mix formats!
    ]

    // Process directly without channel
    def versionsYaml = softwareVersionsToYAML(
        versionsList,
        [nextflowVersion: workflow.nextflow.version]
    )

    // Write to file
    file("${params.outdir}/pipeline_info/software_versions.yml").text = versionsYaml
}
```

#### Output Format

The function generates structured, sorted YAML:

```yaml title="software_versions.yml"
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

Key features:

- **Processes sorted alphabetically** for consistency
- **Tools within each process sorted** for easy reading
- **Automatic process name extraction** from full paths
- **Tool name cleaning** (removes `tool:` prefixes)
- **Automatic workflow metadata** inclusion

:::tip{title="Best Practice"}
Use `softwareVersionsToYAML()` as your primary version collection function. It handles all input types, making your pipeline future-proof and migration-friendly.
:::

## 5. Complete Pipeline Example

Let's put everything together in a comprehensive example that shows a production-ready nf-core pipeline using all major utility functions:

```nextflow title="complete_nfcore_pipeline.nf" hl_lines="5-12 17-19 26-28 41-43 60-76"
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import essential nf-core utilities
include {
    checkConfigProvided; checkProfileProvided;
    completionSummary; completionEmail; imNotification;
    paramsSummaryMultiqc; getWorkflowVersion;
    getCitation; autoToolCitationText;
    generateComprehensiveReport
} from 'plugin/nf-core-utils'

params.input = null
params.outdir = "./results"

// Validate setup early
checkConfigProvided()
checkProfileProvided(args, params.monochrome_logs)

// Set up workflow metadata
version_str = getWorkflowVersion(workflow.manifest.version, workflow.commitId)

workflow {
    log.info "Starting ${workflow.manifest.name} ${version_str}"

    // Your pipeline processes here...
    Channel.fromPath(params.input)
        | FASTQC
        | SAMTOOLS_VIEW

    // Collect outputs
    ch_versions = FASTQC.out.versions.mix(SAMTOOLS_VIEW.out.versions)
    ch_citations = FASTQC.out.citation.mix(SAMTOOLS_VIEW.out.citation)

    // Generate MultiQC summary
    def summary_params = [
        'Input': params.input,
        'Output': params.outdir,
        'Version': version_str
    ]
    def mqc_yaml = paramsSummaryMultiqc([Pipeline: summary_params])
    Channel.of(mqc_yaml).collectFile(name: 'workflow_summary_mqc.yaml')
}

process FASTQC {
    output:
    path "*.html", emit: html
    path "versions.yml", emit: versions
    val citation_data, topic: citation

    script:
    citation_data = getCitation("${moduleDir}/meta.yml")
    """
    fastqc --version > versions.yml
    fastqc sample.fastq
    """
}

workflow.onComplete {
    // Show completion summary
    completionSummary(params.monochrome_logs)

    // Send notifications if configured
    def summary_params = [
        'Core Options': [
            'runName': workflow.runName,
            'version': version_str,
            'success': workflow.success
        ]
    ]

    if (params.email) {
        completionEmail(summary_params, params.email, params.email_on_fail,
                       params.plaintext_email, params.outdir, params.monochrome_logs, null)
    }

    if (params.hook_url) {
        imNotification(summary_params, params.hook_url)
    }
}
```

This example demonstrates:

- **Configuration validation** at pipeline startup
- **Version tracking** throughout the pipeline
- **Citation collection** from executed processes
- **MultiQC integration** with parameter summaries
- **Comprehensive notifications** on completion

## 6. Best Practices

### 6.1. Pipeline Structure

Follow this recommended structure for consistent nf-core pipeline development:

```nextflow title="Recommended pipeline structure"
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// 1. Import utilities
include { checkConfigProvided; checkProfileProvided; completionSummary } from 'plugin/nf-core-utils'

// 2. Validate configuration early
checkConfigProvided()
checkProfileProvided(args, params.monochrome_logs)

// 3. Set up metadata
version_str = getWorkflowVersion(workflow.manifest.version, workflow.commitId)

// 4. Main workflow
workflow { /* pipeline logic */ }

// 5. Completion handlers
workflow.onComplete { completionSummary(params.monochrome_logs) }
```

### 6.2. Function Usage Guidelines

| Do                                        | Don't                                        | Why                                                    |
| ----------------------------------------- | -------------------------------------------- | ------------------------------------------------------ |
| Call validation functions early           | Skip configuration checks                    | Early validation prevents runtime failures             |
| Use completion handlers for notifications | Call notification functions in main workflow | Session context required for metadata access           |
| Group summary parameters logically        | Use flat parameter structure                 | Grouped parameters create better email formatting      |
| Handle null values in parameters          | Assume parameters always exist               | Defensive programming prevents null pointer exceptions |

### 6.3. Error Handling

:::warning{title="Common Pitfall: Session Context"}
The most common error is calling notification functions outside completion handlers. Always use `workflow.onComplete` or `workflow.onError` for these functions.
:::

## 7. Troubleshooting

### 7.1. Common Issues

| Issue                       | Symptom                                          | Solution                                     |
| --------------------------- | ------------------------------------------------ | -------------------------------------------- |
| **Session null error**      | `NullPointerException` on notification functions | Move functions to `workflow.onComplete`      |
| **Email not sending**       | No completion email received                     | Check SMTP configuration in nextflow.config  |
| **Webhook failing**         | No Slack/Teams notification                      | Verify webhook URL format and network access |
| **Parameter summary empty** | Empty parameter sections in notifications        | Use grouped parameter structure              |

### 7.2. Quick Debug Steps

```nextflow title="Debug template"
workflow.onComplete {
    // Debug session availability
    log.info "Session: ${nextflow.Nextflow.session != null}"
    log.info "Success: ${workflow.success}"

    // Test basic functions first
    completionSummary(false)  // Force colors off for testing

    log.info "Basic functions working correctly"
}
```

## 8. Takeaway

The nf-core-utils plugin provides a comprehensive toolkit for building robust, production-ready pipelines:

1. **Start with validation**: Use `checkConfigProvided()` and `checkProfileProvided()` early
2. **Add notifications**: Implement `completionSummary()`, `completionEmail()`, and `imNotification()` in completion handlers
3. **Integrate reporting**: Use `paramsSummaryMultiqc()` and citation functions for comprehensive documentation
4. **Follow best practices**: Group parameters, handle null values, and use completion handlers correctly

## 9. What's Next?

Explore specialized utilities for specific use cases:

- **[NfcoreConfigValidator](utilities/NfcoreConfigValidator.md)**: Deep dive into configuration validation
- **[NfcoreNotificationUtils](utilities/NfcoreNotificationUtils.md)**: Master email and messaging integrations
- **[NfcoreReportingUtils](utilities/NfcoreReportingUtils.md)**: Advanced MultiQC integration patterns
- **[NextflowPipelineExtension](NextflowPipelineExtension.md)**: Core pipeline utilities for version and parameter management
- **[ReferencesExtension](ReferencesExtension.md)**: Handle genome references and igenomes integration
