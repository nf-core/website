---
# Copied from docs/utilities/NfcoreCitationUtils.md of nf-core/nf-core-utils 0.6.0 by .github/workflows/add-nf-core-utils-docs.yml: edit it there
title: "NfcoreCitationUtils"
shortTitle: "NfcoreCitationUtils"
weight: 4
md_github_url: https://github.com/nf-core/nf-core-utils/blob/main/docs/utilities/NfcoreCitationUtils.md
---

# NfcoreCitationUtils

Citation management and topic channel utilities for nf-core pipelines.

## Overview

The `NfcoreCitationUtils` utility provides comprehensive citation management functionality for nf-core pipelines, supporting both traditional file-based approaches and modern topic channel patterns. This utility enables progressive migration from legacy meta.yml citation files to the more efficient and accurate topic channel system.

### Topic Channel Migration Strategy for Citations

Similar to version management, the citation system supports a three-phase migration approach:

1. **Legacy Stage**: Traditional meta.yml file processing with manual collection
2. **Transition Stage**: Mixed approach supporting both files and topic channels
3. **Modern Stage**: Pure topic channel approach with automatic runtime citation collection

This migration strategy ensures backward compatibility while enabling pipelines to benefit from runtime-accurate citations where only tools that actually execute are included in the final reports.

## Available Functions

### Citations on the Fly (Recommended)

The lowest-friction way to make citations reflect the tools a pipeline actually ran — generated _on the fly_ from the run itself, not from a hand-maintained list. It requires **no per-module changes**: every nf-core module already emits a `[process, tool, version]` tuple to the `versions` topic when it executes, so that topic _is_ the list of tools that ran. These functions intersect that list with citation metadata parsed from each module's `meta.yml`, so only tools that actually executed are cited.

Two functions cover the flow:

| Function                                          | Purpose                                                 | When to call        |
| ------------------------------------------------- | ------------------------------------------------------- | ------------------- |
| `toolsFromVersionsTopic(topicVersions)`           | Reduce collected `versions` data to the tool names used | Workflow completion |
| `citationsOnTheFly(topicVersions, metaFilePaths)` | Citations for only the tools that ran, from `meta.yml`  | Workflow completion |

---

#### `citationsOnTheFly()`

`citationsOnTheFly(List topicVersions, List<String> metaFilePaths)`

**Description:**
Builds a citations map for **only the tools that ran**, by intersecting the `versions` topic (what executed) with citations parsed from the supplied `meta.yml` files (the citation source). The returned map plugs directly into `toolCitationText()`, `toolBibliographyText()`, and `methodsDescriptionText()`.

Each **short citation** is formatted for direct copy-paste into a methods paragraph — `tool (version, Author et al. year)` — with the version taken from the `versions` topic and the author/year from `meta.yml`. The reference segment degrades gracefully: short author + year → `doi: …` → homepage url → description. So richer `meta.yml` files produce richer citations:

| `meta.yml` fields present   | Short citation                                       |
| --------------------------- | ---------------------------------------------------- |
| `author`, `year`            | `fastqc (0.12.1, Andrews 2010)`                      |
| `author` (multiple), `year` | `samtools (1.21, Danecek et al. 2021)`               |
| `doi` only                  | `multiqc (1.21, doi: 10.1093/bioinformatics/btw354)` |
| `homepage` only             | `seqtk (1.4, https://github.com/lh3/seqtk)`          |

The full `bibliography` entry (author, year, title, journal, doi, url) is unchanged.

**Function Signature:**

```nextflow
Map citationsOnTheFly(List topicVersions, List<String> metaFilePaths)
```

**Parameters:**

- `topicVersions` (List): Collected `versions` topic data — `[process, tool, version]` tuples. Collect with `.collect(flat: false)` so the tuples are preserved (plain `.collect()` flattens them, and the function logs a warning if it receives flattened data).
- `metaFilePaths` (List<String>): Paths to module `meta.yml` files to resolve citation metadata from.
- `extraTools` (List<String>, optional): Tool names to cite even though they never emit a version — e.g. a MultiQC plugin like `multiqcsav`. They're resolved from `meta.yml` like everything else (with no version segment), so you can still add a tool by hand: `citationsOnTheFly(versions, meta_yml_paths, ['multiqcsav'])`.

**Returns:**

- `Map`: Citations keyed by tool name (`[citation: ..., bibliography: ...]`), restricted to tools that ran.

:::note{title="Tool-name matching"}
Tool names from the `versions` topic are matched against `meta.yml` tool keys **exactly first, then case-insensitively** (so `FastQC` lines up with `fastqc`). A tool that ran but has no matching `meta.yml` citation entry is logged as a warning — add it to that module's `meta.yml` `tools:` section to silence it.
:::

**Pipeline Integration (`utils_nfcore_<pipeline>_pipeline/main.nf`):**

This is the reproducible replacement for a hand-maintained `toolReferencesMap` and a static `toolCitationText()` / `toolBibliographyText()`. Wire it into the methods-description channel:

```nextflow title="subworkflows/local/utils_nfcore_<pipeline>_pipeline/main.nf"
include { citationsOnTheFly; methodsDescriptionText } from 'plugin/nf-core-utils'

// Resolve every module meta.yml once (citation metadata source).
// files() returns an empty list when nothing matches the glob.
def meta_yml_paths = files("${projectDir}/modules/**/meta.yml").collect { it.toString() }

ch_methods_description = channel.topic('versions')   // every module emits here when it runs
    .collect(flat: false)                            // flat: false keeps each [process, tool, version] tuple
    .map { topic_versions ->
        // Citations for ONLY the tools that ran, resolved from meta.yml
        def citations = citationsOnTheFly(topic_versions, meta_yml_paths)
        methodsDescriptionText(
            "${projectDir}/assets/methods_description_template.yml",
            citations,
            [:]
        )
    }

ch_multiqc_files = ch_multiqc_files.mix(
    ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true)
)
```

**Before / After:**

```nextflow title="Before — static, hand-maintained"
// Every tool listed by hand; cited even if it never ran
def citation_text = [
    "Tools used in the workflow included:",
    "BWAMEM2 (Vasimuddin et al. 2019)",
    "FastQC (Andrews 2010),",
    // ...one line per tool, forever...
].join(' ').trim()
```

```nextflow title="After — citations on the fly, zero maintenance"
// Cites exactly the tools that ran, metadata from each module's meta.yml
def citations = citationsOnTheFly(topic_versions, meta_yml_paths)
def citation_text = toolCitationText(citations)
```

---

#### `toolsFromVersionsTopic()`

**Description:**
Extracts the unique, sorted set of tool names that actually executed from collected `versions` topic data. Useful on its own — e.g. for conditional logic like `if ('star' in tools)`. (`citationsOnTheFly()` uses it internally to select which citations to emit.)

**Function Signature:**

```nextflow
List<String> toolsFromVersionsTopic(List topicVersions)
```

**Parameters:**

- `topicVersions` (List): Collected `versions` topic data — `[process, tool, version]` tuples (collect with `.collect(flat: false)`).

**Returns:**

- `List<String>`: Sorted, de-duplicated tool names. Empty list for `null`/empty input.

**Usage Example:**

```nextflow
include { toolsFromVersionsTopic } from 'plugin/nf-core-utils'

channel.topic('versions').collect(flat: false).map { versions ->
    def tools = toolsFromVersionsTopic(versions)
    log.info "Tools that ran: ${tools.join(', ')}"
    // e.g. "Tools that ran: fastqc, multiqc, samtools"
}
```

:::tip{title="meta.yml vs. the citation topic"}
`citationsOnTheFly()` keeps the citation source as canonical `meta.yml` files while gating by what ran — no need to add `getCitation(...)` emissions to every module. If you have already migrated modules to emit a dedicated `citation` topic, the [`getCitation` + `autoToolCitationText`](#getcitation) path below is equivalent and needs no `meta.yml` paths.
:::

---

### Modern Approach Functions (Recommended)

### `getCitation()`

**Description:**  
Extracts citation information from a module's meta.yml file and formats it for topic channel emission. This is the modern approach for automatic citation collection that only includes citations for tools that actually execute.

**Function Signature:**

```nextflow
List getCitation(String metaYmlPath)
```

**Parameters:**

- `metaYmlPath` (String): Path to the module's meta.yml file, typically `"${moduleDir}/meta.yml"`

**Returns:**

- `List`: Citation data formatted for topic channel emission `[module, tool, citation_data]`

**Usage in Process:**

```nextflow
include { getCitation } from 'plugin/nf-core-utils'

process FASTQC {
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip
    val citation_data, topic: citations

    script:
    citation_data = getCitation("${moduleDir}/meta.yml")
    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}
```

**Pipeline Integration:**

```nextflow
workflow {
    FASTQC(samples)
    MULTIQC(reports)

    // Collect all citations automatically - only from processes that executed
    citation_ch = channel.topic('citations').collect()

    // Generate citation text and bibliography
    citation_ch.map { citations ->
        def citationText = autoToolCitationText(citations)
        def bibliography = autoToolBibliographyText(citations)

        // Write to MultiQC template files
        file("${params.outdir}/multiqc_citations.txt").text = citationText
        file("${params.outdir}/multiqc_bibliography.html").text = bibliography
    }
}
```

---

### `autoToolCitationText()`

**Description:**  
Automatically generates formatted citation text from topic channel citation data. This processes the automatic citation collection from executed processes, ensuring runtime accuracy.

**Function Signature:**

```nextflow
String autoToolCitationText(List topicCitations = [])
```

**Parameters:**

- `topicCitations` (List): Citation data collected from topic channel

**Returns:**

- `String`: Formatted citation text for methods descriptions

**Usage Example:**

```nextflow
include { autoToolCitationText } from 'plugin/nf-core-utils'

workflow {
    // Process data and collect citations
    PROCESS_SAMPLES(samples)

    // Collect all citations from topic channel
    def citation_ch = channel.topic('citations').collect()

    citation_ch.map { citations ->
        def citationText = autoToolCitationText(citations)
        log.info "Generated citations: ${citationText}"
        // Returns: "Tools used in the workflow included: fastqc (Andrews et al. 2010), multiqc (Ewels et al. 2016)."
    }
}
```

---

### `autoToolBibliographyText()`

**Description:**  
Automatically generates HTML bibliography from topic channel citation data. Creates bibliography entries for tools that were actually executed.

**Function Signature:**

```nextflow
String autoToolBibliographyText(List topicCitations = [])
```

**Parameters:**

- `topicCitations` (List): Citation data collected from topic channel

**Returns:**

- `String`: HTML bibliography for MultiQC reports

**Usage Example:**

```nextflow
include { autoToolBibliographyText } from 'plugin/nf-core-utils'

workflow {
    // After processing with citation-enabled modules
    def citation_ch = channel.topic('citations').collect()

    citation_ch.map { citations ->
        def bibliography = autoToolBibliographyText(citations)
        // Write bibliography for MultiQC
        file("${params.outdir}/multiqc_bibliography.html").text = bibliography
    }
}
```

---

### Legacy Functions (Still Supported)

### `generateModuleToolCitation()`

**Description:**  
Extracts citation information from a module's meta.yml file using the traditional approach. This function processes all tools in the meta.yml regardless of execution.

**Function Signature:**

```nextflow
Map generateModuleToolCitation(Object metaFilePath)
```

**Parameters:**

- `metaFilePath` (String|File): Path to the meta.yml file

**Returns:**

- `Map`: Citations map with tool names as keys and citation/bibliography data as values

**Usage Example:**

```nextflow
include { generateModuleToolCitation } from 'plugin/nf-core-utils'

// Extract citations from specific module
def fastqcCitations = generateModuleToolCitation('modules/nf-core/fastqc/meta.yml')
// Returns: [fastqc: [citation: 'fastqc (Andrews et al. 2010)', bibliography: '<li>Andrews S...</li>']]

log.info "FASTQC citations: ${fastqcCitations}"
```

---

### `toolCitationText()`

**Description:**  
Generates formatted citation text from collected citations map (legacy approach).

**Function Signature:**

```nextflow
String toolCitationText(Map collectedCitations)
```

**Parameters:**

- `collectedCitations` (Map): Map of tool citations from `generateModuleToolCitation()`

**Returns:**

- `String`: Formatted citation text for use in methods descriptions

**Usage Example:**

```nextflow
include { generateModuleToolCitation; toolCitationText } from 'plugin/nf-core-utils'

// Collect citations from multiple modules
def allCitations = [:]
allCitations.putAll(generateModuleToolCitation('modules/nf-core/fastqc/meta.yml'))
allCitations.putAll(generateModuleToolCitation('modules/nf-core/multiqc/meta.yml'))

def citationText = toolCitationText(allCitations)
log.info "Citations: ${citationText}"
// Returns: "Tools used in the workflow included: fastqc (Andrews et al. 2010), multiqc (Ewels et al. 2016)."
```

---

### `toolBibliographyText()`

**Description:**  
Generates HTML bibliography from collected citations (legacy approach).

**Function Signature:**

```nextflow
String toolBibliographyText(Map collectedCitations)
```

**Parameters:**

- `collectedCitations` (Map): Map of tool citations

**Returns:**

- `String`: HTML bibliography for MultiQC reports

**Usage Example:**

```nextflow
include { collectCitationsFromFiles; toolBibliographyText } from 'plugin/nf-core-utils'

def citationFiles = ['modules/nf-core/fastqc/meta.yml', 'modules/nf-core/multiqc/meta.yml']
def citations = collectCitationsFromFiles(citationFiles)
def bibliography = toolBibliographyText(citations)

file("${params.outdir}/bibliography.html").text = bibliography
```

---

### `methodsDescriptionText()`

**Description:**  
Generates methods description text using collected citations and a MultiQC methods template. This function substitutes citation and bibliography variables into the template.

**Function Signature:**

```nextflow
String methodsDescriptionText(String mqcMethodsYamlPath, Map collectedCitations = [:], Map meta = [:])
```

**Parameters:**

- `mqcMethodsYamlPath` (String): Path to MultiQC methods YAML template file
- `collectedCitations` (Map, optional): Map containing all tool citations
- `meta` (Map, optional): Additional metadata for template substitution

**Returns:**

- `String`: Formatted methods description HTML

**Usage Example:**

```nextflow
include { collectCitationsFromFiles; methodsDescriptionText } from 'plugin/nf-core-utils'

// Collect citations and generate methods description
def citations = collectCitationsFromFiles(['modules/nf-core/fastqc/meta.yml'])
def methodsYaml = "${projectDir}/assets/methods_description_template.yml"

def methodsDescription = methodsDescriptionText(methodsYaml, citations, [
    'pipeline_name': workflow.manifest.name,
    'pipeline_version': workflow.manifest.version
])

file("${params.outdir}/methods_description.html").text = methodsDescription
```

**Template Example (`methods_description_template.yml`):**

```yaml
id: "custom-methods"
section_name: "nf-core/rnaseq Methods Description"
description: |
  <h4>Methods</h4>
  <p>Data was processed using nf-core/rnaseq v${pipeline_version}.
  ${tool_citations_text}</p>

  <h4>References</h4>
  ${tool_bibliography_text}
```

---

### `collectCitationsFromFiles()`

**Description:**  
Collects citations from multiple meta.yml files (legacy batch processing approach).

**Function Signature:**

```nextflow
Map collectCitationsFromFiles(List<String> metaFilePaths)
```

**Parameters:**

- `metaFilePaths` (List<String>): List of paths to module meta.yml files

**Returns:**

- `Map`: Combined citations from all files

**Usage Example:**

```nextflow
include { collectCitationsFromFiles } from 'plugin/nf-core-utils'

// Batch collect citations from multiple modules
def metaFiles = [
    'modules/nf-core/fastqc/meta.yml',
    'modules/nf-core/multiqc/meta.yml',
    'modules/local/custom_tool/meta.yml'
]

def allCitations = collectCitationsFromFiles(metaFiles)
log.info "Collected ${allCitations.size()} tool citations"
```

---

### Migration Functions

### `processMixedCitationSources()`

**Description:**  
Processes citations from both topic channels (modern approach) and legacy files for progressive migration.

**Function Signature:**

```nextflow
Map processMixedCitationSources(List<List> topicCitations, List<String> citationFiles)
```

**Parameters:**

- `topicCitations` (List<List>): Citation data from topic channels
- `citationFiles` (List<String>): List of meta.yml file paths

**Returns:**

- `Map`: Combined citations from both sources

**Migration Example:**

```nextflow
include { processMixedCitationSources } from 'plugin/nf-core-utils'

workflow {
    // Modern processes using topic channels
    MODERN_PROCESS(input_modern)

    // Legacy processes still using file-based citations
    LEGACY_PROCESS(input_legacy)

    // Collect from both sources
    def topicCitations = channel.topic('citations').collect()
    def legacyFiles = ['modules/legacy/old_tool/meta.yml']

    // Combine both approaches during migration
    def allCitations = processMixedCitationSources(topicCitations, legacyFiles)

    // Generate final citation text
    def citationText = toolCitationText(allCitations)
}
```

---

### `convertMetaYamlToTopicFormat()`

**Description:**  
Converts legacy meta.yml data to new topic channel format for migration purposes.

**Function Signature:**

```nextflow
List<List> convertMetaYamlToTopicFormat(String metaFilePath, String moduleName = null)
```

**Parameters:**

- `metaFilePath` (String): Path to meta.yml file
- `moduleName` (String, optional): Name of the module (defaults to filename)

**Returns:**

- `List<List>`: List of `[module, tool, citation_data]` tuples

**Migration Usage:**

```nextflow
include { convertMetaYamlToTopicFormat } from 'plugin/nf-core-utils'

// Convert existing meta.yml to topic format for consistency
def legacyMetaPath = 'modules/legacy/old_tool/meta.yml'
def topicFormatCitations = convertMetaYamlToTopicFormat(legacyMetaPath, 'OLD_TOOL')

// Emit to topic channel for unified processing
channel.fromList(topicFormatCitations)
    .set { converted_citations_ch }
```

---

### Topic Channel Processing Functions

### `processCitationsFromTopic()`

**Description:**  
Processes citations exclusively from topic channel format (modern approach).

**Function Signature:**

```nextflow
Map processCitationsFromTopic(List<List> topicData)
```

**Parameters:**

- `topicData` (List<List>): List of `[module, tool, citation_data]` tuples

**Returns:**

- `Map`: Processed citations map

**Usage Example:**

```nextflow
include { processCitationsFromTopic } from 'plugin/nf-core-utils'

workflow {
    // Pure topic channel approach
    PROCESS_A(input_a)
    PROCESS_B(input_b)

    // Collect topic channel citations
    def citations_ch = channel.topic('citations').collect()

    citations_ch.map { topicCitations ->
        def citations = processCitationsFromTopic(topicCitations)
        def citationText = toolCitationText(citations)

        file("${params.outdir}/citations.txt").text = citationText
    }
}
```

---

### `processCitationsFromFile()`

**Description:**  
Processes citations exclusively from traditional YAML files (legacy approach).

**Function Signature:**

```nextflow
Map processCitationsFromFile(List<String> citationFiles)
```

**Parameters:**

- `citationFiles` (List<String>): List of file paths to meta.yml files

**Returns:**

- `Map`: Tool citations from files

**Legacy Usage:**

```nextflow
include { processCitationsFromFile } from 'plugin/nf-core-utils'

// Legacy file-based processing
def citationFiles = [
    'modules/nf-core/fastqc/meta.yml',
    'modules/nf-core/multiqc/meta.yml'
]

def citations = processCitationsFromFile(citationFiles)
def citationText = toolCitationText(citations)
```

---

## Topic Channel Migration Guide

### Stage 1: Legacy File-Based Approach

Traditional approach collecting citations from all meta.yml files:

```nextflow
include { collectCitationsFromFiles; toolCitationText; toolBibliographyText } from 'plugin/nf-core-utils'

workflow {
    // Process samples with various tools
    FASTQC(samples)
    MULTIQC(reports)

    // Manually collect citations from all potential modules
    def metaFiles = [
        'modules/nf-core/fastqc/meta.yml',
        'modules/nf-core/multiqc/meta.yml'
    ]

    // Problem: Citations included even if process didn't execute
    def citations = collectCitationsFromFiles(metaFiles)
    def citationText = toolCitationText(citations)
    def bibliography = toolBibliographyText(citations)
}
```

### Stage 2: Mixed Migration Approach

Gradual migration supporting both approaches:

```nextflow
include { processMixedCitationSources; toolCitationText } from 'plugin/nf-core-utils'

workflow {
    // Modern processes using topic channels
    MODERN_FASTQC(samples)  // Emits to citations topic

    // Legacy processes still using manual collection
    LEGACY_MULTIQC(reports)  // No topic channel emission

    // Collect from both sources
    def topicCitations = channel.topic('citations').collect()
    def legacyFiles = ['modules/nf-core/multiqc/meta.yml']

    // Combine approaches during migration
    def allCitations = processMixedCitationSources(topicCitations, legacyFiles)
    def citationText = toolCitationText(allCitations)
}
```

### Stage 3: Modern Topic Channel Approach

Fully migrated with automatic runtime-accurate citations:

```nextflow
include { getCitation; autoToolCitationText; autoToolBibliographyText } from 'plugin/nf-core-utils'

process FASTQC {
    output:
    val citation_data, topic: citations

    script:
    citation_data = getCitation("${moduleDir}/meta.yml")
    """
    fastqc ${reads}
    """
}

workflow {
    // All processes automatically emit citations when they execute
    FASTQC(samples)
    MULTIQC(reports)

    // Automatic collection - only includes executed tools
    def citation_ch = channel.topic('citations').collect()

    citation_ch.map { citations ->
        def citationText = autoToolCitationText(citations)
        def bibliography = autoToolBibliographyText(citations)

        // Zero-maintenance: accurate citations without manual collection
        file("${params.outdir}/citations.txt").text = citationText
        file("${params.outdir}/bibliography.html").text = bibliography
    }
}
```

---

## Best Practices

### 1. Progressive Module Migration

Migrate modules gradually to topic channels:

```nextflow
// Phase 1: High-priority modules first
process CRITICAL_ANALYSIS {
    output:
    val citation_data, topic: citations  // Modernize first

    script:
    citation_data = getCitation("${moduleDir}/meta.yml")
    """
    # Critical analysis here
    """
}

// Phase 2: Supporting modules
process SUPPORTING_TOOL {
    // Still using legacy approach temporarily
    script:
    """
    # Supporting analysis
    """
}

// Phase 3: Utility modules last
process UTILITY_PROCESS {
    // Convert when convenient
}
```

### 2. Runtime Citation Accuracy

Ensure citations reflect actual tool execution:

```nextflow
workflow {
    // Conditional processing with accurate citations
    if (params.run_fastqc) {
        FASTQC(samples)  // Only emits citation if it runs
    }

    if (params.run_multiqc) {
        MULTIQC(reports)  // Only emits citation if it runs
    }

    // Citations automatically match actual execution
    channel.topic('citations').collect().map { citations ->
        autoToolCitationText(citations)  // Only includes tools that ran
    }
}
```

### 3. Error Handling and Graceful Degradation

Handle missing citation data gracefully:

```nextflow
include { getCitation; autoToolCitationText } from 'plugin/nf-core-utils'

process ROBUST_PROCESS {
    output:
    val citation_data, topic: citations

    script:
    try {
        citation_data = getCitation("${moduleDir}/meta.yml")
    } catch (Exception e) {
        log.warn "Citation extraction failed for ${task.process}: ${e.message}"
        citation_data = ['UNKNOWN', 'unknown_tool', [description: 'Citation unavailable']]
    }
    """
    # Process execution
    """
}
```

### 4. Testing Citation Functions

Test citation processing across migration stages:

```groovy
test("Modern topic channel citations") {
    when {
        def topicCitations = [
            ['FASTQC', 'fastqc', [doi: '10.1093/bioinformatics/btv033', author: 'Andrews S']],
            ['MULTIQC', 'multiqc', [doi: '10.1093/bioinformatics/btw354', author: 'Ewels P']]
        ]
        def result = autoToolCitationText(topicCitations)
    }

    then {
        assert result.contains('fastqc')
        assert result.contains('multiqc')
        assert result.contains('Andrews')
    }
}

test("Legacy file citations") {
    when {
        def citations = generateModuleToolCitation('test_meta.yml')
        def result = toolCitationText(citations)
    }

    then {
        assert result.contains('Tools used in the workflow')
        assert citations.size() > 0
    }
}

test("Mixed citation processing") {
    when {
        def topicCitations = [['MODERN', 'modern_tool', [doi: 'test']]]
        def legacyFiles = ['test_legacy_meta.yml']
        def result = processMixedCitationSources(topicCitations, legacyFiles)
    }

    then {
        assert result.containsKey('modern_tool')
        // Also contains legacy citations
    }
}
```

---

## Integration Examples

### Complete Modern Pipeline Integration

```nextflow
#!/usr/bin/env nextflow

include { getCitation; autoToolCitationText; autoToolBibliographyText } from 'plugin/nf-core-utils'

// Modern citation-enabled processes
process FASTQC {
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    val citation_data, topic: citations

    script:
    citation_data = getCitation("${moduleDir}/meta.yml")
    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}

process MULTIQC {
    input:
    path('*')

    output:
    path("multiqc_report.html"), emit: report
    val citation_data, topic: citations

    script:
    citation_data = getCitation("${moduleDir}/meta.yml")
    """
    multiqc .
    """
}

workflow {
    // Process samples
    FASTQC(samples)

    // Only run MultiQC if reports exist
    if (params.skip_multiqc == false) {
        MULTIQC(FASTQC.out.html.collect())
    }

    // Automatic citation collection - runtime accurate
    def citation_ch = channel.topic('citations').collect()

    citation_ch.map { citations ->
        // Generate citation components
        def citationText = autoToolCitationText(citations)
        def bibliography = autoToolBibliographyText(citations)

        // Write files for MultiQC template substitution
        file("${params.outdir}/pipeline_info/citations.txt").text = citationText
        file("${params.outdir}/pipeline_info/bibliography.html").text = bibliography

        log.info "Generated citations for ${citations.size()} tools"
    }
}
```

### MultiQC Methods Template Integration

```nextflow
include { autoToolCitationText; methodsDescriptionText } from 'plugin/nf-core-utils'

workflow {
    // Process with citation-enabled modules
    ANALYSIS_PROCESSES(data)

    // Collect citations and generate methods description
    def citation_ch = channel.topic('citations').collect()

    citation_ch.map { citations ->
        // Generate components
        def citationText = autoToolCitationText(citations)
        def bibliography = autoToolBibliographyText(citations)

        // Create comprehensive methods description
        def methodsYaml = "${projectDir}/assets/methods_description_template.yml"
        def methodsDescription = methodsDescriptionText(
            methodsYaml,
            [
                tool_citations_text: citationText,
                tool_bibliography_text: bibliography
            ],
            [
                pipeline_name: workflow.manifest.name,
                pipeline_version: workflow.manifest.version,
                pipeline_doi: workflow.manifest.doi ?: 'N/A'
            ]
        )

        // Write for MultiQC inclusion
        file("${params.outdir}/multiqc_methods_description.html").text = methodsDescription
    }
}
```

---

## Common Issues and Solutions

### Issue: Empty Citation Output

**Problem:**

```nextflow
def citations = autoToolCitationText([])
// Result: "Tools used in the workflow included: ."
```

**Solution:**

```nextflow
// Verify citations are being emitted
channel.topic('citations').view { "Citation collected: $it" }

// Check citation collection
channel.topic('citations').collect().view { "All citations: $it" }

// Handle empty case gracefully
def citations = channel.topic('citations').collect().map { citationList ->
    if (citationList.isEmpty()) {
        return "No tool citations available."
    } else {
        return autoToolCitationText(citationList)
    }
}
```

### Issue: Missing Module Citations

**Problem:**
Process doesn't emit citations even with `getCitation()` call

**Solution:**

```nextflow
process DEBUG_CITATIONS {
    output:
    val citation_data, topic: citations

    script:
    citation_data = getCitation("${moduleDir}/meta.yml")

    // Debug: Log citation data
    log.info "Citation data for ${task.process}: ${citation_data}"

    """
    # Check if meta.yml exists
    ls -la ${moduleDir}/
    """
}
```

### Issue: Legacy to Modern Migration Errors

**Problem:**
Mixed processing produces inconsistent citation formats

**Solution:**

```nextflow
// Validate citation data format before processing
def validatedCitations = topicCitations.findAll { citation ->
    citation instanceof List &&
    citation.size() == 3 &&
    citation[0] && citation[1] && citation[2]
}

def citations = processCitationsFromTopic(validatedCitations)
```

---

## Performance Considerations

### Topic Channel Benefits

- **Runtime Accuracy**: Only tools that execute appear in citations
- **Zero Maintenance**: No manual citation file management
- **Better Performance**: No file I/O for citation collection
- **Automatic Updates**: Citations update when modules change

### Migration Performance

```nextflow
// Optimize citation processing during migration
workflow {
    // Collect in parallel
    def topicFuture = channel.topic('citations').collect()
    def legacyFiles = legacy_citation_files.collect()

    // Process efficiently
    tuple(topicFuture, legacyFiles).map { topic, files ->
        processMixedCitationSources(topic, files)
    }.set { final_citations }
}
```

For more information on citation management and topic channels, see the [Nextflow documentation](https://www.nextflow.io/docs/latest/channel.html#topic) and [nf-core modules documentation](https://nf-co.re/docs/specifications/components/modules/documentation#documenting-of-all-tools).
