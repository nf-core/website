---
# Copied from docs/ReferencesExtension.md of nf-core/nf-core-utils 0.6.0 by .github/workflows/add-nf-core-utils-docs.yml: edit it there
title: "References Extension Guide"
shortTitle: "References Extension Guide"
weight: 3
parentWeight: 3
md_github_url: https://github.com/nf-core/nf-core-utils/blob/main/docs/ReferencesExtension.md
---

# References Extension Guide

The References Extension simplifies genome reference and annotation file management in Nextflow pipelines. It provides intelligent parameter resolution that seamlessly handles both user-provided files and standardized reference collections like iGenomes, making your pipelines more flexible and user-friendly.

## 1. Overview

Managing reference files is a common challenge in bioinformatics pipelines. Users might provide:

- Custom reference files via parameters
- References from standardized collections (iGenomes)
- Mixed approaches depending on the analysis

The References Extension solves this by providing a unified interface that automatically resolves the appropriate reference source based on user input and pipeline configuration.

### 1.1. Core Capabilities

The extension provides two essential functions:

- **File Resolution**: `getReferencesFile()` - Resolves file paths from parameters or reference metadata
- **Value Resolution**: `getReferencesValue()` - Retrieves metadata values with parameter override support

### 1.2. Migration from Legacy Systems

This function is used to retrieve genome attributes in the nf-core TEMPLATE.

It retrieves a specific attribute (such as `fasta`, `gtf`, or index paths) for the selected genome from the `params.genomes` map. It is useful for pipelines that support multiple genomes and need to access reference files or metadata for the currently selected genome.

**Function Signature**

```groovy
Object getGenomeAttribute(String attribute)
```

or (static utility):

```groovy
Object ReferencesUtils.getGenomeAttribute(Map params, String attribute)
```

**Parameters**

| Parameter   | Type   | Required | Description                                                        |
| ----------- | ------ | -------- | ------------------------------------------------------------------ |
| `attribute` | String | Yes      | The attribute name to retrieve (e.g. 'fasta', 'gtf', 'star')       |
| `params`    | Map    | Yes      | (static) The Nextflow params map containing `genome` and `genomes` |

**Return Value**

Returns the value of the requested attribute for the selected genome, or `null` if not found.

**Practical Example**

```groovy
// Example params structure
def params = [
    genome: 'GRCh38',
    genomes: [
        GRCh38: [
            fasta: 's3://bucket/genome.fa',
            gtf: 's3://bucket/genes.gtf',
        ],
        GRCh37: [
            fasta: 's3://bucket/genome37.fa'
            star: 's3://bucket/star_index/'
        ]
    ]
]

// Retrieve the FASTA file for the selected genome
def fasta = ReferencesUtils.getGenomeAttribute(params, 'fasta')
// Returns: 's3://bucket/genome.fa'

// Retrieve the GTF file
def gtf = ReferencesUtils.getGenomeAttribute(params, 'gtf')
// Returns: 's3://bucket/genes.gtf'

// If the attribute or genome is missing, returns null
def missing = ReferencesUtils.getGenomeAttribute(params, 'star')
// Returns: null
```

**Usage in Nextflow**

```nextflow
include { getGenomeAttribute } from 'plugin/nf-core-utils'

workflow {
    // Example: get the FASTA file for the selected genome
    genome_fasta = getGenomeAttribute('fasta')
    log.info "Selected genome FASTA: ${genome_fasta}"
}
```

## 2. Getting Started

### 2.1. Basic Import and Setup

Let's start with a simple example that demonstrates the core concept:

```nextflow title="basic_references.nf"
#!/usr/bin/env nextflow

// Import reference utilities
include { getReferencesFile  } from 'plugin/nf-core-utils'
include { getReferencesValue } from 'plugin/nf-core-utils'

// Pipeline parameters
params.fasta = null          // User can override with custom file
params.genome = 'GRCh38'     // Default to standard reference
params.igenomes_base = 's3://ngi-igenomes/igenomes'

workflow {
    log.info "Setting up genome references for ${params.genome}"

    // Create references channel (example structure)
    references_ch = Channel.of([
        genome: params.genome,
        fasta: "${params.igenomes_base}/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
    ])

    // Get reference file - user parameter takes precedence
    genome_fasta = getReferencesFile(
        references_ch,           // Reference metadata channel
        params.fasta,           // User-provided parameter (null = use metadata)
        'fasta',               // Attribute to look for in metadata
        params.igenomes_base   // Base path for reference resolution
    )

    genome_fasta.view { "Using genome: ${it}" }
}
```

```console title="Output"
N E X T F L O W  ~  version 25.04.0
Launching `basic_references.nf` [peaceful-darwin] - revision: abc1234

INFO: Setting up genome references for GRCh38
Using genome: s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
```

### 2.2. Understanding Parameter Precedence

The extension follows a clear precedence hierarchy:

1. **User Parameters** (highest priority) - Direct file paths provided by users
2. **Reference Metadata** (fallback) - Files from reference collections
3. **Default Values** (lowest priority) - Pipeline defaults

```nextflow title="parameter_precedence.nf"
#!/usr/bin/env nextflow

params.fasta = "/custom/path/genome.fa"  // User override
params.genome = 'GRCh38'

workflow {
    references_ch = Channel.of([
        genome: 'GRCh38',
        fasta: 's3://igenomes/GRCh38/genome.fa'  // This will be ignored
    ])

    // User parameter takes precedence
    genome_fasta = getReferencesFile(references_ch, params.fasta, 'fasta', null)

    genome_fasta.view { "Selected: ${it}" }  // Shows custom path
}
```

## 3. Core Functions Reference

### 3.1. getReferencesFile - File Path Resolution

This function intelligently resolves file paths based on user parameters and reference metadata.

**Parameters**

| Parameter    | Type        | Required | Description                                   |
| ------------ | ----------- | -------- | --------------------------------------------- |
| `references` | Channel     | Yes      | Channel containing reference metadata         |
| `param`      | String/null | Yes      | User-provided file path (null = use metadata) |
| `attribute`  | String      | Yes      | Metadata attribute name to extract            |
| `basepath`   | String      | No       | Base path for relative path resolution        |

**Practical Example**

```nextflow title="file_resolution_example.nf"
#!/usr/bin/env nextflow

include { getReferencesFile } from 'plugin/nf-core-utils'

params.fasta = null
params.gtf = "/custom/annotations.gtf"
params.igenomes_base = 's3://ngi-igenomes/igenomes'

workflow {
    // Create comprehensive reference metadata
    references = Channel.of([
        genome: 'GRCh38',
        fasta: 'Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa',
        gtf: 'Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf',
        readme: 'Homo_sapiens/NCBI/GRCh38/README.txt'
    ])

    // Resolve multiple reference files
    genome_fasta = getReferencesFile(references, params.fasta, 'fasta', params.igenomes_base)
    genome_gtf = getReferencesFile(references, params.gtf, 'gtf', params.igenomes_base)

    // Combine for downstream processing
    references_ready = genome_fasta.combine(genome_gtf)

    references_ready.view { fasta, gtf ->
        """
        Reference files ready:
        FASTA: ${fasta}
        GTF: ${gtf}
        """
    }
}
```

**Practical Example**

```nextflow title="file_resolution_example.nf"
#!/usr/bin/env nextflow

include { getReferencesFile } from 'plugin/nf-core-utils'

params.fasta = null
params.gtf = "/custom/annotations.gtf"
params.igenomes_base = 's3://ngi-igenomes/igenomes'

workflow {
    // Create comprehensive reference metadata
    references = Channel.of([
        genome: 'GRCh38',
        fasta: 'Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa',
        gtf: 'Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf',
        readme: 'Homo_sapiens/NCBI/GRCh38/README.txt'
    ])

    // Resolve multiple reference files
    genome_fasta = getReferencesFile(references, params.fasta, 'fasta', params.igenomes_base)
    genome_gtf = getReferencesFile(references, params.gtf, 'gtf', params.igenomes_base)

    // Combine for downstream processing
    references_ready = genome_fasta.combine(genome_gtf)

    references_ready.view { fasta, gtf ->
        """
        Reference files ready:
        FASTA: ${fasta}
        GTF: ${gtf}
        """
    }
}
```

### 3.2. getReferencesValue - Metadata Value Resolution

This function extracts metadata values with user parameter override support.

**Function Signature**

```groovy
Object getReferencesValue(
    Channel references,    // Reference metadata channel
    Object param,         // User parameter value or null
    String attribute     // Metadata attribute name
)
```

**Parameters**

| Parameter    | Type     | Required | Description                               |
| ------------ | -------- | -------- | ----------------------------------------- |
| `references` | Channel  | Yes      | Channel containing reference metadata     |
| `param`      | Any/null | Yes      | User-provided value (null = use metadata) |
| `attribute`  | String   | Yes      | Metadata attribute name to extract        |

**Practical Example**

```nextflow title="value_resolution_example.nf"
#!/usr/bin/env nextflow

include { getReferencesValue } from 'plugin/nf-core-utils'

params.species = null        // Use metadata default
params.build = "custom_v2"   // Override metadata

workflow {
    references = Channel.of([
        genome: 'GRCh38',
        species: 'Homo sapiens',
        build: 'GRCh38.p13',
        assembly_date: '2020-12-01'
    ])

    // Extract various metadata values
    species_name = getReferencesValue(references, params.species, 'species')
    genome_build = getReferencesValue(references, params.build, 'build')
    assembly_date = getReferencesValue(references, null, 'assembly_date')

    // Combine all metadata
    metadata = species_name.combine(genome_build).combine(assembly_date)

    metadata.view { species, build, date ->
        """
        Genome Metadata:
        Species: ${species}        // From metadata (params.species = null)
        Build: ${build}           // From params.build override
        Date: ${date}            // From metadata
        """
    }
}
```

## 4. Integration Patterns

### 4.1. Complete Reference Pipeline

Here's a comprehensive example showing how to integrate reference resolution into a real bioinformatics pipeline:

```nextflow title="complete_reference_pipeline.nf" hl_lines="8-11 18-27 35-42"
#!/usr/bin/env nextflow

include { getReferencesFile  } from 'plugin/nf-core-utils'
include { getReferencesValue } from 'plugin/nf-core-utils'

// Pipeline parameters with sensible defaults
params.input = 'samples.csv'
params.fasta = null                    // User can override
params.gtf = null                     // User can override
params.genome = 'GRCh38'              // Default reference
params.igenomes_base = 's3://ngi-igenomes/igenomes'

workflow {
    log.info "Starting analysis with genome: ${params.genome}"

    // Load reference metadata (in practice, this might come from a YAML file)
    references = Channel.of([
        genome: params.genome,
        species: 'Homo sapiens',
        fasta: "Homo_sapiens/NCBI/${params.genome}/Sequence/WholeGenomeFasta/genome.fa",
        gtf: "Homo_sapiens/NCBI/${params.genome}/Annotation/Genes/genes.gtf",
        star_index: "Homo_sapiens/NCBI/${params.genome}/Sequence/STARIndex/",
        build: "${params.genome}.p13"
    ])

    // Resolve reference files intelligently
    genome_fasta = getReferencesFile(references, params.fasta, 'fasta', params.igenomes_base)
    genome_gtf = getReferencesFile(references, params.gtf, 'gtf', params.igenomes_base)
    star_index = getReferencesFile(references, null, 'star_index', params.igenomes_base)

    // Extract metadata values
    species_name = getReferencesValue(references, null, 'species')
    genome_build = getReferencesValue(references, null, 'build')

    // Use in downstream processes
    ALIGNMENT(Channel.fromPath(params.input), genome_fasta, star_index)
    ANNOTATION(ALIGNMENT.out.bam, genome_gtf)

    // Create analysis report with metadata
    species_name.combine(genome_build).view { species, build ->
        log.info "Analysis completed for ${species} (${build})"
    }
}

process ALIGNMENT {
    input:
    path samples
    path fasta
    path index

    output:
    path "*.bam", emit: bam

    script:
    """
    echo "Aligning samples using:"
    echo "Reference: ${fasta}"
    echo "Index: ${index}"
    # STAR alignment commands would go here
    touch aligned.bam
    """
}

process ANNOTATION {
    input:
    path bam
    path gtf

    output:
    path "*.counts", emit: counts

    script:
    """
    echo "Counting features using: ${gtf}"
    # featureCounts commands would go here
    touch counts.txt
    """
}
```

### 4.2. Working with Reference Collections

For pipelines using standardized reference collections, create a systematic approach:

```nextflow title="igenomes_integration.nf"
#!/usr/bin/env nextflow

include { getReferencesFile  } from 'plugin/nf-core-utils'
include { getReferencesValue } from 'plugin/nf-core-utils'

params.genome = 'GRCh38'
params.igenomes_base = 's3://ngi-igenomes/igenomes'

// User overrides (any can be null to use defaults)
params.fasta = null
params.gtf = null
params.bed12 = null

workflow {
    // Define comprehensive iGenomes structure
    igenomes_references = Channel.of([
        genome: params.genome,
        species: 'Homo sapiens',
        provider: 'NCBI',
        build: 'GRCh38.p13',
        fasta: "Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa",
        fasta_fai: "Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa.fai",
        gtf: "Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf",
        bed12: "Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.bed",
        star_index: "Homo_sapiens/NCBI/GRCh38/Sequence/STARIndex/",
        bowtie2_index: "Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2Index/"
    ])

    // Resolve all required references
    genome_fasta = getReferencesFile(igenomes_references, params.fasta, 'fasta', params.igenomes_base)
    genome_fasta_fai = getReferencesFile(igenomes_references, null, 'fasta_fai', params.igenomes_base)
    genome_gtf = getReferencesFile(igenomes_references, params.gtf, 'gtf', params.igenomes_base)
    genome_bed12 = getReferencesFile(igenomes_references, params.bed12, 'bed12', params.igenomes_base)

    // Create reference bundle for processes
    reference_bundle = genome_fasta
        .combine(genome_fasta_fai)
        .combine(genome_gtf)
        .combine(genome_bed12)

    reference_bundle.view { fasta, fai, gtf, bed ->
        """
        Reference Bundle Ready:
        - FASTA: ${fasta}
        - Index: ${fai}
        - GTF: ${gtf}
        - BED12: ${bed}
        """
    }
}
```

## 5. Advanced Usage Patterns

### 5.1. Multiple Genome Support

For pipelines supporting multiple genomes:

```nextflow title="multi_genome_support.nf"
#!/usr/bin/env nextflow

include { getReferencesFile  } from 'plugin/nf-core-utils'
include { getReferencesValue } from 'plugin/nf-core-utils'

// Support for multiple genomes
params.genomes = ['GRCh38', 'mm10']
params.igenomes_base = 's3://ngi-igenomes/igenomes'

workflow {
    // Create references for multiple genomes
    genome_configs = Channel.fromList([
        [genome: 'GRCh38', species: 'Homo sapiens', provider: 'NCBI'],
        [genome: 'mm10', species: 'Mus musculus', provider: 'UCSC']
    ])

    // Add file paths to each genome config
    references = genome_configs.map { config ->
        config + [
            fasta: "${config.species.replace(' ', '_')}/${config.provider}/${config.genome}/Sequence/WholeGenomeFasta/genome.fa",
            gtf: "${config.species.replace(' ', '_')}/${config.provider}/${config.genome}/Annotation/Genes/genes.gtf"
        ]
    }

    // Resolve references for each genome using proper channel operations
    resolved_references = references
        .map { ref ->
            tuple(ref,
                getReferencesFile(Channel.of(ref), null, 'fasta', params.igenomes_base),
                getReferencesFile(Channel.of(ref), null, 'gtf', params.igenomes_base)
            )
        }
        .flatMap { ref, fasta_ch, gtf_ch ->
            fasta_ch
                .combine(gtf_ch)
                .map { fasta, gtf ->
                    [
                        genome: ref.genome,
                        species: ref.species,
                        fasta: fasta,
                        gtf: gtf
                    ]
                }
        }

    resolved_references.view { ref ->
        "Ready: ${ref.genome} (${ref.species}) - FASTA: ${ref.fasta}, GTF: ${ref.gtf}"
    }
}
```

### 5.2. Custom Reference Validation

Add validation to ensure reference files exist and are compatible:

```nextflow title="reference_validation.nf"
#!/usr/bin/env nextflow

include { getReferencesFile } from 'plugin/nf-core-utils'

params.fasta = null
params.igenomes_base = 's3://ngi-igenomes/igenomes'

workflow {
    references = Channel.of([
        genome: 'GRCh38',
        fasta: 'Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa'
    ])

    genome_fasta = getReferencesFile(references, params.fasta, 'fasta', params.igenomes_base)

    // Validate reference file
    VALIDATE_REFERENCE(genome_fasta)
}

process VALIDATE_REFERENCE {
    input:
    path fasta

    output:
    path fasta, emit: validated_fasta
    stdout emit: validation_report

    script:
    """
    # Check if file exists and is not empty
    if [[ ! -s "${fasta}" ]]; then
        echo "ERROR: Reference file ${fasta} is empty or doesn't exist"
        exit 1
    fi

    # Check FASTA format
    if ! grep -q "^>" "${fasta}"; then
        echo "ERROR: ${fasta} doesn't appear to be a valid FASTA file"
        exit 1
    fi

    # Count sequences
    num_sequences=\$(grep -c "^>" "${fasta}")
    echo "Validated FASTA with \${num_sequences} sequences"
    """
}
```

## 6. Best Practices

### 6.1. Reference Organization

:::tip{title="Parameter Naming"}
Use consistent parameter names across your pipeline:

- `params.fasta` for genome sequences
- `params.gtf` for gene annotations
- `params.bed12` for BED format annotations
- `params.{tool}_index` for tool-specific indices

:::

### 6.2. Error Handling

```nextflow title="robust_error_handling.nf"
#!/usr/bin/env nextflow

include { getReferencesFile } from 'plugin/nf-core-utils'

params.fasta = null
params.genome = 'GRCh38'

workflow {
    // Validate required parameters
    if (!params.genome && !params.fasta) {
        error "Either --genome or --fasta must be provided!"
    }

    references = Channel.of([
        genome: params.genome,
        fasta: params.genome ? "genomes/${params.genome}/genome.fa" : null
    ])

    // Handle missing reference gracefully
    try {
        genome_fasta = getReferencesFile(references, params.fasta, 'fasta', null)
        genome_fasta.view { "Using reference: ${it}" }
    } catch (Exception e) {
        log.error "Failed to resolve reference: ${e.message}"
        log.error "Please check --genome parameter or provide --fasta directly"
        System.exit(1)
    }
}
```

### 6.3. Documentation

Always document your reference requirements:

```nextflow title="documented_references.nf"
#!/usr/bin/env nextflow

/*
 * REFERENCE FILES
 *
 * This pipeline supports flexible reference file specification:
 *
 * Option 1: Use standard genome (automatic file resolution)
 *   --genome GRCh38
 *
 * Option 2: Provide custom files (override defaults)
 *   --fasta /path/to/genome.fa
 *   --gtf /path/to/annotations.gtf
 *
 * Option 3: Mix standard and custom
 *   --genome GRCh38 --gtf /custom/annotations.gtf
 *
 * Supported genomes: GRCh38, GRCh37, mm10, mm9
 */

include { getReferencesFile  } from 'plugin/nf-core-utils'
include { getReferencesValue } from 'plugin/nf-core-utils'

// Reference parameters with documentation
params.genome = null          // Standard genome name (e.g., 'GRCh38')
params.fasta = null          // Custom genome FASTA file
params.gtf = null           // Custom gene annotation file
params.igenomes_base = 's3://ngi-igenomes/igenomes'

workflow {
    // Implementation here...
}
```

## 7. Takeaway

The References Extension provides a powerful, flexible system for managing genome references in Nextflow pipelines:

1. **Unified Interface**: Single functions handle both custom files and reference collections
2. **Smart Resolution**: Automatic parameter precedence with user override support
3. **iGenomes Integration**: Seamless integration with standardized reference collections
4. **Pipeline Flexibility**: Users can mix custom and standard references as needed

## 8. What's Next?

- Explore **[NfCore Utilities](NfCoreUtilities.md)** for comprehensive pipeline utilities
- Learn about **[NextflowPipelineExtension](NextflowPipelineExtension.md)** for core pipeline functions
- Check out the utility documentation in `utilities/` for specialized functions
