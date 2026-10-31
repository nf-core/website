---
# Copied from docs/ReferencesExtension.md of nf-core/nf-core-utils 0.7.0 by .github/workflows/add-nf-core-utils-docs.yml: edit it there
title: "References Extension Guide"
shortTitle: "References Extension Guide"
weight: 3
parentWeight: 3
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

The extension provides:

- `getGenomeAttribute()`: reads one attribute of the selected genome from `params.genomes` (see [1.2](#12-migration-from-legacy-systems)).
- `getReferencePath()`, `getReferenceMetadata()` and `getReferenceAttribute()`: resolve references from a references datasheet (see [2.1](#21-reference-datasheet-resolution)).

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

## 2. Core Functions Reference

### 2.1. Reference Datasheet Resolution

A [references datasheet](https://github.com/nf-core/references-datasheets) describes one genome. Its keys are pipeline reference parameter names.

```groovy
String getReferencePath(Map args)
// args: references, genome, igenomesIgnore, datasheetPath

Map getReferenceMetadata(Map args)
// args: rows (required)

Object getReferenceAttribute(Map args)
// args: param (required), datasheet, genomeMap, genomeKey
```

Missing arguments default to `null`. An unknown argument name throws an `IllegalArgumentException`.

`getReferencePath()` returns the datasheet to load:

| Input                                 | Result                         |
| ------------------------------------- | ------------------------------ |
| `references` set                      | `references`                   |
| `genome` set and `igenomesIgnore` set | `<datasheetPath>/<genome>.yml` |
| otherwise                             | `null`                         |

It throws if `references` and `genome` are both set, or if it needs `datasheetPath` and none is given. A `genome` of `null`, `''` or `'null'` counts as unset, so `--genome null` turns off a default genome.

`getReferenceMetadata()` takes the output of `samplesheetToList()` and returns the genome's non-empty values as a map. It throws unless the datasheet holds exactly one genome. Map every property of the datasheet schema to `meta`; other properties are dropped.

`getReferenceAttribute()` returns `datasheet[param]`, else `genomeMap[genomeKey][param]`, else `null`.

```nextflow
include { getReferencePath      } from 'plugin/nf-core-utils'
include { getReferenceMetadata  } from 'plugin/nf-core-utils'
include { getReferenceAttribute } from 'plugin/nf-core-utils'
include { samplesheetToList     } from 'plugin/nf-schema'

workflow {
    reference_path = getReferencePath(
        references: params.references,
        genome: params.genome,
        igenomesIgnore: params.igenomes_ignore,
        datasheetPath: params.references_datasheets_path,
    )
    datasheet = reference_path
        ? getReferenceMetadata(rows: samplesheetToList(reference_path, "${projectDir}/assets/schema_references.json"))
        : null
    fasta = getReferenceAttribute(param: 'fasta', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
}
```

Resolve references once, during pipeline initialisation, and pass the values to workflows as inputs. A `params.x = ...` assignment inside a workflow does not reach included files or config closures.

## 3. Best Practices

### 3.1. Reference Organization

:::tip{title="Parameter Naming"}
Use consistent parameter names across your pipeline:

- `params.fasta` for genome sequences
- `params.gtf` for gene annotations
- `params.bed12` for BED format annotations
- `params.{tool}_index` for tool-specific indices

:::

## 4. Takeaway

The References Extension provides a powerful, flexible system for managing genome references in Nextflow pipelines:

1. **Unified Interface**: Single functions handle both custom files and reference collections
2. **Smart Resolution**: Automatic parameter precedence with user override support
3. **iGenomes Integration**: Seamless integration with standardized reference collections
4. **Pipeline Flexibility**: Users can mix custom and standard references as needed

## 5. What's Next?

- Explore **[NfCore Utilities](NfCoreUtilities.md)** for comprehensive pipeline utilities
- Learn about **[NextflowPipelineExtension](NextflowPipelineExtension.md)** for core pipeline functions
- Check out the utility documentation in `utilities/` for specialized functions
