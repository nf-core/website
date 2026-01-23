---
title: Documentation
subtitle: Requirements for documenting nf-core modules
markdownPlugin: addNumbersToHeadings
shortTitle: Documentation
weight: 4
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Module documentation is required

The module MUST have a `meta.yaml` in the same directory as the module `main.nf`.

## Number of keywords

Keywords SHOULD be sufficient to make the module findable through research domain, data types, and tool function keywords.
Keywords MUST NOT be solely the (sub)tool name.

:::info
For multi-tool modules, add the keyword `multi-tool` and all the (sub)tools.
:::

## Keyword formatting

Keywords MUST be all lower case.

## Documenting of all tools

The tools section MUST list every tool used in the module.
For example:

```yml
tools:
  - bowtie2: <....>
  - samtools: <....>
```

## Documentation of args of each piped or multiple command

The tools section MUST have a `args_id:` field for every tool in the module that describes which `$args` (`$args2`, `$args3`) variable is used for.
A single tool module will only have `args_id: "$args"`.

```yml
tools:
  - bowtie2:
      <...>
      args_id: "$args"
  - samtools:
      <...>
      args_id: "$args2"
```

## Required channel documentation

Only include entries of input and output channels in the Input and Output sections of the `meta.yaml`.

## Documentation of tuples

Split input and output tuples into separate entries.
That is, `meta` should be a separate entry to the `file` it is associated with.

## Input and output channel types

Use only the following categories for input/output types:

- `map`
- `file`
- `directory`
- `string`
- `boolean`
- `integer`
- `float`
- `boolean`
- `list`

## Correspondence of input/outputs entries to channels

Your input/output entries MUST match a corresponding channel in the module itself:

- There should be a one-to-one relationship between the module's inputs and outputs and those described in `meta.yml`.
- You MUST NOT combine multiple output channels in input/output entries.

## Useful input/output descriptions

Your input/output descriptions SHOULD be descriptive of the contents of file.
For example, not only 'A TSV file' but describe what the file contains.

## Input/output glob pattern

Your input/output patterns (if present) MUST follow a [Java glob pattern](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob).

## Ontology

Your module `meta.yml` files SHOULD contain a [bio.tools](https://bio.tools/) ID when available.

Your module `meta.yml` files SHOULD contain ontology URLs for files when relevant.

Some bioinformatics tools listed on `bio.tools` already have a list of inputs and outputs with their format that you can use (e.g., [FastQC](https://bio.tools/fastqc)), the EDAM ontology term can be obtained by selecting the relevant format in the input section of the diagrams.
Otherwise, you can get the ontology terms for a given format by searching the term in the EBI's [Ontology Lookup Service](https://www.ebi.ac.uk/ols4/ontologies/edam) (recommended), or the [EDAM browser](https://edamontology.github.io/edam-browser/#topic_0091).

## Indication of input channel requirement

You should mark input entries as Mandatory or Optional
