---
title: Cross-organisational components
subtitle: Use modules from multiple organisations in subworkflows
shortTitle: Cross-organisational components
weight: 3
---

You can specify individual remotes for each module within a subworkflow's `meta.yml` file.
This allows you to combine components from your own organisation with modules from other sources like nf-core.

This guide explains how to use modules from different organisations within a single subworkflow using `nf-core/tools`.

## Implementation

To mix modules from different sources, configure the `meta.yml` file with specific remotes for each component.

For example, the `fastq_trim_fastp_fastqc` subworkflow from [nf-core-test](https://github.com/nf-core-test/modules/tree/main/subworkflows/nf-core-test/fastq_trim_fastp_fastqc) organization combines:

- FastQC from the local organisation
- fastp from the nf-core repository

The `meta.yml` configuration uses this structure:

```yaml title="meta.yml" {11-12}
name: "fastq_trim_fastp_fastqc"
description: Read QC, fastp trimming and read qc
keywords:
  - qc
  - quality_control
  - adapters
  - trimming
  - fastq
components:
  - fastqc
  - fastp:
      git_remote: https://github.com/nf-core/modules.git
```

This tells `nf-core/tools` to install FastQC from the same repository as the subworkflow, while fastp comes from the `nf-core/modules` remote.

## Installation

When you install subworkflows from non-nf-core organisations, include the `--git-remote` parameter in your command.

## Schema configuration

For cross-organisational subworkflows, update your JSON schema to allow both object and string types in the components section.
This enables the flexible syntax shown above.

```json title="yaml-schema.json" {5}
"components": {
  "type": "array",
  "description": "Modules and subworkflows used in the subworkflow",
  "items": {
    "type": ["object", "string"]
  },
  "minItems": 0
},
```

:::tip
You can find a complete schema example in the [nf-core-test repository](https://github.com/nf-core/nf-core-test).
:::
