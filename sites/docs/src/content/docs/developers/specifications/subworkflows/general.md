---
title: General
subtitle: General specifications for nf-core Nextflow DSL2 subworkflows
weight: 1
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Minimum subworkflow size

Subworkflows SHOULD combine tools that make up a logical unit in an analysis step.
A subworkflow MUST contain at least two modules.

## Version reporting channel

Each subworkflow emits a channel that contains all tool versions from `versions.yml` files.
Collect the versions within the workflow and add them to the output as `versions`:

```bash
take:
  input

main:

  ch_versions = Channel.empty()

  FASTQC(input)

  ch_versions = ch_versions.mix(FASTQC.out.versions())

emit:
  versions = ch_versions
```
