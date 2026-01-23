---
title: Documentation
subtitle: Documentation specifications for nf-core Nextflow DSL2 subworkflows
weight: 5
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Code comment of channel structure

Each input and output channel SHOULD have a comment that describes the output structure of the channel.
For example:

```nextflow
input:
ch_reads // channel: [mandatory] meta, reads
val_sort // boolean: [mandatory] false
<...>

emit:
bam = SAMTOOLS_VIEW.out.bam // channel: [ val(meta), path(bam) ]
versions = ch_versions      // channel: [ path(versions.yml) ]
```

## Meta.yml documentation of channel structure

Each input and output channel structure SHOULD also be described in the `meta.yml` file's description entry.

```text
description: |
  Structure: [ val(meta), path(tsv)]
  (Sub)contig coverage table
```
