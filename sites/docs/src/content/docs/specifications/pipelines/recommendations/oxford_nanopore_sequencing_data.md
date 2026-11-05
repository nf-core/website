---
title: Oxford Nanopore sequencing data
subtitle: Pipelines that analyse Oxford Nanopore sequencing data should take basecalled reads as inputs and avoid duplicating EPI2ME workflows
shortTitle: Oxford Nanopore sequencing data
menu:
  main:
    weight: 240
---

The keywords "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

Pipelines that analyse Oxford Nanopore data SHOULD take basecalled reads as input, and SHOULD NOT run basecalling themselves.

:::info{title="Accepted basecalled formats"}
- FASTQ (`.fastq.gz`)
- Unaligned or aligned BAM
- Unaligned or aligned modBAM (BAM with `MM`/`ML` modification tags), for modification-aware analyses
:::

### Dorado - Oxford Nanopore's multipurpose tool

ONT's recommended basecalling and alignment tool, [dorado](https://github.com/nanoporetech/dorado), is distributed under the
[Oxford Nanopore Technologies PLC Public Licence 1.0](https://github.com/nanoporetech/dorado/blob/master/LICENCE.txt),
which permits research use only. nf-core pipelines must be usable freely by everyone, so basecalling is
left to users as an upstream step.

## Overlap with EPI2ME workflows

Proposed pipelines SHOULD check for overlap with existing [EPI2ME workflows](https://github.com/epi2me-labs)
and discuss it with the community before development.
