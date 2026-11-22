---
title: Avoid duplicating industry-standard workflows
subtitle: Pipelines should add value beyond established industry-standard workflows
shortTitle: Industry-standard workflows
menu:
  main:
    weight: 240
---

The keywords "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

An industry-standard workflow is an openly available, actively maintained workflow that is widely adopted for a given data or analysis type.
Examples include workflows maintained by instrument vendors or technology providers.

This recommendation applies to reimplementing a workflow.
nf-core pipelines SHOULD still use established tools and follow established analysis protocols.

## Pipelines

A pipeline SHOULD NOT reproduce an industry-standard workflow with the same scope, inputs, and outputs.

A pipeline that overlaps an industry-standard workflow MAY be accepted if it offers functionality that the industry-standard workflow does not.
You SHOULD describe that functionality in your pipeline proposal and agree it with the community before development.

## Pipeline steps

When an industry-standard workflow commonly produces the input for a step, the pipeline SHOULD accept the output of that workflow as input.
The pipeline MAY also run the step itself as an optional entry point.

:::info{title="Rationale" collapse}
Reproducing an industry-standard workflow splits users and development effort between two solutions to the same problem.
It also adds maintenance work to keep the nf-core implementation in step with upstream tools and recommendations.
Accepting standard outputs as input lets users start from the data they already have.
:::

## Example: Oxford Nanopore basecalling

Oxford Nanopore Technologies maintains the [EPI2ME workflows](https://github.com/epi2me-labs), including [wf-basecalling](https://github.com/epi2me-labs/wf-basecalling), which basecalls raw signal data.
This recommendation applies to basecalling as follows:

- A pipeline that only basecalls a directory of signal data for a single sample SHOULD NOT be developed, as it reproduces wf-basecalling.
- A pipeline that analyses Oxford Nanopore data SHOULD accept basecalled reads, such as FASTQ or BAM files, as input.
- A pipeline MAY optionally basecall from raw signal data before its analysis. For example, [nf-core/methylong](https://nf-co.re/methylong) accepts either modification-basecalled BAM files or raw POD5 files.
- A dedicated basecalling pipeline MAY be proposed if it adds functionality that wf-basecalling does not provide, such as processing multiple samples or sequencing runs from a samplesheet.

:::tip
See [Workflow specificity](../requirements/workflow_specificity) for the equivalent requirement on overlap between nf-core pipelines, and [Build with the community](build_with_community) for how to propose a pipeline before you build it.
:::
