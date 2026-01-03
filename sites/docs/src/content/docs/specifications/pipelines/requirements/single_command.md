---
title: Single command
subtitle: Pipelines should run in a single command.
menu:
  main:
    weight: 110
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

Every nf-core pipeline repository MUST contain a single pipeline.
That is, there MUST be a `main.nf` file that is the single way to launch a pipeline.

- It is acceptable to have multiple 'tracks' within the pipeline, selectable with configuration options.
- It is acceptable to have workflows that use the output of another nf-core pipeline as input.

It MUST be possible to run all parts of the workflow using `nextflow run nf-core/<pipeline>`, without any specific filename.
