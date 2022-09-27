---
title: Single command
subtitle: Pipelines should run in a single command.
menu:
  main:
    weight: 110
---

Every nf-core pipeline repository should contain a single pipeline.
That is, there should be a `main.nf` file that is the single way to launch a pipeline.

- It is ok to have multiple 'tracks' within the pipeline, selectable with configuration options.
- It is ok to have workflows that use the output of _another_ nf-core pipeline as input

It should be possible to run all parts of the workflow using `nextflow run nf-core/<pipeline>`, without any specific filename.
