---
title: Minimum inputs
subtitle: Pipelines should be able to run with as little input as possible.
menu:
  main:
    weight: 150
---

Pipelines can accept as many input files as you like, but it should be possible to run with as few as possible.

For example, pipelines should auto-generate missing reference files, where possible.
So given a reference genome Fasta file a pipeline would build the reference index files.
The pipeline should also be able to optionally accept the reference index files in this case, if available.
