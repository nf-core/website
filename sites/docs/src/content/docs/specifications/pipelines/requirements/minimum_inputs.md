---
title: Minimum inputs
subtitle: Pipelines should be able to run with as little input as possible.
menu:
  main:
    weight: 150
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

Pipelines MAY accept as many input files as needed, but it SHOULD be possible to run with as few as possible.

Pipelines SHOULD auto-generate missing reference files where possible.
For example, given a reference genome Fasta file a pipeline would build the reference index files.
The pipeline SHOULD also be able to optionally accept the reference index files if available.
