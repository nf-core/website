---
title: Input/output options
subtitle: Input and output specifications for nf-core Nextflow DSL2 subworkflows
weight: 3
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Required input channels

Input channel declarations MUST be defined for all _possible_ input files that the subworkflow requires (i.e., both required and optional files) within the `take` block.

## Required output channels

Named file extensions MUST be emitted for ALL output channels (e.g., `path "*.txt", emit: txt`).

## Optional inputs

Nextflow does not currently support optional inputs.
Pass an empty list (`[]`) instead of a file as a subworkflow parameter to work around this limitation.
