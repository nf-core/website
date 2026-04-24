---
title: Documentation
subtitle: Follow config documentation guidelines
markdownPlugin: addNumbersToHeadings
shortTitle: Documentation
weight: 1
---

The keywords "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## General

A config MUST have a Markdown based documentation file describing the usage of the config.

- **`docs/<config name>.md`** for an institutional config
- **`docs/pipeline/<pipeline>/<config name>`** for an institutional pipeline specific config
- **`docs/pipeline/<pipeline>/<pipeline>.md`** for a global pipeline specific config

A pipeline specific config MUST be linked to within the institutional level config documentation.

## Scope of documentation

A config documentation Markdown file MUST contain information on how to use the config in the context of the infrastructure.
For example, the documentation SHOULD describe how to load Nextflow and other requirements for execution, such as `module load <...>` or `conda activate <...>`.

The documentation MAY also provide additional recommendations such as how to specify scratch or temporary directories using Nextflow-related [environment variables](https://www.nextflow.io/docs/latest/reference/env-vars.html#nextflow-settings).

The config documentation Markdown file MUST NOT contain generic information, such as installing Nextflow.

For offline or sensitive data infrastructures, the documentation MUST describe
any special steps required (e.g., pre-pulling containers, local reference paths).

## Description of environmental variables

Any mandatory parameter or environment variable required for pipeline execution MUST be described in the documentation.
Non-mandatory but common parameters or environment variables SHOULD be described.

## Examples in documentation

It is RECOMMENDED to include copy-paste-able command examples of launching a nf-core workflow with the `test` profile (e.g. using `$USER` rather than placeholders).
