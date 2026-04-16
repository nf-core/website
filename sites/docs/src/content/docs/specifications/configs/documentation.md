---
title: Documentation
subtitle: Follow config documentation guidelines
markdownPlugin: addNumbersToHeadings
shortTitle: Documentation
weight: 1
---

## General

A config MUST have a Markdown based documentation file describing the usage of the config.

- **`docs/<config name>.md`** for an institutional config
- **`docs/pipeline/<pipeline>/<config name>`** for a pipeline specific config

A pipeline specific config MUST be linked to within the institutional level config documentation.

## Scope of documentation

A config documentation Markdown file MUST contain information on how to use the config in the context of the infrastructure.
For example, the documentation SHOULD describe how to load Nextflow and other requirements for execution, such as `module load <...>` or `conda activate <...>`.

The documentation MAY also provide additional recommendations such as how to specify scratch or temporary directories using Nextflow-related [environment variables](https://www.nextflow.io/docs/latest/reference/env-vars.html#nextflow-settings).

The config documentation Markdown file MUST NOT contain generic information, such as installing Nextflow.

## Description of environmental variables

Any mandatory parameter or environment variable required for pipeline execution MUST be described in the documentation.
Non-mandatory but common parameters or environment variables SHOULD be described.

## Examples in documentation

It is RECOMMENDED to include copy-paste-able command examples of launching a nf-core workflow with the `test` profile (e.g. using `$USER` rather than placeholders).
