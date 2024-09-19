---
title: Overview
subtitle: nf-core pipeline configuration
shortTitle: Overview
weight: 1
parentWeight: 20
---

## How to configure nf-core pipelines

Each nf-core pipeline comes with a set of “sensible defaults” for a "typical" analysis of a full size dataset.
While the defaults are a great place to start, you will certainly want to modify these to fit your own data and system requirements. For example, modifying a tool flag of compute resources allocated for a process.

When a pipeline is launched, Nextflow will look for config files in several locations.
As each source can contain conflicting settings, the sources are ranked to decide which settings to apply.

Configuration sources are reported below and listed in order of priority:

1. Parameters specified on the command line (`--parameter`)
2. Parameters that are provided using the `-params-file` option
3. Config file that are provided using the `-c` option
4. The config file named `nextflow.config` in the current directory
5. The config file named `nextflow.config` in the pipeline project directory
6. The config file `$HOME/.nextflow/config`
7. Values defined within the pipeline script itself (e.g., `main.nf`)

While some of these files are already included in the nf-core pipeline repository (e.g., the `nextflow.config` file in the nf-core pipeline repository), some are automatically identified on your local system (e.g., the `nextflow.config` in the launch directory), and others are only included if they are specified using run options (e.g., `-params-file`, and `-c`).

:::warning
You should not clone and manually edit an nf-core pipeline. Manually edited nf-core pipelines cannot be updated to more recent versions of the pipeline without overwriting your changes. You also risk moving away from the canonical pipeline and losing reproducibility.
:::
