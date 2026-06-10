---
title: Resource requirements
subtitle: Specify resource requirements
markdownPlugin: addNumbersToHeadings
shortTitle: Resource requirements
weight: 6
---

The keywords "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Use of labels in modules

An appropriate resource `label` MUST be provided for the module as listed in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/conf/base.config).

The template (since nf-core/tools:4.1.0) defines two kinds of labels:

- **Bundled labels** set CPU, memory and time together. Use exactly one per module: `process_single`, `process_low`, `process_medium`, or `process_high`.
- **Modifier labels** override a single axis and SHOULD be stacked on top of a bundled label: `process_long` (extends time), `process_low_memory` (drops memory), `process_high_memory` (raises memory).

When a process matches more than one `withLabel:` block, Nextflow applies them in the order they appear in the config file. The order of the `label` directives in the process itself does not affect precedence.
The nf-core template defines bundled labels first and modifier labels after, so modifier labels always overwrite the corresponding bundled value.

## Source of multiple threads or cores value

If the tool supports multi-threading, provide the appropriate parameter using the Nextflow `task` variable.
For example, `--threads $task.cpus`.

If the tool does not support multi-threading, consider `process_single` unless large amounts of RAM are required.

## GPU acceleration

Modules that support GPU acceleration SHOULD use `task.accelerator{:groovy}` to detect whether a GPU has been requested.
Pipelines control GPU allocation by setting `accelerator = 1{:groovy}` in their process config (e.g., via a `process_gpu` label or a `withName` block).

The module SHOULD NOT set the `accelerator` directive itself.

:::info{title="Rationale" collapse}
Placing GPU allocation in the pipeline config lets users control it through their pipeline config or profiles.
A label-only alternative (e.g., requiring a `process_gpu` label) would not work for modules that support both CPU and GPU modes (e.g., [`ribodetector`](https://github.com/nf-core/modules/tree/master/modules/nf-core/ribodetector)), so the specification leaves this to the pipeline author.
:::

See [Software requirements: GPU-capable modules](/docs/specifications/components/modules/software-requirements#gpu-capable-modules) for container patterns based on `task.accelerator`.

:::tip{title="Pipeline-side GPU configuration"}
Pipelines set `accelerator = 1{:groovy}` and GPU container flags via `containerOptions` in their process config.
Use `containerOptions` (not global `docker.runOptions`) to scope GPU flags to GPU processes only.
:::

:::caution{title="GPU concurrency under Singularity"}
Multiple concurrent GPU processes sharing a single GPU can deadlock under Singularity.
Docker's NVIDIA runtime handles GPU memory arbitration, but Singularity does not.
When GPU tasks may land on the same machine (CI, local executor, shared HPC nodes), set `maxForks = 1` on GPU processes to serialise access.
:::

## Specifying multiple threads for piped commands

If a module contains _multiple_ tools that support multi-threading (e.g., [piping output into a samtools command](https://github.com/nf-core/modules/blob/c4cc1db284faba9fc4896f64bddf7703cedc7430/modules/nf-core/bowtie2/align/main.nf#L47-L54)), assign CPUs per tool:

- [`task.cpus`] is supplied unchanged when a process uses multiple cores.
- If one tool is multi-threaded and another uses a single thread, specify directly in the command itself. For example, with [`${task.cpus}`](https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/sampe/main.nf#L34).
