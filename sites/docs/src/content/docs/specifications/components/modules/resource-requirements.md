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

When labels are stacked, later labels override earlier ones for any axis they set. This lets a module express a shape that none of the bundled labels covers on its own. For example, a CPU-bound but memory-light tool:

process FOO {
    label 'process_high'
    label 'process_low_memory'
    ...
}
```

resolves to the CPU and time defaults of `process_high` with the memory of `process_low_memory`.

## Source of multiple threads or cores value

If the tool supports multi-threading, provide the appropriate parameter using the Nextflow `task` variable.
For example, `--threads $task.cpus`.

If the tool does not support multi-threading, consider `process_single` unless large amounts of RAM are required.

## Specifying multiple threads for piped commands

If a module contains _multiple_ tools that support multi-threading (e.g., [piping output into a samtools command](https://github.com/nf-core/modules/blob/c4cc1db284faba9fc4896f64bddf7703cedc7430/modules/nf-core/bowtie2/align/main.nf#L47-L54)), assign CPUs per tool:

- [`task.cpus`] is supplied unchanged when a process uses multiple cores.
- If one tool is multi-threaded and another uses a single thread, specify directly in the command itself. For example, with [`${task.cpus}`](https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/sampe/main.nf#L34).
