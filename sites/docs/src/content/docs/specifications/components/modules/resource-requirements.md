---
title: Resource requirements
subtitle: Specify resource requirements
markdownPlugin: addNumbersToHeadings
shortTitle: Resource requirements
weight: 6
---

The keywords "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Use of labels in modules

An appropriate set of resource `label` directives MUST be provided for the module as listed in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/conf/base.config).

The template defines two label schemes. New and updated modules SHOULD use the axis-decomposed labels; the bundled labels are retained for backwards compatibility and are scheduled for deprecation.

### Axis-decomposed labels (preferred)

Three independent axes cover CPU, memory and time. A module picks at most one label per axis and combines them to describe its resource shape:

| Axis   | Labels (low to high)                                                                  |
| ------ | ------------------------------------------------------------------------------------- |
| CPU    | `process_cpus_single`, `process_cpus_low`, `process_cpus_medium`, `process_cpus_high` |
| Memory | `process_mem_low`, `process_mem_medium`, `process_mem_high`                           |
| Time   | `process_time_short`, `process_time_medium`, `process_time_long`                      |

A CPU-bound but memory-light short job declares:

```groovy title="main.nf"
process FOO {
    label 'process_cpus_high'
    label 'process_mem_low'
    label 'process_time_short'
    ...
}
```

A module MAY omit a label on an axis it does not need to tune. Axes without an explicit label receive the per-process defaults defined at the top of `base.config`.

A module MUST NOT declare more than one label on the same axis. `process_cpus_low` and `process_cpus_high` together is invalid and is rejected by the linter.

### Bundled labels (legacy)

The original label set ties CPU, memory and time together. Use exactly one per module: `process_single`, `process_low`, `process_medium`, or `process_high`. Two modifier labels override a single axis and SHOULD be stacked on top of a bundled label: `process_long` (extends time), `process_low_memory` (drops memory), `process_high_memory` (raises memory).

These labels remain functional but emit a lint warning, and SHOULD be migrated to the axis-decomposed scheme. See the [migration guide](/docs/contributing/migrating-resource-labels) for the legacy-to-axis mapping table and a recipe for site-config maintainers.

### Label precedence

When a process matches more than one `withLabel:` block, Nextflow applies them in the order they appear in the config file. The order of the `label` directives in the process itself does not affect precedence.

The nf-core template defines bundled labels first, axis-decomposed labels next, and the legacy modifier labels (`process_long`, `process_low_memory`, `process_high_memory`) last. Stacking a modifier on top of a bundled label resolves to the modifier's value on the axis it overrides, and the bundled values on the others.

## Source of multiple threads or cores value

If the tool supports multi-threading, provide the appropriate parameter using the Nextflow `task` variable.
For example, `--threads $task.cpus`.

If the tool does not support multi-threading, consider `process_cpus_single` unless large amounts of RAM are required.

## Specifying multiple threads for piped commands

If a module contains _multiple_ tools that support multi-threading (e.g., [piping output into a samtools command](https://github.com/nf-core/modules/blob/c4cc1db284faba9fc4896f64bddf7703cedc7430/modules/nf-core/bowtie2/align/main.nf#L47-L54)), assign CPUs per tool:

- [`task.cpus`] is supplied unchanged when a process uses multiple cores.
- If one tool is multi-threaded and another uses a single thread, specify directly in the command itself. For example, with [`${task.cpus}`](https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/sampe/main.nf#L34).
