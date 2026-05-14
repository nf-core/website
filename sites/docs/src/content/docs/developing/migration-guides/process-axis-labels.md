---
title: Migrating to axis-decomposed resource labels
subtitle: Move pipelines, modules and site configs from bundled to axis-decomposed labels
description: How to update modules, pipelines and site configs from the bundled resource labels (process_low etc.) to the axis-decomposed labels (process_cpus_*, process_mem_*, process_time_*).
shortTitle: Migrating resource labels
---

The bundled resource labels (`process_low`, `process_medium`, `process_high`, `process_long`, `process_high_memory`) are being replaced by axis-decomposed labels (`process_cpus_*`, `process_mem_*`, `process_time_*`) that tune CPU, memory and time independently. Bundled labels still work and emit a lint warning; they are removed after a two-release deprecation window.

The full scheme is defined in [Resource requirements](/docs/specifications/components/modules/resource-requirements); the rationale is in [nf-core/proposals#139](https://github.com/nf-core/proposals/issues/139).

## Why migrate

- Most tools scale on only one axis; bundled labels over-provision the others, often by hundreds of times.
- Common shapes like "many CPUs, little memory" or "moderate CPUs, high memory" cannot be expressed cleanly with bundled labels.
- Pipelines have been working around this with locally invented labels (`process_high_cpu` etc.), which diverge across the ecosystem.

## Default mapping

Start from this table. The `process_low` row is opinionated: its 12 GB legacy default is closer to `process_mem_medium` than `process_mem_low`, so module authors should sanity-check the choice rather than apply it mechanically.

| Legacy label          | CPU                   | Memory               | Time                  |
| --------------------- | --------------------- | -------------------- | --------------------- |
| `process_single`      | `process_cpus_single` | (template default)   | (template default)    |
| `process_low`         | `process_cpus_low`    | `process_mem_medium` | (template default)    |
| `process_medium`      | `process_cpus_medium` | `process_mem_medium` | `process_time_medium` |
| `process_high`        | `process_cpus_high`   | `process_mem_high`   | `process_time_long`   |
| `process_long`        | (unchanged)           | (unchanged)          | `process_time_long`   |
| `process_high_memory` | (unchanged)           | `process_mem_high`   | (unchanged)           |
| `process_low_memory`  | (unchanged)           | `process_mem_low`    | (unchanged)           |

A `process_high` module becomes:

```diff groovy title="main.nf"
 process FOO {
-    label 'process_high'
+    label 'process_cpus_high'
+    label 'process_mem_high'
+    label 'process_time_long'
     ...
 }
```

A `process_medium` + `process_high_memory` stack becomes:

```diff groovy title="main.nf"
 process FOO {
-    label 'process_medium'
-    label 'process_high_memory'
+    label 'process_cpus_medium'
+    label 'process_mem_high'
+    label 'process_time_medium'
     ...
 }
```

Omit any axis label that matches the template default; that is the canonical way to say "this axis does not need to be tuned".

## Modules and pipelines

1. Replace each legacy label on the process with the axis labels from the table, adjusting per the tool's real scaling profile.
2. Apply the same change to any `conf/base.config` overrides in the pipeline, and drop pipeline-local label definitions that just duplicate the template.
3. Run `nf-core modules lint` to confirm the new labels are recognised.

The bundled-label `withLabel:` blocks in the template's `base.config` are removed by the nf-core template sync; they should not be edited by hand.

## Institutional site configs

Sites that override `withLabel: process_*` need parallel overrides for the new axis labels. The default translation is the mapping table above; institutions with custom tiering may want different values per axis.

A `sed` recipe to seed the CPU axis from existing overrides (always review the diff before committing):

```bash
sed -E '
  s/withLabel: ?process_low\b/withLabel:process_cpus_low/;
  s/withLabel: ?process_medium\b/withLabel:process_cpus_medium/;
  s/withLabel: ?process_high\b/withLabel:process_cpus_high/;
' conf/site.config
```

The recipe rewrites the CPU axis only. Memory and time overrides do not translate one-to-one and need review by hand.

## Deprecation timeline

| Stage    | Status                                                                                                |
| -------- | ----------------------------------------------------------------------------------------------------- |
| Now      | Both schemes supported. Lint warns on bundled labels.                                                 |
| +1 minor | Both schemes supported. nf-core/modules and nf-core/configs migrated.                                 |
| +2 minor | Bundled labels removed from the template. Pipelines that still use them must define them locally.     |

## See also

- [Resource requirements](/docs/specifications/components/modules/resource-requirements)
- RFC: [nf-core/proposals#139](https://github.com/nf-core/proposals/issues/139)
- Template change: [nf-core/tools#4265](https://github.com/nf-core/tools/pull/4265)
