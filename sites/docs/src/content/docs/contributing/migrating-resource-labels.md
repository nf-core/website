---
title: Migrating to axis-decomposed resource labels
subtitle: Move pipelines, modules and site configs from bundled to axis-decomposed labels
shortTitle: Migrating resource labels
---

The nf-core pipeline template historically shipped a small set of bundled resource labels (`process_low`, `process_medium`, `process_high`, `process_long`, `process_high_memory`) that set CPU, memory and time together. These are being replaced by axis-decomposed labels (`process_cpus_*`, `process_mem_*`, `process_time_*`) that tune each axis independently.

The bundled labels remain functional for now. Linting emits a warning when a module uses them, and they will be removed after a two-release deprecation window.

This guide covers the migration steps for module authors, pipeline maintainers, and institutional site-config maintainers. The canonical definition of the label scheme lives in [Resource requirements](/docs/specifications/components/modules/resource-requirements).

## Why migrate

- Many tools scale on only one axis (CPU-bound, memory-bound, or wall-clock-bound). Bundled labels force authors to pick a tier that gives appropriate values on one axis and wildly over-provisions the others.
- Common shapes like "many CPUs, little memory" cannot be expressed cleanly with bundled labels.
- Pipelines have already invented ad-hoc additions (`process_high_cpu`, `process_high_memory`) to work around this, with no community-wide convention.

## Default mapping

Use this table as the starting point when migrating a module. The values are the closest equivalents in the current template; per-module judgment may override them, especially for `process_low` whose 12 GB memory default is closer to `process_mem_medium` than `process_mem_low`.

| Legacy label          | CPU                   | Memory               | Time                  |
| --------------------- | --------------------- | -------------------- | --------------------- |
| `process_single`      | `process_cpus_single` | (template default)   | (template default)    |
| `process_low`         | `process_cpus_low`    | `process_mem_medium` | (template default)    |
| `process_medium`      | `process_cpus_medium` | `process_mem_medium` | `process_time_medium` |
| `process_high`        | `process_cpus_high`   | `process_mem_high`   | `process_time_long`   |
| `process_long`        | (unchanged)           | (unchanged)          | `process_time_long`   |
| `process_high_memory` | (unchanged)           | `process_mem_high`   | (unchanged)           |
| `process_low_memory`  | (unchanged)           | `process_mem_low`    | (unchanged)           |

A module that previously declared `label 'process_high'` becomes:

```groovy title="main.nf"
process FOO {
    label 'process_cpus_high'
    label 'process_mem_high'
    label 'process_time_long'
    ...
}
```

A module that previously stacked `label 'process_medium'` with `label 'process_high_memory'` becomes:

```groovy title="main.nf"
process FOO {
    label 'process_cpus_medium'
    label 'process_mem_high'
    label 'process_time_medium'
    ...
}
```

## For module authors

1. Review the tool's actual scaling behaviour. Bundled labels often over-provision on at least two axes; this is the chance to right-size the request.
2. Replace each legacy label on the process with the equivalent axis labels from the mapping table, adjusting per the tool's profile.
3. Drop any axis label that matches the template default - omitting a label is the canonical way to say "this axis does not need to be tuned".
4. Run `nf-core modules lint` to confirm no warnings remain.

## For pipeline maintainers

1. Update local modules in `modules/local/` using the same approach as above.
2. Override values in the pipeline's `conf/base.config` only where the tool genuinely needs different resources from the template defaults.
3. Drop pipeline-local label definitions that simply duplicate the template (the template's labels are inherited from the nf-core template sync).
4. Keep the bundled-label `withLabel:` blocks in `conf/base.config` until upstream tools-template sync removes them.

## For institutional config maintainers

Sites that ship overrides for `withLabel: process_*` in their config bundle need parallel overrides for the new axis labels. The mapping table above gives the default translation; institutions with custom tiering may want different values per axis.

A minimal `sed` recipe to seed the new labels from existing overrides (always review the diff before committing):

```bash
# duplicate each bundled-label override into axis overrides
sed -E '
  s/withLabel: ?process_low\b/withLabel:process_cpus_low/;
  s/withLabel: ?process_medium\b/withLabel:process_cpus_medium/;
  s/withLabel: ?process_high\b/withLabel:process_cpus_high/;
' conf/site.config
```

The recipe deliberately rewrites to CPU axes only; memory and time should be reviewed by hand because the bundled labels conflate values in ways that do not translate one-to-one.

## Deprecation timeline

| Stage    | Status                                                                                                |
| -------- | ----------------------------------------------------------------------------------------------------- |
| Now      | Axis labels and bundled labels both supported. Lint emits a warning on bundled labels.                |
| +1 minor | Bundled labels still supported. nf-core/modules migrated. nf-core/configs migrated.                   |
| +2 minor | Bundled labels removed from the template. Pipelines that still reference them must define them locally or migrate. |

## See also

- [Resource requirements](/docs/specifications/components/modules/resource-requirements)
- RFC: [nf-core/proposals#139](https://github.com/nf-core/proposals/issues/139)
- Template change: [nf-core/tools#4265](https://github.com/nf-core/tools/pull/4265)
