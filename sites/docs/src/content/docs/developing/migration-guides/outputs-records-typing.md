---
title: Migrating to workflow outputs, records and static typing
subtitle: Plan a migration to the output block, records and typed Nextflow code
description: A draft guide to adopting Nextflow workflow outputs, records and static typing together in nf-core pipelines
shortTitle: Outputs, records and typing
---

:::warning{title="Draft guidance"}
These patterns are still being proven in [nf-core/rnaseq](https://github.com/nf-core/rnaseq/pull/1945) (work in progress, tracked in [#1931](https://github.com/nf-core/rnaseq/issues/1931) and [#1935](https://github.com/nf-core/rnaseq/issues/1935)) and may change.
Some behaviour depends on Nextflow features merged but not yet released ([#7646](https://github.com/nextflow-io/nextflow/pull/7646), [#7674](https://github.com/nextflow-io/nextflow/pull/7674), [#7692](https://github.com/nextflow-io/nextflow/pull/7692)).
:::

Three newer Nextflow features change how a pipeline is written:

1. **Workflow outputs:** a top-level `output {}` block that replaces `publishDir`.
2. **Records:** `record(...)` values with named fields, instead of tuples.
3. **Static typing:** typed processes and workflows, enabled with `nextflow.enable.types`.

## Why adopt them together

Each feature pulls you towards the other two.

- The `output {}` block works best with one record per sample per stage, because its path closures read named fields (`r.quant_dir`) instead of tuple positions. Outputs pull you towards records.
- An untyped record is an unchecked bag of fields, and `join(by: 'id')` on records needs an `id` field. Typed processes want a single record output, and typed workflows drop `.out`. Records pull you towards typing.
- Typing removes or discourages `.out`, `branch`, `multiMap`, `groupTuple`, `collectFile` and tuple plumbing, which records replace. Typing pulls you towards records, and towards the output block as the place they end up.

Doing one alone leaves glue such as typed processes that emit tuples, or an output block fed by hand-rolled joins.
Plan the migration as one program of work, delivered in small steps.

## What a migrated pipeline looks like

- Every per-sample value is a record with an `id` (the join key). Keep an open `meta` map too where config closures read `meta.x` (see [meta maps](../components/meta-map)).
- Stages share field names (`reads`, `bam`, `bai`), so one stage's result record is the next stage's input with no `.map`.
- Each process takes one record and emits one record. Reference files are plain `Path` inputs.
- Workflows assign calls to variables instead of reading `.out`. `take:` and `emit:` carry named record types.
- One `output {}` block routes record fields to published locations.
- A module's own types live in its file. Shared and pipeline-level types live in one shared file.

```nextflow
nextflow.enable.types = true

process FASTQC {
    input:
    record(id: String, meta: Map, reads: List<Path>)

    output:
    record(id: id, meta: meta, html: files('*.html', optional: true), zip: files('*.zip', optional: true))

    script:
    """
    touch ${id}.html ${id}.zip
    """
}

workflow QC_QUANT {
    take:
    ch_samples: Channel<Sample>
    index:      Value<Path>

    main:
    ch_fastqc  = FASTQC(ch_samples)
    ch_quant   = SALMON_QUANT(ch_samples, index)
    ch_results = ch_fastqc.join(ch_quant, by: 'id')

    emit:
    ch_results
}

output {
    results: Channel<QuantResult> {
        path { r -> [([r.quant_dir]): "quant/${r.id}/"] }
    }
}
```

:::note{title="Needs Nextflow 26.10"}
The map form of `path` (a closure returning files mapped to targets) needs [nextflow#7646](https://github.com/nextflow-io/nextflow/pull/7646). Without it, use one `>>` per target.
:::

## Migration path

A suggested order, from experience so far.
All three steps have been applied across one full pipeline, but its full test suite had not been shown green at the time of writing, so runtime claims not marked verified are unproven.

1. **Result records first.** Give each stage a result record and feed the `output {}` block from them, keeping the published layout identical to the old `publishDir` layout.
2. **Type the processes.** Convert them one at a time to emit one record each.
3. **Records in, typed workflows.** Move inputs to records and type the workflows. Replace `.out`, `set`, `tap`, `|` and tuple plumbing with assigned calls and core operators (`join`, `map`, `filter`, `flatMap`, `mix`, `collect`, `groupBy`).
4. **Replace `collectFile`** with small modules or record-based collection. For samplesheets, the `index` directive works (field order is column order, but it quotes every CSV field).
5. **Verify** against the pre-migration pipeline (see below), then update tests and snapshots.

## Conventions that keep it simple

Several of these are proposals from a single proof of concept.

**Records and joins**

- Join with `join(by: 'id')`, never `mix` plus `groupTuple`. An unsized `groupTuple` waits for the channel to close and delays every downstream release to the end of the run. The same applies to joining an empty channel with `remainder: true`, so join optional stages conditionally.
- On a shared field the right-hand record wins, also for `left + record(...)`. Join a minimal partial record (`record(id: r.id, field: r.field)`) to add one stage's result.
- The typed `join` has no `failOnMismatch`. Where an unmatched sample can only be a bug, use `remainder: true`, a `subscribe` that calls `error`, then `.filter` the complete rows.
- Build each stage record right after its process call, then join records to records. A record-to-record join with `remainder: true` leaves unmatched fields absent instead of `null`, which changes `index` JSON.
- Mark optional files `Path?` and filter on `r.field != null` where an empty channel used to do the job.
- Field names are stable API: they appear in output docs and `index` headers. Reuse a module's own emit name.

**Processes**

- Use one unnamed record output. `file(glob)` needs `optional: true` if it may not match, `files()` is unordered so sort it, and a typed `List<Path>` interpolates as `[a, b]` so join it explicitly.
- A named record input (`sample: BamBaiInput`) binds one variable, so scripts, `task.ext` closures and templates read `sample.meta.id`. A destructured `record(id: String, ...)` keeps bare names. Use one input variable name across modules so shared config selectors still match.
- The output cast (`record(...) as Type`) makes lint check missing fields and wrong types, but not nullability or extra fields. Check your Nextflow includes the fix for [nextflow#7680](https://github.com/nextflow-io/nextflow/issues/7680) before casting remote files.
- Nextflow validates record inputs at runtime ([nextflow#7692](https://github.com/nextflow-io/nextflow/pull/7692)): a null field not marked `?` aborts the run.
- Emit versions on the topic channel (see [migrating to topic channels](update-pipelines)).

**Workflows and outputs**

- To `emit:` one result, assign the variable without `def`, or emit an expression.
- A `take:` typed `Value<...>` should receive a real value channel, including in tests (`channel.value(...)`).
- Everything published lives in a record: per-sample files in stage records, run-level files in small records grouped by output directory.
- Every target needs one live `>>`. If all resolve to `null` for a record, the target throws an NPE and the run can still exit 0 ([nextflow#7669](https://github.com/nextflow-io/nextflow/issues/7669), closed, so check your version). Use an anchor field, an `enabled` gate written as a top-level helper function, or a `.filter()` before `publish:`.
- Reference fields hold task outputs only. Files that did not come from a task need a static `path` ([nextflow#7667](https://github.com/nextflow-io/nextflow/issues/7667)), so set user-supplied pass-through files to `null`.
- Do not list a file twice in one closure (later entries win). A directory target keeps task-relative subfolders. Same-named files from different processes need separate outputs ([nextflow#6617](https://github.com/nextflow-io/nextflow/issues/6617)).
- Put shared prefix logic in top-level helper functions below the output block.
- Sort file lists passed to MultiQC so report inputs are deterministic.
- Keep `nextflow.enable.types = true` in each script that uses typed processes or workflows.

## Checking a migration

- Run `nextflow lint`, then the real tests early. Lint does not check that closures, templates and config expressions resolve at runtime.
- Compare with the pre-migration pipeline: same published file tree, same md5sums. Explain every difference by the migration, and update snapshots only afterwards.
- Keep edits to vendored nf-core modules minimal so the diff against upstream stays reviewable.
- After scripted renames, grep for leftovers. Configs and templates are not linted.

## Practical notes

- Plugin functions such as nf-schema's `samplesheetToList` run in typed scripts, but `nextflow lint` rejects them ([nextflow#7720](https://github.com/nextflow-io/nextflow/issues/7720)). A small untyped helper file is a pragmatic home.
- nf-core tooling infers dependencies from `include` paths and does not yet recognise type imports. Expect a tools update before types can be shared across installed components.
- Vendored components edited in place show as diverged in `nf-core` lint, and `nf-core modules update` must not be run until the changes are upstreamed (related: [nf-core/tools#3157](https://github.com/nf-core/tools/issues/3157)).
- `task.ext` and `meta` interact with record input naming. An upstream experiment replaces them with plain process inputs ([nf-core/methylseq#625](https://github.com/nf-core/methylseq/pull/625)).
- Declare numeric field types to match what tools produce, and leave values unconverted if their text appears in an output file.
- Workflow-level type checks are still catching up ([nextflow#7693](https://github.com/nextflow-io/nextflow/issues/7693)), so rely on runtime validation and tests. Publishing directories to object storage ([nextflow#7398](https://github.com/nextflow-io/nextflow/issues/7398)) is not yet verified.

Plan for these user-visible changes: `process.withName:<NAME>.publishDir` overrides in user configs stop working (breaking), directories are only created when a file is published into them, and the minimum Nextflow version rises.

## Open questions

- How should nf-core handle `task.ext.when` alongside typed `when:`?
- How do record-typed modules land in nf-core/modules, and where do shared types live?
- Should record types be cast at all, given [nextflow#7680](https://github.com/nextflow-io/nextflow/issues/7680)?
- How should new files in vendored components (such as `types.nf`) be tracked without patch files?
- Should `task.ext` and `meta` give way to plain process inputs?
- Which conventions above should become specifications, and which stay guidance?
