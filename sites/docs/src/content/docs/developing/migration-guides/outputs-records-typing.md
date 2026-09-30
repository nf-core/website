---
title: Migrating to workflow outputs, records and static typing
subtitle: Plan a migration to the output block, records and typed Nextflow code
description: A draft guide to adopting Nextflow workflow outputs, records and static typing together in nf-core pipelines
shortTitle: Outputs, records and typing
---

:::warning{title="Draft guidance"}
These patterns are still being proven in [nf-core/rnaseq](https://github.com/nf-core/rnaseq/pull/1945) (work in progress, tracked in [#1931](https://github.com/nf-core/rnaseq/issues/1931) and [#1935](https://github.com/nf-core/rnaseq/issues/1935)) and may change.
Some behaviour depends on Nextflow features merged but not yet released ([#7646](https://github.com/nextflow-io/nextflow/pull/7646), [#7674](https://github.com/nextflow-io/nextflow/pull/7674), [#7692](https://github.com/nextflow-io/nextflow/pull/7692)).
Snippets are trimmed from that work in progress.
:::

Three newer Nextflow features change how a pipeline is written:

1. **Workflow outputs:** a top-level `output {}` block replaces `publishDir`.
2. **Records:** `record(...)` values with named fields replace tuples.
3. **Static typing:** typed processes and workflows, enabled with `nextflow.enable.types`.

## Why adopt them together

Each feature pulls you towards the other two. Here is a small pipeline before and after (the "before" is an illustrative sketch of typical current nf-core style):

```nextflow
// Before: tuples, publishDir and .out
process SALMON_QUANT {
    publishDir "${params.outdir}/quant", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path('quant'), emit: quant
}

workflow QC_QUANT {
    take: ch_reads, index
    main:
    FASTQC(ch_reads)
    SALMON_QUANT(ch_reads, index)
    ch_results = FASTQC.out.zip.join(SALMON_QUANT.out.quant)   // positional: [meta, zip, quant]
    emit: results = ch_results
}
```

```nextflow
// After: records, types and an output block
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

**Outputs pull you towards records.**
Publishing needs to know which file is which.
A tuple such as `[meta, zip, quant]` is read by position, and every added file shifts the positions.
A record is read by name (`r.quant_dir`), so the output block stays readable as stages grow.

**Records pull you towards typing.**
An untyped record is an unchecked bag of fields.
A typo like `r.quant_dri` is `null` at runtime, whereas a declared `QuantResult` type lets the language server report it while you edit.
Typed processes also want a single record output, and `join(by: 'id')` needs an `id` field.

**Typing pulls you towards records and the output block.**
A typed workflow has no `.out`, so `FASTQC.out.zip` becomes `ch_fastqc = FASTQC(...)`, and a call has to return one thing.
Operators such as `groupTuple`, `branch` and `multiMap` are discouraged because the type checker cannot validate them.
Records with `join` and `map` replace them, and the output block is where those records end up.

Adopt only one and the glue stays.
Typed processes that emit tuples still need positional `.map` calls.
An output block fed by untyped tuples still breaks when a stage adds a file.
Together, the `publishDir` lines, `.out` access, positional joins and renaming `.map` calls in the "before" disappear.
Plan the migration as one program of work, delivered in small steps.

## What a migrated pipeline looks like

A typed process takes one record and emits one record. Its input and result types live in the module file:

```nextflow
nextflow.enable.types = true

record SamtoolsSortResult {
    id:   String
    meta: Map
    bam:  Path?
    bai:  Path?
}

process SAMTOOLS_SORT {
    tag "${sample.meta.id}"

    input:
    sample: RawBams
    fasta: Path?
    index_format: String

    output:
    record(
        id:   sample.id,
        meta: sample.meta,
        bam:  file("${prefix}.bam", optional: true),
        bai:  file("${prefix}.{bam,cram,sam}.bai", optional: true)
    ) as SamtoolsSortResult

    script:
    prefix = task.ext.prefix ?: "${sample.meta.id}"
    // ...
}
```

A typed workflow assigns calls to variables (no `.out`) and joins stage records on `id`:

```nextflow
workflow BAM_SORT_STATS_SAMTOOLS {
    take:
    ch_bam: Channel<RawBams>
    ch_fasta: Value<Path?>

    main:
    ch_sorted  = SAMTOOLS_SORT(ch_bam, ch_fasta, '')
    ch_results = ch_sorted.join(SAMTOOLS_STATS(ch_sorted), by: 'id')

    emit:
    ch_results
}
```

Design rules:

- Every per-sample value is a record with an `id`. Keep an open `meta` map where config closures read `meta.x` (see [meta maps](../components/meta-map)).
- Stages share field names (`reads`, `bam`, `bai`), so one stage's result is the next stage's input with no `.map`.
- Shared and pipeline-level types live in one file. A module's own types live in its `main.nf`.

## Migration path

All three steps have been applied across one full pipeline, but its full test suite had not been shown green at the time of writing, so treat runtime claims as unproven unless marked verified.

1. **Result records first.** Feed the `output {}` block from per-stage records, keeping the published layout identical to the old `publishDir` layout.
2. **Type the processes,** one at a time, each emitting one record.
3. **Records in, typed workflows.** Replace `.out`, `set`, `tap`, `|` and tuple plumbing with assigned calls and core operators (`join`, `map`, `filter`, `flatMap`, `mix`, `collect`, `groupBy`).
4. **Replace `collectFile`** with small modules, record-based collection, or the `index` directive.
5. **Verify** against the pre-migration pipeline, then update tests and snapshots.

## Conventions

Several of these are proposals from a single proof of concept.

### Joins

Join stage records with `join(by: 'id')`. Do not `mix` then `groupTuple`: an unsized `groupTuple` waits for the channel to close, which delays every downstream release to the end of the run.

The typed `join` has no `failOnMismatch`. Where an unmatched sample can only be a bug:

```nextflow
ch_sorted_indexed = ch_sorted.join(ch_index, by: 'id', remainder: true)
ch_sorted_indexed.subscribe { r ->
    if( r.bam == null || r.bai == null ) {
        error "Sample '${r.id}' is missing its samtools index result"
    }
}
ch_indexed = ch_sorted_indexed.filter { r -> r.bam != null && r.bai != null }
```

On a shared field the right-hand record wins. To add one stage's result to a running record, join a minimal partial record:

```nextflow
.map { r -> record(id: r.id, meta: r.meta, bam: r.bam, bai: r.bai, samtools: r.samtools) }
```

Build each stage record right after its process call, then join records to records. A record-to-record join with `remainder: true` leaves unmatched fields absent instead of `null`, which changes `index` JSON.

### Optional values

Mark optional files `Path?`. Filter where an empty channel used to do the job:

```nextflow
ch_sorted = SAMTOOLS_SORT(ch_bam, ch_fasta, ch_fai, '').filter { r -> r.bam != null }
```

Nextflow validates record inputs at runtime ([nextflow#7692](https://github.com/nextflow-io/nextflow/pull/7692)), so a null in a field not marked `?` aborts the run.

### Process outputs

- `file(glob)` needs `optional: true` if it may not match.
- `files()` is unordered, so sort it. A typed `List<Path>` interpolates as `[a, b]`, so join it explicitly.
- The output cast (`record(...) as Type`) makes lint check missing fields and wrong types, but not nullability or extra fields. Check that your Nextflow includes the fix for [nextflow#7680](https://github.com/nextflow-io/nextflow/issues/7680) before casting remote files.
- A named record input (`sample: BamBaiInput`) means scripts, `task.ext` closures and templates read `sample.meta.id`. A destructured `record(id: String, ...)` keeps bare names. Use one variable name across modules so shared config selectors still match.
- Emit versions on the topic channel (see [migrating to topic channels](update-pipelines)).

### Where record types live

A module declares its own input and result types, and imports shared ones:

```nextflow
// modules/nf-core/samtools/sort/main.nf
include { RawBams } from '../../types'

record SamtoolsSortResult {
    id:   String
    meta: Map
    bam:  Path?
    bai:  Path?
}
```

Shared vocabulary types live in one file that modules, subworkflows and the pipeline import from:

```nextflow
// modules/nf-core/types.nf
record Sample {
    id:   String
    meta: Map
}

record ReadsInput {
    id:    String
    meta:  Map
    reads: List<Path>
}
```

Workflows import types by path and use them to annotate channels:

```nextflow
// subworkflows/nf-core/bam_sort_stats_samtools/main.nf
include { RawBams; Bam } from '../../../modules/nf-core/types'
include { SamtoolsSortResult } from '../../../modules/nf-core/samtools/sort/main'

def ch_sorted: Channel<SamtoolsSortResult> = SAMTOOLS_SORT(ch_bam, ch_fasta, ch_fai, '')
```

- Declare a type used by several components once, in the component the others already depend on, and `include` it. Do not copy it.
- Open the shared file with a comment saying what belongs there and what stays in module files.
- An alternative layout gives each subworkflow a sibling `types.nf`, importing from its neighbours (`include { SamtoolsStatsFiles } from '../bam_stats_samtools/types'`). Which layout suits nf-core is an open question.
- Records that are never cast (see [nextflow#7680](https://github.com/nextflow-io/nextflow/issues/7680)) still document the expected shape.

### Parameters and output directory

Declare parameter types in a typed `params {}` block. A parameter that `nextflow.config` or `conf/` reads takes its default from the config, because the config is resolved before the block, so the block only declares its type. Point `outputDir` at the existing parameter:

```nextflow
// main.nf
params {
    input:   Path? = null
    outdir:  String?
    trimmer: String = 'trimgalore'
}

// nextflow.config
outputDir = params.outdir
```

### Output block

Locals hold shared prefixes, and a list key groups files that share a destination. Top-level helpers keep prefixes consistent:

```nextflow
def alignedDir(id: String) -> String { "${samplePrefix(id)}${params.aligner}/" }

output {
    stringtie: Channel<StringtieSample> {
        enabled !params.skip_stringtie
        path { s ->
            def dir = "${alignedDir(s.id)}stringtie/"
            [
                ([s.transcript_gtf, s.abundance, s.coverage_gtf]): dir,
                (s.ballgown):                                      dir,
            ]
        }
    }
}
```

:::note{title="Needs Nextflow 26.10"}
The map form of `path` (a closure returning files mapped to targets) needs [nextflow#7646](https://github.com/nextflow-io/nextflow/pull/7646). Without it, use one `>>` per target.
:::

An `index` target replaces `collectFile` for samplesheets. Field order is the column order, and every CSV field is quoted:

```nextflow
samplesheet: Channel<SamplesheetRow> {
    enabled params.save_align_intermeds && !params.skip_alignment
    index {
        path 'samplesheets/samplesheet_with_bams.csv'
        header true
    }
}
```

Rules that avoid silent failures:

- **Every target needs one live `>>`.** If all resolve to `null` for a record, the target throws an NPE and the run can still exit 0 ([nextflow#7669](https://github.com/nextflow-io/nextflow/issues/7669), closed, so check your version). Use an anchor field, an `enabled` gate, or a `.filter()` before `publish:`. Write a gate as a top-level helper so the gate and closure cannot drift apart.
- **Publish task outputs only.** Set user-supplied pass-through files to `null`, since files that did not come from a task need a static `path` ([nextflow#7667](https://github.com/nextflow-io/nextflow/issues/7667)).
- **List a file once per closure.** Later entries win. Same-named files from different processes need separate outputs ([nextflow#6617](https://github.com/nextflow-io/nextflow/issues/6617)).
- **Everything published lives in a record.** Per-sample files go in stage records, run-level files in small records grouped by output directory.
- **Field names are stable API.** They appear in output docs and `index` headers.

### Replacing `collectFile`

Write files with a small module that takes and returns a record, using `exec:` so no container is needed:

```nextflow
record MultiqcWriteFileInput {
    id:      String
    meta:    Map
    name:    String
    content: String
}

record MultiqcWriteFileResult {
    id:   String
    meta: Map
    file: Path
}

process MULTIQC_WRITE_FILE {
    input:
    sample: MultiqcWriteFileInput

    output:
    record(id: sample.id, meta: sample.meta, file: file(sample.name)) as MultiqcWriteFileResult

    exec:
    task.workDir.resolve(sample.name) << sample.content
}
```

A sibling module concatenates a list of tables the same way, skipping repeated header rows.

### Workflows

- To `emit:` one result, assign the variable without `def`, or emit an expression.
- A `take:` typed `Value<...>` should receive a real value channel, including in tests (`channel.value(...)`).
- Annotate process call results with the result type, so the shape is visible at the call site: `def ch_sorted: Channel<SamtoolsSortResult> = SAMTOOLS_SORT(...)`.
- For cross-cutting collectors such as a report, have each stage emit a small record of the files it contributes (`MultiqcFiles`), and let the collector take the channel of them.
- Keep `nextflow.enable.types = true` in each script that uses typed processes or workflows.
- Sort file lists passed to MultiQC so report inputs are deterministic.

### Tests

Build test inputs with `record(...)`, use `channel.value(...)` for `Value` inputs, and assert on named fields:

```groovy
when {
    workflow {
        """
        input[0] = Channel.of(record(
            id: 'test',
            meta: [ id:'test', single_end:false ],
            raw_bams: file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/bam/test.single_end.bam', checkIfExists: true)
        ))
        input[1] = channel.value(file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/genome.fasta', checkIfExists: true))
        """
    }
}
then {
    assert workflow.success
    assert workflow.out[0][0].bam ==~ ".*.bam"
}
```

## Checking a migration

- Run `nextflow lint`, then the real tests early. Lint does not check that closures, templates and config expressions resolve at runtime.
- Compare with the pre-migration pipeline: same published file tree, same md5sums. Explain every difference by the migration, and update snapshots only afterwards.
- Keep edits to vendored nf-core modules minimal so the diff against upstream stays reviewable.
- After scripted renames, grep for leftovers. Configs and templates are not linted.

## Practical notes

- Static typing needs the strict parser, which is the default from Nextflow 26.04, as well as `nextflow.enable.types`.
- nf-schema's `samplesheetToList` runs in typed scripts, but `nextflow lint` rejects it ([nextflow#7720](https://github.com/nextflow-io/nextflow/issues/7720)). A small untyped helper file is a pragmatic home.
- Do not run `nf-core modules update` or `subworkflows update` on edited vendored components until the changes are upstreamed.
- `task.ext` and `meta` interact with record input naming. An upstream experiment replaces them with plain process inputs ([nf-core/methylseq#625](https://github.com/nf-core/methylseq/pull/625)).
- Declare numeric field types to match what tools produce, and leave values unconverted if their text appears in an output file.
- Workflow-level type checks are still catching up ([nextflow#7693](https://github.com/nextflow-io/nextflow/issues/7693)). Publishing directories to object storage ([nextflow#7398](https://github.com/nextflow-io/nextflow/issues/7398)) is not yet verified.

User-visible changes to plan for: `process.withName:<NAME>.publishDir` overrides in user configs stop working (breaking), directories are only created when a file is published into them, and the minimum Nextflow version rises.

## Changes needed in nf-core/tools

Vendored components edited for this migration expose gaps in current tooling:

- **Dependency inference:** `nf-core subworkflows install` infers dependencies from `include` paths and lower-cases the symbol, so `include { SamtoolsSortResult } from '.../samtools/sort/main'` is read as a module dependency. Type imports need to be recognised and ignored, or resolved as types.
- **Patching new files:** `nf-core modules patch` and `subworkflows patch` record a new file (such as `types.nf`) as a bare "was created" note without its content. A later `update` could delete the file and leave the `include` behind (related: [nf-core/tools#3157](https://github.com/nf-core/tools/issues/3157)).
- **Local-copy lint:** components edited without patch files fail `check_local_copy`. This is expected until the changes are upstreamed.
- **Shared types:** a shared `types.nf` is not a module or subworkflow, so there is no way yet to install, version or update it.
- **Component metadata:** `meta.yml` needs a way to describe record inputs and outputs. Whether lint checks these against typed process signatures is not yet confirmed.

## Open questions

- How should nf-core handle `task.ext.when` alongside typed `when:`?
- How do record-typed modules land in nf-core/modules, and where do shared types live: one shared file or a `types.nf` per component?
- Should record types be cast at all, given [nextflow#7680](https://github.com/nextflow-io/nextflow/issues/7680)?
- How should new files in vendored components (such as `types.nf`) be tracked without patch files?
- Should `task.ext` and `meta` give way to plain process inputs?
- Which conventions above should become specifications, and which stay guidance?
