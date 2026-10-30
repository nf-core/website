---
title: Static types and records
subtitle: Write modules with static types and records
markdownPlugin: addNumbersToHeadings
shortTitle: Static types and records
weight: 10
---

The keywords "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

These specifications apply to modules that enable [static types](https://docs.seqera.io/nextflow/process-typed) with `nextflow.enable.types = true`.
Static types and records require Nextflow 26.10.0 or later.
Where a rule on this page conflicts with the [general module specifications](/docs/specifications/components/modules/general), this page takes precedence for modules that use static types.

To convert an existing module, see [Migrating modules to static types and records](/docs/developing/migration-guides/static-types).

## One input record and one output record

All processes MUST have one record as the input and MUST have one record as the output.

## Input record structure

The input record MUST contain the required fields and MAY contain the optional fields listed below.

### Required fields

The input record MUST include the following required fields:

- Any required input values used to control internal module logic, including the `id` field
- Any required input files of the process

### Optional fields

The following fields SHOULD be part of the input record, and SHOULD be optional for the module user.
Mark a field as optional by adding `?` after its type, for example `meta: Record?`.

- `meta`: A `Record` that passes metadata from the input of the module to the output of the module.
  The process code MUST NOT use the content of `meta`.
  Pipeline users MAY use it in `ext` configuration options.
- `<tool_name>_args` or `<tool_name>_<subtool_name>_args`: A `String` that contains the non-mandatory arguments of the tool.
  Specify one field per tool, replacing each `ext.argsN`.
  Name the field following the [module naming conventions](/docs/specifications/components/modules/naming-conventions).
  For repeated calls of the same tool, add a number, for example `samtools_view1_args`.
- `output_prefix`: A `String` that specifies the output file prefix.
  In the rare case where a module needs a second prefix, add `output_prefix2`.

### File types

File fields in the input record SHOULD use one of the following types:

| Type         | Description                                       |
| ------------ | ------------------------------------------------- |
| `Path`       | A single file                                     |
| `List<Path>` | One or more files, where the order matters        |
| `Bag<Path>`  | One or more files, where the order doesn't matter |

## Output record structure

The output record MUST contain the following fields:

- `id`: The input identifier from the input record
- `meta`: The optional `meta` record from the input record, used to supply sample-level metadata beyond `id`.
  The output record MUST always emit `meta`.
- All possible file outputs of the process

## Version reporting

All modules MUST report versions in the `topic:` block using a record with three fields: `process`, `tool`, and `version`.
Send the record to the `versions` topic:

```groovy title="main.nf"
topic:
record(process: "${task.process}", tool: 'samtools', version: eval('samtools version | sed "1!d;s/.* //"')) >> 'versions'
```

## `ext` directive

Modules MUST specify all `ext` variables as an independent `ext` directive, instead of with `def prefix = ...` or `def args = ...` in the `script:` block.
Each `ext` variable SHOULD be defined as a closure, so that it evaluates correctly in all cases.
Use `task.ext.args` and `task.ext.prefix` in the `script:` and `stub:` blocks.

```groovy title="main.nf"
ext [
    prefix: { output_prefix ?: "${id}" },
    args: { samtools_stats_args ?: '' }
]
```

See [nf-core/modules#8808](https://github.com/nf-core/modules/pull/8808/changes) for the origin of this pattern.

## Conditional execution

Modules MUST NOT include a `when:` block that reads `task.ext.when`.
If a module needs conditional execution, the calling workflow decides whether to run the module, either with an `if` statement around the call or with a `filter` on the input channel.
Pipelines MUST NOT set `process.when` in `modules.config`.

## Arguments and prefixes in `modules.config`

Pipelines MAY move `ext.args` from `modules.config` to the `<tool_name>_args` field of the input record.
Pipeline users can still override the arguments by supplying their own config with `ext.args` set.

Similarly, pipelines MAY move `ext.prefix` from `modules.config` to the `output_prefix` field of the input record.

## Examples

The following examples show a module before and after conversion to static types and records, and a module that stages a script file.

### Module without static types and records

The following `samtools/stats` module uses tuple inputs and outputs, `def` variables for `ext` values, and a `when:` block:

```groovy title="main.nf"
process SAMTOOLS_STATS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(input), path(input_index)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    samtools \\
        stats \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reference} \\
        ${input} \\
        > ${prefix}.stats
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stats
    """
}
```

### Module with static types and records

The same module with static types uses one input record, one output record, an `ext` directive, and a `topic:` block:

```groovy title="main.nf"
nextflow.enable.types = true

process SAMTOOLS_STATS {
    tag "${id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    ext [
        prefix: { output_prefix ?: "${id}" },
        args: { samtools_stats_args ?: '' }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        input: Path,
        input_index: Path,
        fasta: Path?,
        fai: Path?,
        samtools_stats_args: String?,
        output_prefix: String?
    )

    output:
    record(
        id: id,
        meta: meta,
        stats: file("*.stats", optional: true)
    )

    topic:
    record(process: "${task.process}", tool: 'samtools', version: eval('samtools version | sed "1!d;s/.* //"')) >> 'versions'

    script:
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    samtools \\
        stats \\
        ${task.ext.args} \\
        --threads ${task.cpus} \\
        ${reference} \\
        ${input} \\
        > ${task.ext.prefix}.stats
    """

    stub:
    """
    touch ${task.ext.prefix}.stats
    """
}
```

### Module with a staged script

Modules with static types replace the `template` function with a `stage:` block that stages the script into the task directory.
The following `eautils/gtf2bed` module stages a Perl script from its `templates/` directory and runs it directly:

```groovy title="main.nf"
nextflow.enable.types = true

process EAUTILS_GTF2BED {
    tag "${id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/perl:5.26.2'
        : 'quay.io/biocontainers/perl:5.26.2'}"

    ext [
        prefix: { output_prefix ?: "${id}" },
        args: { eautils_gtf2bed_args ?: '' }
    ]

    input:
    record(
        id: String,
        meta: Record?,
        gtf: Path,
        eautils_gtf2bed_args: String?,
        output_prefix: String?
    )

    stage:
    stageAs file("${moduleDir}/templates/gtf2bed.pl"), '*'

    output:
    record(
        id: id,
        meta: meta,
        bed: file("*.bed", optional: true)
    )

    topic:
    record(process: "${task.process}", tool: 'perl', version: eval("perl -e 'print substr(\$^V, 1)'")) >> 'versions'

    script:
    """
    ./gtf2bed.pl ${task.ext.args} ${gtf} ${task.ext.prefix}
    """

    stub:
    """
    touch ${task.ext.prefix}.bed
    """
}
```
