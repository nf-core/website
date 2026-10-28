---
title: "Upcoming changes to nf-core/modules"
subtitle: "Announcement of the upcoming changes to the guidelines when migrating to types and records in nf-core/modules"
pubDate: 2026-10-28T09:00:00+01:00
headerImage: "https://images.unsplash.com/photo-1582730147924-d92f4da00252"
headerImageAlt: "A bunch of records laid out onto a black surface."
authors:
  - nvnieuwk
label:
  - modules
---

TODO add more memes

During the nf-core hackathon in October 2026, two hackathon groups merged together to form one big group to tackle the upcoming changes caused by the introduction of static typing and record types and decide on the guidelines related to these topics.

This blog post gives an overview of all the upcoming changes that are coming to nf-core/modules along with the reasoning behind every change.

## Singular inputs and outputs

Previously, multiple lines were allowed in the input and output block of a process. In fact we even enforced that the output block must contain one output per file produced by the process. This will change as of now by enforcing each module to take one [Record](https://docs.seqera.io/nextflow/reference/stdlib-types/record) as an input and to emit one [Record](https://docs.seqera.io/nextflow/reference/stdlib-types/record) as the output.

This change will be enforced due to many reasons, these are the most prominent ones:

1. Having a singular input to the module will allow full flexibility for each process run. Due to this change we will never have to do massive multiMap operations to make sure that secondary inputs are in sync with the actual sample data.
2. Having a singular output with record types will make it easier to subset the record type and merge it back into another record type. This is especially useful for processes that emit a lot of files as there will be no need for multi-line joins anymore.
3. Due to this the whole idea of having a dataflow becomes even more clear.

## Composition of the input record type

The input record should consist of several required and optional typed fields. A record type is a `Map`-like structure that consists of a predefined set of keys. The record is not limited to only these keys, but it will validate if these keys exists (if needed) and if they have the right type. To type a specific field, put a `:` followed by the expected type after the field key name. A record with typed fields can look like this:

```groovy
input:
record(
	id: String,
	input_file: Path,
	input_index: Path?
)
```

To annotate a field as optional, put a `?` after the type definition of the specific field.

Read more about typing in the [Nextflow docs](https://docs.seqera.io/nextflow/reference/stdlib-types).

One special case of typing is file inputs:

- In case that only one file is expected use the `Path` type
- In case multiple files are expected or possible, use `List<Path>` if the ordering matters and use `Bag<Path>` if the ordering doesn’t matter. `List<Path>` however is going to be the most used type so use that one if in doubt.

The input record types should consist of the following fields:

1. All required input values needed to run the process. One big change on this front is that mandatory fields from the meta map should now also be values of the record type. This is the `id` field for every process and optionally `single_end` wherever it is needed. This will increase the transparency of each process meaning that it will no longer have a silent requirement on meta fields.
2. All files (required and optional) that the process can use to run. These should all have some variation of the `Path` type as described above.
3. The `meta` [Map](https://docs.seqera.io/nextflow/reference/stdlib-types/map) will become the `meta` [Record](https://docs.seqera.io/nextflow/reference/stdlib-types/record). This field will now be fully optional and is mainly there to help pipeline developers pass additional pipeline-related information to the process. This can be used by pipeline users in `ext.args` and will still be passed to the output so it can also be used for downstream logic in the pipeline (like joining)
4. Optional arguments strings. There should be one argument string for every tool used by the process, for which the name should follow the following structure: `<tool_name(_subtool_name)>_args` (this is the same notation as for the module name). This will replace the usage of `ext.args` for pipeline **developers**. `ext.args` will still be available for pipeline users. In case the same tool is used twice or more in a process, a number can be added after the tool name in the field name.
5. The optional `output_prefix` field. This field will replace the usage of `ext.prefix` for pipeline **developers**. `ext.prefix` will still be available for pipeline users. In the rare case that multiple prefixes are needed an additional `output_prefix2` field can also be added.

The main logic behind most of these changes is to make the modules compatible with new and upcoming Nextflow features. Due to these changes a config file is no longer needed to configure the process properly. This makes the creation of meta pipelines a lot easier. Additionally by removing the requirement of the meta [Record](https://docs.seqera.io/nextflow/reference/stdlib-types/record), the usage of `nextflow module run` will also be a lot easier.

In case the inputs need to be staged, please move these to the new `stage` block of the process. More information on this can be found in the [docs](https://docs.seqera.io/nextflow/reference/process/inputs-outputs-typed#stage-directives)

### Example

This example will show the changes that need to be done for the `SAMTOOLS_STATS` module.

This module takes two required files called `input` and `input_index` and two optional files called `fasta` and `fai`. These files are split over two input lines, that need to be merged into one record. `meta2` is removed because the two optional files are now part of the main input entry. The `id` field has been added as a required value and the required `meta` Map was turned into an optional Record. Finally, `samtools_stats_args` and `output_prefix` have been added as optional String values.

The input section of a module will change like this:

```diff lang="groovy"
- input:
- tuple val(meta), path(input), path(input_index)
- tuple val(meta2), path(fasta), path(fai)
+ input:
+ record(
+ 	id: String,
+ 	meta: Record?,
+ 	input: Path,
+ 	input_index: Path,
+ 	fasta: Path?,
+ 	fai: Path?,
+ 	samtools_stats_args: String?,
+ 	output_prefix: String?
+ )
```

So, with inputs out of the way, let's move on to the module outputs!

## Composition of the output record

The record definition in the output block is a bit different when compared to the record definition in the input block. There is no option to provide the expected type of the output since Nextflow already derives this based on the type of the value passed to each record field. Read more about this in the [docs](https://docs.seqera.io/nextflow/reference/process/inputs-outputs-typed#outputs)

The output record should consist of the following fields:

1. The `id` field given via the input record
2. The `meta` Record given via the input record (this will be `null` in case no meta was passed)
3. All possible file outputs of the process correctly configured using the `file()` or `files()` functions (this includes required and optional files, use `optional: true` for optional files)
4. All version outputs should be moved to the `topics` block

### Example

Let's return to our previous example of `SAMTOOLS_STATS` and take a look at what needs to be done to update the output block of that module.

The module outputs one file called `stats` and one version of the `samtools stats` tool (as an output and a topic). In this case no multiple outputs should be merged because the module only outputs one file. The `id` field was added to the output record, which directly uses the `id` value passed to the input of the module. Finally the versions output has also been converted to a record type and has been moved to the `topic` block.

The output section of the module will change like this:

```diff lang="groovy"
- output:
- tuple val(meta), path("*.stats"), emit: stats
- tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), emit: versions_samtools, topic: versions
+ output:
+ record(
+ 	id: id,
+ 	meta: meta,
+ 	stats: file("*.stats")
+ )
+ 
+ topic:
+ record(process: "${task.process}", tool: 'samtools', version: eval('samtools version | sed "1!d;s/.* //"')) >> 'versions'
```

As seen in the input of the module, some extra inputs were added that have a similar function to the current `ext` variables. Let's take a look how to properly handle these variables after the migration.

## Changes to `ext` variables

Instead of defining the `prefix` and `args` variables in the `script` and `stub` block separately, they must now be defined under the container declaration in the process code. This way the declaration happens only once in the process, making sure that there is more parity between the `script` and `stub` block. The definition of these variables now also needs to take into account that the `output_prefix` and `<tool>_args` options are given in the input. This is the order of priority for each option where one of these variables can be declared (with the first one taking priority over all others, etc):

1. The `ext.<var>` declaration in a configuration file passed by the pipeline **user**
2. The value passed in the input record of the process by the pipeline **developer**
3. The default of the value (for example `id` for the `prefix`)

`ext.when` should also be removed since it's no longer supported in statically typed modules, the new `process.when` option can be used instead. This configuration option should only be reserved for pipeline **users**. Pipeline developers should setup their logic in the workflow itself using `if` statements and channel operator logic.

### Example

Let's return to `SAMTOOLS_STATS` once again and take a look at how to handle the `ext` variables there.

The handling of `ext` variables will change like this:

```diff lang="groovy"
process SAMTOOLS_STATS {
	...

+    ext [
+        prefix: { output_prefix ?: "${id}" },
+        args: { samtools_stats_args ?: '' }
+    ]

	input:
	...

	output:
	...

+   topic:
+   ...

-	when: task.ext.when == null || task.ext.when

	script:
-	def args = task.ext.args ?: ''
-	def prefix = task.ext.prefix ?: "${meta.id}"
-	def reference = fasta ? "--reference ${fasta}" : ""
	"""
	samtools \\
		stats \\
-		${args} \\
+       ${task.ext.args} \\
		--threads ${task.cpus} \\
		${reference} \\
		${input} \\
-		> ${prefix}.stats
+       > ${task.ext.prefix}.stats
	"""

	stub:
-	def prefix = task.ext.prefix ?: "${meta.id}"
	"""
-	touch ${prefix}.stats
+   touch ${task.ext.prefix}.stats
	"""
}
```

## Full example

This is full example of the `SAMTOOLS_STATS` module.

```groovy title="modules/nf-core/samtools/stats/main.nf"
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
