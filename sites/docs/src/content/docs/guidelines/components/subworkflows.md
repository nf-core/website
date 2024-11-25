---
title: Subworkflow Specifications
subtitle: Specifications for writing nf-core Nextflow DSL2 subworkflows
markdownPlugin: addNumbersToHeadings
weight: 20
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## General

### Minimum subworkflow size

Subworkflows should combine tools that make up a logical unit in an analysis step.
A subworkflow must contain at least two modules.

### Version reporting channel

Each `subworkflow` emits a channel containing all `versions.yml` collecting the tool(s) versions.
They MUST be collected within the workflow and added to the output as `versions` :

```bash
take:
  input

main:

  ch_versions = Channel.empty()

  FASTQC(input)

  ch_versions = ch_versions.mix(FASTQC.out.versions())

emit:
  versions = ch_versions
```

## Naming conventions

### Name format of subworkflow files

Choose an appropriate name for your subworkflow related to its module composition.

For short chains of modules without conditional logic, the naming convention should be of the format `<file_type>_<operation_1>_<operation_n>_<tool_1>_<tool_n>` e.g. `bam_sort_stats_samtools` where `bam` = `<file_type>`, `sort` = `<operation>` and `samtools` = `<tool>`. Not all operations are required in the name if they are routine (e.g. indexing after creation of a BAM). Operations can be collapsed to a general name if the steps are directly related to each other. For example if in a subworkflow, a binning tool has three required steps (e.g. `<tool> split`, `<tool> calculate`, `<tool> merge`) to perform an operation (contig binning) these can be collapsed into one (e.g. `fasta_binning_concoct`, rather than `fasta_split_calculate_merge_concoct`).

If a subworkflow has a large number of steps discounting routine operations, if the sequence of steps differs dependent on input arguments, or if the module complement is likely to change over time, the above naming scheme will not be appropriate. In this case it will be more useful to potential users of your subworkflow to name the it according to its purpose and logical operations, rather than the module complement. For example, a subworkflow that takes FASTQ files, peforms multiple QC checks, applies a user defined trimming operation, filters and sets a strandedness, would be named something like 'fastq_qc_trim_filter_setstrandedness'. This tells users what the the input is, and the logical steps involved, without trying to shoehorn the conditional logic or very long sequences of modules into the name.

Whatever name is used, the directory structure for the subworkflow name must be all lowercase e.g. [`subworkflows/nf-core/bam_sort_stats_samtools/`](https://github.com/nf-core/modules/tree/master/subworkflows/nf-core/bam_sort_stats_samtools/).

If in doubt regarding what to name your subworkflow, and always for the more complex type of subworkflow described above, please contact us on the [nf-core Slack `#subworkflows` channel](https://nfcore.slack.com/channels/subworkflows) (you can join with [this invite](https://nf-co.re/join/slack)) to discuss possible options.

### Name format of subworkflow parameters

All parameter names MUST follow the `snake_case` convention.

### Name format subworkflow functions

All function names MUST follow the `camelCase` convention.

### Name format subworkflow channels

Channel names MUST follow `snake_case` convention and be all lower case.

### Input channel name structure

Input channel names SHOULD signify the input object type.
For example, a single value input channel will be prefixed with `val_`, whereas input channels with multiple elements (e.g. meta map + file) should be prefixed with `ch_`.

### Output channel name structure

Output channel names SHOULD only be named based on the major output file of that channel (i.e, an output channel of `[[meta], bam]` should be emitted as `bam`, not `ch_bam`).
This is for more intuitive use of these output objects downstream with the `.out` attribute.

## Input/output options

### Required input channels

Input channel declarations MUST be defined for all _possible_ input files that will be required by the subworkflow (i.e. both required and optional files) within the `take` block.

### Required output channels

Named file extensions MUST be emitted for ALL output channels e.g. `path "*.txt", emit: txt`.

### Optional inputs

Optional inputs are not currently supported by Nextflow.
However, passing an empty list (`[]`) instead of a file as a subworkflow parameter can be used to work around this issue.

## Subworkflow parameters

### Usage of parameters

Named `params` defined in the parent workflow MUST NOT be assumed to be passed to the subworkflow to allow developers to call their parameters whatever they want.
In general, it may be more suitable to use additional `input` value channels to cater for such scenarios.

## Documentation

### Code comment of channel structure

Each input and output channel SHOULD have a comment describing the output structure of the channel e.g

```nextflow
input:
ch_reads // channel: [mandatory] meta, reads
val_sort // boolean: [mandatory] false
<...>

emit:
bam = SAMTOOLS_VIEW.out.bam // channel: [ val(meta), path(bam) ]
versions = ch_versions      // channel: [ path(versions.yml) ]
```

### Meta.yml documentation of channel structure

Each input and output channel structure SHOULD also be described in the `meta.yml` in the description entry.

```text
description: |
  Structure: [ val(meta), path(tsv)]
  (Sub)contig coverage table
```

## Testing

### All output channels must be tested

All output channels SHOULD be present in the nf-test snapshot file, or at a minimum, it MUST be verified that the files exist.

### Tags

Tags for any dependent modules MUST be specified to ensure changes to upstream modules will re-trigger tests for the current subworkflow.

```groovy
tag "subworkflows"
tag "subworkflows_nfcore"
tag "<subworkflow_name>"
tag "<tool>" // Add each tool as a separate tag
tag "<tool>/<subtool>" // Add each subtool as a separate tag
```

### `assertAll()`

The `assertAll()` function MUST be used to specify an assertion, and there MUST be a minimum of one success assertion and versions in the snapshot.

### Assert each type of input and output

There SHOULD be a test and assertions for each type of input and output.

[Different assertion types](/docs/contributing/nf-test/assertions) should be used if a straightforward `workflow.out` snapshot is not feasible.

:::tip
Always check the snapshot to ensure that all outputs are correct!
For exmaple, make sure there are no md5sums representing empty files.
:::

### Test names

Test names SHOULD describe the test dataset and configuration used. some examples below:

```groovy
test("homo_sapiens - [fastq1, fastq2] - bam")
test("sarscov2 - [ cram, crai ] - fasta - fai")
test("Should search for zipped protein hits against a DIAMOND db and return a tab separated output file of hits")
```

### Input data

Input data SHOULD be referenced with the `modules_testdata_base_path` parameter:

```groovy
file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam', checkIfExists: true)
```

### Configuration

Subworkflow nf-tests SHOULD use a single `nextflow.config` to supply `ext.args` to a subworkflow. They can be defined in the `when` block of a test under the `params` scope.

```groovy {4-7} title="main.nf.test"
config './nextflow.config'

when {
  params {
    moduleA_args = '--extra_opt1 --extra_opt2'
    moduleB_args = '--extra_optX'
  }
  process {
    """
    input[0] = [
      [ id:'test1', single_end:false ], // meta map
      file(params.modules_testdata_base_path + 'genomics/prokaryotes/bacteroides_fragilis/genome/genome.fna.gz', checkIfExists: true)
    ]
    """
  }
}
```

```groovy {3,6} title="nextflow.config"
process {
  withName: 'MODULEA' {
    ext.args = params.moduleA_args
  }
  withName: 'MODULEB' {
    ext.args = params.moduleB_args
  }
}
```

No other settings should go into this file.

:::tip
Supply the config only to the tests that use `params`, otherwise define `params` for every test including the stub test.
:::

## Misc

### General module code formatting

All code MUST be aligned to follow the '[Harshil Alignment™️](/docs/contributing/code_editors_and_styling/harshil_alignment)' format.
