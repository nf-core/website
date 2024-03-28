---
title: Subworkflow Specifications
subtitle: Specifications for writing nf-core Nextflow DSL2 subworkflows
---

## New subworkflow specifications

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

### General

1. Subworkflows should combine tools that make up a logical unit in an analysis step. A subworkflow must contain at least two modules.

2. Each `subworkflow` emits a channel containing all `versions.yml` collecting the tool(s) versions. They MUST be collected within the workflow and added to the output as `versions` :

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

### Naming conventions

1. The directory structure for the subworkflow name must be all lowercase e.g. [`subworkflows/nf-core/bam_sort_stats_samtools/`](https://github.com/nf-core/modules/tree/master/subworkflows/nf-core/bam_sort_stats_samtools/). The naming convention should be of the format `<file_type>_<operation_1>_<operation_n>_<tool_1>_<tool_n>` e.g. `bam_sort_stats_samtools` where `bam` = `<file_type>`, `sort` = `<operation>` and `samtools` = `<tool>`. Not all operations are required in the name if they are routine (e.g. indexing after creation of a BAM). Operations can be collapsed to a general name if the steps are directly related to each other. For example if in a subworkflow, a binning tool has three required steps (e.g. `<tool> split`, `<tool> calculate`, `<tool> merge`) to perform an operation (contig binning) these can be collapsed into one (e.g. `fasta_binning_concoct`, rather than `fasta_split_calculate_merge_concoct`). If in doubt regarding what to name your subworkflow, please contact us on the [nf-core Slack `#subworkflows` channel](https://nfcore.slack.com/channels/subworkflows) (you can join with [this invite](https://nf-co.re/join/slack)) to discuss possible options.

2. All parameter names MUST follow the `snake_case` convention.

3. All function names MUST follow the `camelCase` convention.

4. Channel names MUST follow `snake_case` convention and be all lower case.

5. Input channel names SHOULD signify the input object type. For example, a single value input channel will be prefixed with `val_`, whereas input channels with multiple elements (e.g. meta map + file) should be prefixed with `ch_`.

6. Output channel names SHOULD only be named based on the major output file of that channel (i.e, an output channel of `[[meta], bam]` should be emitted as `bam`, not `ch_bam`). This is for more intuitive use of these output objects downstream with the `.out` attribute.

### Input/output options

1. Input channel declarations MUST be defined for all _possible_ input files that will be required by the subworkflow (i.e. both required and optional files) within the `take` block.

2. Named file extensions MUST be emitted for ALL output channels e.g. `path "*.txt", emit: txt`.

3. Optional inputs are not currently supported by Nextflow. However, passing an empty list (`[]`) instead of a file as a subworkflow parameter can be used to work around this issue.

### Subworkflow parameters

1. Named `params` defined in the parent workflow MUST NOT be assumed to be passed to the subworkflow to allow developers to call their parameters whatever they want. In general, it may be more suitable to use additional `input` value channels to cater for such scenarios.

### Documentation

1. Each input and output channel SHOULD have a comment describing the output structure of the channel e.g

   ```nextflow
   input:
   ch_reads // channel: [mandatory] meta, reads
   val_sort // boolean: [mandatory] false
   <...>

   emit:
   bam = SAMTOOLS_VIEW.out.bam // channel: [ val(meta), path(bam) ]
   versions = ch_versions      // channel: [ path(versions.yml) ]
   ```

2. Each input and output channel structure SHOULD also be described in the `meta.yml` in the description entry.

   ```text
   description: |
     Structure: [ val(meta), path(tsv)]
     (Sub)contig coverage table
   ```
