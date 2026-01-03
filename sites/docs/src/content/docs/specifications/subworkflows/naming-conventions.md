---
title: Naming conventions
subtitle: Naming conventions for nf-core Nextflow DSL2 subworkflows
weight: 2
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Name format of subworkflow files

Subworkflow names SHOULD be based on their module composition.

For short chains of modules without conditional logic, use the format `<file_type>_<operation_1>_<operation_n>_<tool_1>_<tool_n>`. For example, `bam_sort_stats_samtools` where `bam` = `<file_type>`, `sort` = `<operation>` and `samtools` = `<tool>`. Exclude routine operations from the name (for example, indexing after BAM creation). Collapse directly related steps into a general name. For example, if a binning tool has three required steps (`<tool> split`, `<tool> calculate`, `<tool> merge`) to perform contig binning, collapse these into one (for example, `fasta_binning_concoct`, rather than `fasta_split_calculate_merge_concoct`).

The above naming scheme may not be appropriate if:

- The subworkflow has a large number of steps (excluding routine operations)
- The sequence of steps differs based on input arguments
- The module complement is likely to change over time

In these cases, name the subworkflow according to its purpose and logical operations rather than the module complement. For example, a subworkflow that takes FASTQ files, performs multiple QC checks, applies a trimming operation, filters and sets strandedness would be named `fastq_qc_trim_filter_setstrandedness`. This tells users the input and logical steps without including every conditional or module in the name.

Use all lowercase for the directory structure. For example, [`subworkflows/nf-core/bam_sort_stats_samtools/`](https://github.com/nf-core/modules/tree/master/subworkflows/nf-core/bam_sort_stats_samtools/).

For assistance with naming subworkflows, particularly complex subworkflows, post on the [nf-core Slack `#subworkflows` channel](https://nfcore.slack.com/channels/subworkflows) (you can join with [this invite](https://nf-co.re/join/slack)) to discuss options.

## Name format of subworkflow parameters

Parameter names MUST use the `snake_case` convention.

## Name format subworkflow functions

Function names MUST use the `camelCase` convention.

## Name format subworkflow channels

Channel names MUST use the `snake_case` convention and be lowercase.

## Input channel name structure

Input channel names SHOULD signify the input object type.
Prefix single value input channels with `val_`. Prefix input channels with multiple elements (e.g., meta map and file) with `ch_`.

## Output channel name structure

Output channel names SHOULD be based on the major output file of that channel.
For example, emit an output channel of `[[meta], bam]` as `bam`, not `ch_bam`.
This provides more intuitive use of these output objects downstream with the `.out` attribute.
