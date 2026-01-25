---
title: ext arguments
subtitle: Configure tool-specific options in nf-core components
shortTitle: ext arguments
---

The `ext` directive in Nextflow allows you to configure tool-specific options for nf-core modules and subworkflows.
These special process directives let you insert strings directly into module scripts without modifying the code.

## Understanding ext properties

The most common `ext` property is `ext.args`, which enables you to inject command-line arguments into module scripts.
For modules that run multiple tools, you can use numbered variants like `ext.args2` and `ext.args3`.

This approach provides flexibility whilst maintaining the integrity of the original module code.

## Basic example

To add custom arguments to a module, configure the process in your `nextflow.config`:

```groovy
process {
  withName: 'TOOL_SUBTOOL' {
    ext.args = '-T -K'
  }
}
```

This configuration transforms the module script from:

```bash
tool subtool $args $bam > ${prefix}.log
```

Into the executed command:

```bash
tool subtool -T -K test.bam > test.log
```

The `ext.args` value replaces the `$args` variable in the module script, and the `$bam` and `${prefix}` variables are populated from the module inputs and configuration.

## Common `ext` keys

Different `ext` keys serve specific purposes in module configuration:

| Key          | Purpose                                                          |
| ------------ | ---------------------------------------------------------------- |
| `ext.args`   | Additional command arguments for the primary tool                |
| `ext.args2`  | Secondary argument set for the second tool in multi-tool modules |
| `ext.args3`  | Tertiary argument set for the third tool in multi-tool modules   |
| `ext.prefix` | Custom prefix for output file names                              |

:::warning
The numeric order of `args` keys must match the order of tools as they appear in the module script.
:::

## Advanced usage

You can use Groovy expressions within `ext.args` to conditionally set arguments based on pipeline parameters or task inputs.

### Parameter inputs

Set `ext.args` based on parameter settings:

```groovy
process {
    withName: '.*:FASTQC_UMITOOLS_TRIMGALORE:UMITOOLS_EXTRACT' {
        ext.args   = [
                params.umitools_extract_method ? "--extract-method=${params.umitools_extract_method}" : '',
                params.umitools_bc_pattern     ? "--bc-pattern='${params.umitools_bc_pattern}'" : '',
                params.umitools_bc_pattern2    ? "--bc-pattern2='${params.umitools_bc_pattern2}'" : ''
            ].join(' ').trim()
        ]
    }
}
```

### Task inputs

Set `ext.prefix` based on task inputs:

```groovy
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:.*:BAM_SORT_SAMTOOLS:BAM_STATS_SAMTOOLS:.*' {
        ext.prefix = { "${meta.id}.sorted.bam" }
        publishDir = [
            path: { "${params.outdir}/${params.aligner}/samtools_stats" },
            mode: params.publish_dir_mode,
            pattern: "*.{stats,flagstat,idxstats}"
        ]
    }
}
```

### Parameter and task inputs

Set `ext.args` based on both parameters and task inputs:

```groovy
process {
    withName: '.*:DEDUP_UMI_UMITOOLS_GENOME:UMITOOLS_DEDUP' {
        ext.args   = { [
            meta.single_end                 ? '' : '--unpaired-reads=discard --chimeric-pairs=discard',
            params.umitools_grouping_method ? "--method='${params.umitools_grouping_method}'" : '',
            params.umitools_umi_separator   ? "--umi-separator='${params.umitools_umi_separator}'" : ''
        ].join(' ').trim() }
        ext.prefix = { "${meta.id}.umi_dedup.sorted" }
    }
}
```
