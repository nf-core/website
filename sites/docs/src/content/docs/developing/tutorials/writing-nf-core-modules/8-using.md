---
title: "Chapter 8: Using nf-core modules in pipelines"
subtitle: "How to use nf-core modules in a pipeline"
shortTitle: "Chapter 8: Using modules"
---

This chapter explains how to install an nf-core module in a pipeline and what else to consider when using it.

It assumes you have already contributed your module to the nf-core/modules repository or a personal or organisation-specific module repository, and that the module has been merged in.

## Installing nf-core modules

### Pipelines using the nf-core template

If you generated your pipeline with [nf-core tools](https://nf-co.re/docs/nf-core-tools/pipelines/create) ([whether the pipeline is official or not](https://nf-co.re/docs/developing/pipelines/external-use)), you can install a module directly. Run the following from the root of the repository, or pass `--dir`:

```bash
nf-core modules install <toolname>/<subcommand>
```

This downloads the module files into the existing `modules/nf-core` directory. The console output prints an `include` line you can paste into your pipeline script.

### Custom pipelines

You can use the same command in a custom pipeline. On the first run you will be prompted to:

1. Select `pipeline` when asked whether this is a `pipeline` or `modules` repository.
2. Accept creation of an `nf-core.yaml` file. This supports adding more modules later.
3. Accept creation of a `modules.json` file. This tracks installed module versions for future updates.

The command creates `modules/nf-core/<toolname>/<subcommand>`, a `modules.json`, and an `.nf-core.yml` configuration file. The console output prints an `include` line for your pipeline script.

:::warning
If installing a second module returns `ERROR 'manifest.name'`, add a `nextflow.config` with a [`manifest` scope](https://www.nextflow.io/docs/latest/reference/config.html#manifest) that includes `name`, `description`, and `version`.
:::

## Using nf-core modules

You can use an installed nf-core module as you would any Nextflow module.

Paste the `include` line printed by the install command into the `.nf` file where you want to use the module:

```nextflow
include { SAMTOOLS_FASTA } from '../modules/nf-core/samtools/fasta/main'
```

Adjust the path if you are using the module from a subworkflow in a nested directory.

Invoke it as you would any Nextflow module:

```nextflow
SAMTOOLS_FASTA (ch_input_for_samtoolsfasta, val_interleave)
```

Most pipelines need a few additional adjustments, covered below.

### Module channel structures

Reconfigure your input channels to match nf-core conventions.

The most important convention is the use of [meta maps](https://nf-co.re/docs/developing/components/meta-map). If your pipeline does not use meta maps, add one before passing data to the module:

```nextflow
def val_interleave = false
ch_input_for_samtoolsfasta = ch_input
                              .map {fasta -> [[id: fasta.simpleName] ,fasta]}

SAMTOOLS_FASTA (ch_input_for_samtoolsfasta, val_interleave)
```

Here, the original pipeline channel contained only FASTA files. The module needs a meta map, so `.map {}` creates one with a single `id` attribute assigned from the file's [`simpleName`](https://www.nextflow.io/docs/latest/reference/stdlib.html#stdlib-types-path).

### Tags and labels

#### Tags

Most nf-core modules set a `tag` directive to produce a human-readable process description in the run log.

By default, all nf-core modules use the `id` attribute from the meta map as the tag value. Set `id` on every meta map passed to an nf-core module.

#### Labels

nf-core modules use standardised `label` directives that map to default memory, CPU, and wall time.

- **Pipelines from the nf-core template**: defaults are defined in [`conf/base.config`](https://github.com/nf-core/tools/blob/52e810986e382972ffad0aab28e94f828ffd509b/nf_core/pipeline-template/conf/base.config#L22-L54). You do not need to change anything unless you want to override defaults for your pipeline.
- **Custom pipelines**: add a Nextflow config file that defines default resources for the labels used by your installed modules.

For example, if an installed module uses `process_low`:

```nextflow
process {
  withLabel:process_low {
      cpus   = { 2     * task.attempt }
      memory = { 12.GB * task.attempt }
      time   = { 4.h   * task.attempt }
  }
}
```

:::tip
You can override default resources for a specific module using `withName:` in the same process scope.
:::

### The `modules.conf` file

Pass optional tool arguments to modules through a [`modules.config` file](https://github.com/nf-core/tools/blob/52e810986e382972ffad0aab28e94f828ffd509b/nf_core/pipeline-template/conf/modules.config), which the nf-core pipeline template uses to inject pipeline-level parameters.

For each installed module, specify a `withName:` entry with at least `ext.args` and `publishDir`.

For example:

```nextflow {3}
process{
    withName: FASTQC {
        ext.args = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc" },
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
}
```

This hardcodes FastQC's `--quiet` flag for every invocation and sets the default output directory under `--outdir`. It also excludes the required `versions.yml` file from publication, because nf-core pipelines render versions through [MultiQC](https://multiqc.info/).

Pass arguments dynamically through a closure to expose pipeline parameters to users:

```nextflow {3}
process{
    withName: FASTQC {
        ext.args = { params.fastqc_quiet_mode ? '--quiet' : '' }
        publishDir = [
            path: { "${params.outdir}/fastqc" },
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
}
```

Combine multiple values into the `args` string by joining a list with spaces:

```nextflow
process{
    withName: FASTQC {
        ext.args = { [
            params.fastqc_quiet_mode ? '--quiet' : '',
            "--kmers ${params.fastqc_kmers}"
          ].join(' ').trim() }
        publishDir = [
            path: { "${params.outdir}/fastqc" },
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
}
```

With these steps, you can now use nf-core modules in your pipeline and benefit from the standardised, reproducible, and tested modules the community maintains.
