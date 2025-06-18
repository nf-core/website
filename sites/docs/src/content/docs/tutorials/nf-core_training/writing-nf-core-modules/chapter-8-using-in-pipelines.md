---
title: "Chapter 8: Using nf-core modules in pipelines"
subtitle: "How to effectively install and use an nf-core module in nf-core and non-nf-core pipelines"
shortTitle: "Chapter 8: Using"
---

## Introduction

By this final chapter, you have written and tested your nf-core module and are ready to use it in a pipeline.
We assume you have contributed your module(s) to the nf-core/modules repository or a personal/organisation specific module GitHub repository, and it has been merged in and ready for use.

This chapter will explain how to install your module in an nf-core or your own custom pipeline, and what other wider considerations to take when using them in a pipeline context.

## Installing nf-core/modules

In this section we will describe what you need to install an nf-core module for use within a pipeline.

### Pipeline using the nf-core template

If you have generated your pipeline skeleton using [nf-core tools](https://nf-co.re/docs/nf-core-tools/pipelines/create) ([regardless if this is an official nf-core pipeline or not](https://nf-co.re/docs/guidelines/external_use)) you can directly install your module with no extra setup required!

Simply run the following command from the root of the repository (or use `--dir` to specify the root):

```bash
nf-core modules install <toolname>/<subcommand>
```

This will download the modules files from the nf-core/modules repository directly into the existing `modules/nf-core` directory of your pipeline.

The console output will also handily give you a line you can copy and paste into your pipeline script to include the module 💪.

### Custom pipelines

If you are using a custom pipeline, you can still use the `nf-core modules install` command to download the module files into your pipeline.

If you do this, the first time you will be asked:

1. Whether this is a `pipeline` or `modules` repository: select `pipeline`
2. If you are OK with creating an `nf-core.yaml` file: select `yes`
   - This helps in the future for adding more modules to your pipeline
3. If you are OK with creating a `modules.json`: select `yes`
   - This file is used by nf-core tools to inform you of any updates to the modules you have installed

This will create a directory and subdirectories within your pipeline repository called `modules/nf-core/<toolname>/<subcommand>` where you can see all the files you created.
You will also see the `modules.json` for recording the module versions of all nf-core modules you've installed, and also a `.nf-core.yml` configuration file required by nf-core tools.

The console output will also handily give you a line you can copy and paste into your pipeline script to include the module 💪.

:::warning
If you install a second module, and get an error of `ERROR    'manifest.name`, ensure you have a `nextflow.config` file with the [`manifest` scope](https://www.nextflow.io/docs/latest/reference/config.html#manifest) specified and that this includes the `name`, `description`, and `version` attributes.
:::

## Utilising nf-core/modules

Once your module is created, you can pretty much use it as with any form of Nextflow module.

Copy and paste the line that the install command printed to console into the relevant `.nf` file you want to use the module in:

```nextflow
include { SAMTOOLS_FASTA } from '../modules/nf-core/samtools/fasta/main'
```

You may need to update the path in cases where you are using subworkflow scripts nested in other directories.

You can invoke the module as you would with any Nextflow module:

```nextflow
SAMTOOLS_FASTA (ch_input_for_samtoolsfasta, val_interleave)
```

However, you may need to make some additional additional tweaks to your existing pipeline logic and files to make the module work correctly. These changes are described in the remaining sections of this chapter.

### Module channel structures

One of the most important things to remember when using an nf-core module is to reconfigure your input channels to ensure they match nf-core conventions.

The most important convention is the widespread use of [meta maps](https://nf-co.re/docs/contributing/components/meta_map).
If you do not already use meta maps in your pipeline, you might need to update the input channel going into the nf-core module.

A very simple example would be:

```nextflow
def val_interleave = false
ch_input_for_samtoolsfasta = ch_input
                              .map {fasta -> [[id: fasta.simpleName] ,fasta]}

SAMTOOLS_FASTA (ch_input_for_samtoolsfasta, val_interleave)
```

Here, the pipeline originally had input channel that just contains FASTA files.
However the module requires a meta map to be passed with the FASTA file.

To add a meta map, we use the `.map {}` operator to make a new meta map and emit this with the FASTA file.
When making the meta map, we make it with a single nf-core standard attribute `id`, which used in almost all nf-core modules.
In this example, we simply assign to the `id` attribute the [`simpleName`](https://www.nextflow.io/docs/latest/reference/stdlib.html#stdlib-types-path) of the file (i.e., the file name without the full path and without any file suffix).

By reconfiguring the channel to include a valid meta map, we ensure the input channel is compatible with the nf-core module.

### Tags and labels

#### Tags

Most nf-core modules also include a `tag` directive as you will have seen in the `main.nf` module.
This is used for providing a more 'human readable' process description when Nextflow prints the progress of a pipeline run.

By default, all nf-core modules assume there is a meta attribute called `id` (see above), and this is used as the `tag` value.
We highly recommend specifying the `id` attribute in all meta maps of all the input channels for all nf-core modules that use meta map.

#### Label

The other important component of nf-core modules is the standardised `labels` directive.
These are used to specify default memory, CPU, and wall time resources required for the process.

As we saw in the 'generating boilerplate files' chapter, we will have already specified which label we should use in the process.

In pipelines using the nf-core template, the specifications for each of these labels are already defined with defaults within the pipeline template (see the [`conf/base.config` file](https://github.com/nf-core/tools/blob/52e810986e382972ffad0aab28e94f828ffd509b/nf_core/pipeline-template/conf/base.config#L22-L54)).
Therefore you don't have to make any additional changes to your pipeline, except if you want to change these default pipelines for the specific context of your pipeline.

If you are using a custom pipeline, you will need to ensure you have a Nextflow config file somewhere in your repository, and that it is loaded somewhere in your pipeline script (either implicitly such as a `nextflow.config` file or explicitly imported via an `include` as appropriate).
This config needs to have a process scope with default resources for the label of your installed module specified.

For example, if your installed nf-core module has a label of `process_low`, one of the configs in your pipeline should have:

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
Remember that in addition to the 'defaults' provided by the nf-core module labels, you can always override default resources for a specific module using the `withName` scope within the same process scope as `withLabel`.
:::

### The `modules.conf` file

A final consideration is how to pass optional tool options and arguments to the command within the module.

In the nf-core pipeline template, we utilise a [`modules.config` file](https://github.com/nf-core/tools/blob/52e810986e382972ffad0aab28e94f828ffd509b/nf_core/pipeline-template/conf/modules.config) to 'inject' pipeline level parameters into a given module.

The defintion of the arguments to be injected can be placed in the same `process` block as the labels in any config file using or loaded in the pipeline.
We recommend that for any nf-core module you install, you should specify a `withName:` label for each process, and at a minimum specify a `ext.args` and `publishDir` directives for this option injection and default output directory publishing.

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

Here we specify `ext.args` to 'hardcode' FastQC's `--quiet` flag into all invocation of the FastQC process within the pipeline.
We also specify the default filepath where the output files will be stored in the results directory (here specified with `--outdir`), but to exclude the nf-core module required `versions.yml` file from being copied to the results directory (as this is normally rendered in nf-core pipelines with [MultiQC](https://multiqc.info/)).

If you want to give a user more freedom if what options and arguments are passed to the tool in the module such as via a pipeline-level `params.` parameter, you can dynamically pass these to `ext.args` within a closure.

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

Or even the following, to inject a multiple specific values into the `args` string (which is concatenated with spaces).

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

## Summary

In summary, at it's core using nf-core modules in your pipeline is as simple as installing the module and then invoking it in your pipeline script.
That said, there are a couple of additional considerations for non-nf-core template pipelines to take into account regarding input channels and additional configuration options to specify in Nextflow configuration files, which we have described here.

With these steps, you can now use nf-core modules in your pipeline, and benefit from the standardised, reproducible, and tested modules that the nf-core community has to offer.
