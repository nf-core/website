---
title: "Tutorial: nf-core DSL2 pipeline tips and tricks"
subtitle: Tips and tricks for writing or converting to nf-core DSL2 pipelines
---

## Introduction

nf-core adds a lot of additional 'infrastructure' and guidelines around standard Nextflow pipelines. Particularly with the conversion of DSL1 to the new DSL2 version of Nextflow, this can be at points quite daunting.

This page as a reference for a variety of tips, tricks, and guidance for writing new or converting to nf-core DSL2 pipelines.

## General

### Planning

Plan ahead! Check whether the modules you need for your pipeline are already available on [nf-core/modules](https://github.com/nf-core/modules)! If not, create the modules in advance. A tutorial on how to do this is available [here](dsl2_modules_tutorial.md).

## Modules tips

### Reusing test data

When testing your module, you will need to define test data check the module code is working.

You should as _far as possible_ re-use already existing data. As long a tool can accept the input and produce successful output with exit code `0`, this is sufficient - even if the output itself is nonsense.

> **Example**: [damageprofiler](https://github.com/Integrative-Transcriptomics/DamageProfiler) is a tool to generate DNA damage profiles that are specific to ancient DNA samples from BAM files. Therefore it is normally only run on ancient DNA samples. However, you do _not_ need to make new test data specific for the module that includes damage. Damageprofiler accepts standard BAM files (via htsjdk), and we can therefore supply the already existing _Homo sapiens_ modules test data even if this is from modern samples without damage. This is because damageprofiler will still successfully run and produce output even if there are no damage in the data.

### Debugging test data location

You will define your test commands for a DSL2 module under `tests/software/<module>/main.yml`. You specify which input data you will use for the test using the _nested_ `test_data` param, e.g.

```nextflow
file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]
```

This `params.test_data` is defined in the file [`modules/tests/config/test_data.config`](https://github.com/nf-core/modules/blob/master/tests/config/test_data.config).

If you have errors such as `input not found` or similar during testing, you can debug whether the defined nested param is correct by checking the corresponding config.

## Pipeline

### Converting existing DSL1 pipeline to DSL2

When converting an existing nf-core pipeline from DSL1 to DSL2, it's important that you don't remove the previous `master` that contains the current DSL1 release, nor work on your DSL2 port on `dev`, as you may need to make patch releases to the DSL1 version in the meantime.

Instead, create a new branch called e.g. `dsl2`. Pull this branch to your local machine, and then while in that branch delete the entire contents of the directory to start from scratch.

Once deleted, run `nf-core create` (ensuring you're using at least version 1.15, which includes the new DSL2 template), filling in the required information as with the original pipeline.

Finally, compare the `tree` for both DSL1 and DSL2 branches, and do `diff` on files individually to add in custom code from the pipeline. If in doubt, you can always check for any missing files or content by looking directly at the nf-core pipeline template

### Reusing a module multiple times in one workflow

To re-use a module multiple times, you can use Nextflow's [alias](https://www.nextflow.io/docs/latest/dsl2.html#module-aliases) functionality to specify different names of a different module.

You can do this in the `include` statements of the `IMPORT {NF-CORE,LOCAL} MODULES/WORKFLOWS` of the nf-core pipeline template.

For example to use the `SAMTOOLS FLAGSTAT` module twice:

```nextflow
/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_PRE } from '../modules/nf-core/modules/samtools/flagstat/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_POST } from '../modules/nf-core/modules/samtools/flagstat/main'
```

Where you then use the corresponding modules in the workflow itself using their aliases.

### Pipeline arguments and options

Specifying parameters to customise particular processes in nf-core DSL2 is rather more complex than DSL1, but at the same time is much more complex.

The following diagram describes how various parameters and options are inherited.

> !PLACEHOLDER!

In general there are two main entry points by users for customising options/parameters:

- A custom nextflow.config file
- Command line parameters

Command line parameters are generally should be used for parameters that are _required_ for the pipeline to work and have a explicit default but can be 'tweaked' (e.g. numeric values, such as modifying mapping sensitivity).  These are inserted into the workflow itself in the `<tool>_options.args +=` declarations that made prior to the `include` step of loading modules into the workflow (see [here](https://github.com/nf-core/viralrecon/blob/2ebae61442598302c64916bd5127cf23c8ab5611/workflows/illumina.nf#L50-L60) for an example in viralrecon), and can be evaluated for correctness by the pipeline prior passing to the module. This declaration _extends_ whatever tool command-line defaults defined in the `modules` block.

The `params` block of a config is much like in DSL1, and is used to modify top-level parameters that can be modified using the equivalent command-line flag.

In contrast, the modules block of a config-file is for more complex/in depth (and somewhat more risky) modifications to pipeline tool commands. Parameters in the `modules` block are not modifiable by a user on the command line. This is genereally used by a pipeline developer to hard code arguments of a given module that will always be run by the pipeline. However, users can provide a custom config with a modules block and therefore with a custom `args` parameter. With a this a user directly modifies the command-line arguments of a given tool, and therefore can include parameters that are _not_ evaluated by the pipeline itself for correctness. Importantly, by making these modifications a user can potentially 'break' the pipeline as they may add unsupported output that the pipeline does not export or send downstream. However, this does give a user a lot more freedom to do precisely what they want, if not (yet) directly supported by the pipeline.
