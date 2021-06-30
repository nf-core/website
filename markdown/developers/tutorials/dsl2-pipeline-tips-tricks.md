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

If you have errors such as `input not found` or similar during testing, you can debug whether the defined nested param is correct by checking the corresponding config.
This `params.test_data` is defined in the file [`modules/tests/config/test_data.config`](https://github.com/nf-core/modules/blob/master/tests/config/test_data.config).

## Pipeline

### Converting existing DSL1 pipeline to DSL2

When converting an existing nf-core pipeline from DSL1 to DSL2, it's important that you don't remove the previous `master` that contains the current DSL1 release, nor work on your DSL2 port on `dev`, as you may need to make patch releases to the DSL1 version in the meantime.

Instead, create a new branch called e.g. `dsl2`. Pull this branch to your local machine, and then while in that branch delete the entire contents of the directory to start from scratch.

Once deleted, run `nf-core create` (ensuring you're using at least version 1.15, which includes the new DSL2 template), filling in the required information as with the original pipeline.

Finally, compare the `tree` for both DSL1 and DSL2 branches, and do `diff` on files individually to add in custom code from the pipeline. If in doubt, you can always check for any missing files or content by looking directly at the nf-core pipeline template
