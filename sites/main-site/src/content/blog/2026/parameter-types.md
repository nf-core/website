---
title: "Why parameters are strings all of a sudden"
subtitle: "A guide to mitigating typing issues in strict syntax"
pubDate: 2026-03-05T09:00:00+01:00
headerImage: "assets/images/blog/parameter-types/keyboard.jpg"
headerImageAlt: "An RGB keyboard, because it's cool"
authors:
  - nvnieuwk
label:
  - pipelines
---

## Introduction

With the release of Nextflow 26.04.0, the default syntax parser has been set to v2. The new syntax parser will make it easier for Nextflow to throw useful errors whenever something breaks in your pipeline. Naturally with any change of this magnitude, stuff will start breaking. One of the most common issues reported in the last couple of weeks is that boolean parameter suddenly are only recognized as string values. In this blog post I will go over some ways to overcome this issue, from simple plugin updates to full parameter typing.

## Background

In syntax parser v1, Nextflow would automatically infer the type of all parameters given via the CLI. For example:

- `--answer_to_everything 42` will result in `params.answer_to_everything` to be an integer
- `--am_i_a_teapot false` will result in `params.am_i_a_teapot` to be a boolean
- `--hotel trivago` will result in `params.hotel` to be a string

In syntax parser v2, this system has been removed and all parameters given via the CLI are now always string values. This is causing some unexpected behaviour in pipelines (especially when using type validation in nf-schema for example). The following error is thus quite common after updating Nextflow to 26.04.0 or higher:

```bash
* --am_i_a_teapot (false): Value is [string] but should be [boolean]
```

So next, I will explain what to do as a pipeline user and what you can do as a pipeline developer to mitigate these issues.

## What to do as a pipeline user

As a user you have two options to mitigate the error when bumping Nextflow versions:

### Stop using CLI parameters

The best option as a pipeline user is to not pass any parameters to Nextflow via the CLI, but use a parameters file instead. A parameters file is a JSON or YAML file in which you can specify the parameters to be used by the pipeline. These files have support for simple types (strings, integers, booleans, floats...), thus making sure that Nextflow uses the expected types out of the box.

You can read more about the parameters file in the [Nextflow documentation](https://docs.seqera.io/nextflow/cli#pipeline-parameters)

### Bump nf-schema

nf-schema 2.7.2 and higher automatically converts CLI parameters to their correct types as it used to be in syntax parser v1. Keep in mind however that this conversion will only happen for the validation of the parameters and can still cause some unexpected issues in the pipeline. 

You can bump the plugin via a configuration file:

```groovy
plugins {
    id 'nf-schema@2.7.2'
}
```

or via the CLI with the following option:

```bash
-plugins nf-schema@2.7.2
```

## What to do as a pipeline developer

As a developer you can migrate your pipeline to start using [parameter types](https://docs.seqera.io/nextflow/workflow-typed#typed-parameters). Follow these steps to get an overview of how to do the migration:

1. Open your pipeline directory in VScode
1. Make sure the [Nextflow extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow) is installed
1. Open the `main.nf` file located in the root of the repository
1. Open the command options (CTRL+SHIFT+P or CMD+SHIFT+P for mac users), search for the `Nextflow: Convert script to static types` options and run it. This will create a `params` block with all types inferred from the `nextflow_schema.json` file
1. Check that all types are correct and that all defaults have been correctly filled in. Boolean values don't need a default as missing booleans always will be `false`. 
1. Optionally: Convert the type of all file parameters from `String` to `Path` to let Nextflow automatically convert these parameters to file objects. (This will probably need some tweaking in your pipeline code to remove unnecessary `file()` functions)
1. Make sure all parameters without a default that are optional have a `?` after the type. e.g. for an optional string parameter you would use the `String?` type if it has no default value. This will automatically assign `null` to that parameter.
1. Remove all parameters that are not used in configs or to define defaults for other parameters from the `nextflow.config` file. Read more about this in it's [chapter](#remove-parameters-from-nextflowconfig).

Ideally, the conversion is done now and your pipeline will be fully working again when users use syntax parser v2. There are however some caveats that you will need to account for to make sure everything works as expected. The following chapters explain these caveats and how to resolve them

### StackoverflowError

When using the `Path` type, you will most likely start seeing a `StackoverlflowError` when running your pipeline after the initial conversion. This is caused by a bug in older versions of the nf-schema plugin and `utils_nextflow_pipeline` subworkflow. You can fix this issue by bumping nf-schema to 2.7.1 or higher and by updating the `utils_nextflow_pipeline` and `utils_nfschema_plugin` subworkflows:

```bash
nf-core subworkflows update utils_nextflow_pipeline
nf-core subworkflows update utils_nfschema_plugin
```

You will see that `utils_nfschema_plugin` has a new input called `cli_typecast`. Set this to `false` to make sure nf-schema no longer typecasts parameters given via the CLI.

This should resolve the Stackoverflow errors!

### Remove parameters from `nextflow.config`

All parameter defaults should be removed from the `nextflow.config` file with a few exceptions listed below:

:::note
Don't specify defaults in `main.nf` for these parameters since these will never be used.
:::

#### Parameters that are used for configuration options

These parameters should still be defined in `nextflow.config` as parameter types are only applied after configuration resolution. Keep in mind that no typecasting has been done during config resolution to you should add `.toBoolean()`, `.toInteger()` or `.toFloat()` to all parameters that you do not expect to be string values.

This a non-exhaustive list of parameters that belong to this list. This depends a lot from pipeline to pipeline of course:

- `outdir`: used to define the `outputDir` option
- `publish_dir_mode`: used to define the `workflow.output.mode` and `publishDir mode` options
- `pipelines_testdata_base_path`: used in profiles to define test data
- `trace_report_suffix`: used to define the name of the pipeline reports
- `config_profile_name`, `config_profile_description`, `config_profile_contact` and `config_profile_url`: used by institutional configs
- `custom_config_version` and `custom_config_base`: used to initialize institutional configs
- `igenomes_ignore`: used to fetch the `genomes` block
- `monochrome_logs`: used for the `validation.monochromeLogs` option

:::note
Parameters used in `modules.config` should not have defaults in `nextflow.config` as these can be accessed during the pipeline run using closures (`{}`). 
:::

#### Parameters that are used to define defaults for other parameters

These parameters need to be set before the params block in `main.nf` is resolved, otherwise all values will be `null`. 

This a non-exhaustive list of parameters that belong to this list. This depends a lot from pipeline to pipeline of course:

- `igenomes_base`: used to set the base of the igenomes references
- `genome`: used to define from what organisms the igenome references should be used

### Suddenly I can't use igenomes references anymore

Support for nested parameters has been silently 'deprecated' with the introduction of parameter types. This issue can be resolved by migrating `genomes` parameter in `conf/igenomes.config` to a `Map` structure instead of using nested parameters. e.g.:

```groovy
params.genomes {
    'GRCh38' {
        fasta = "...",
        fai = "..."
    }
    'hg19' {
        fasta = "...",
        fai = "..."
    }
}
```

should become:

```groovy
params.genomes = [
    'GRCh38': [
        fasta: "...",
        fai: "..."
    ],
    'hg19': [
        fasta: "...",
        fai: "..."
    ]
]
```