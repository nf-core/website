---
title: Configuration options
subtitle: How to configure nf-core pipelines
shortTitle: Configuration options
weight: 2
---

## Configuration options

nf-core pipelines may utilize any of these configuration options and files.

### Parameters

Parameters are pipeline specific settings that can be used to customize the execution of a pipeline.

At the highest level, parameters can be customized using the command line. Any parameter can be configured on the command line by prefixing the parameter name with a double dash (--):

```bash
--<parameter>
```

Depending on the parameter type, you may be required to add additional information after your parameter flag.
For example, you would add string parameter after the parameter flag for the `nf-core/rnaseq` `--input` and `--output` parameters.

```bash
nextflow nf-core/rnaseq --input <path/to/input> --outdir <path/to/results>
```

Every nf-core pipeline has a full list of parameters on the nf-core website. You will be shown a description and the type of the parameter when viewing these parameters. Some parameters will also have additional text to help you understand how a parameter should be used. See the [parameters page of the nf-core rnaseq pipeline](https://nf-co.re/rnaseq/3.14.0/parameters/).

### Default configuration files

All parameters have a default configuration that is defined using the `nextflow.config` file in the root of the pipeline directory. Many parameters are set to `null` or `false` by default and are only activated by a profile or config file.

nf-core pipelines also include additional config files from the `conf/` folder of a pipeline repository. Each additional `.config` file contains categorized configuration information for your pipeline execution, some of which can be optionally included as profiles:

- `base.config`
  - Included by the pipeline by default
  - Generous resource allocations using labels
  - Does not specify any method for software dependencies and expects software to be available (or specified elsewhere)
- `igenomes.config`
  - Included by the pipeline by default
  - Default configuration to access reference files stored on AWS iGenomes
- `modules.config`
  - Included by the pipeline by default
  - Module-specific configuration options (both mandatory and optional)
- `test.config`
  - Only included if specified as a profile
  - A configuration profile to test the pipeline with a small test dataset
- `test_full.config`
  - Only included if specified as a profile
  - A configuration profile to test the pipeline with a full-size test dataset

:::note
Some configuration files contain the definition of profiles that can be flexibly applied. For example, the `docker`, `singularity`, and `conda` profiles are defined in the `nextflow.config` file in the pipeline project directory. You should not need to manually edit any of these configuration files.
:::

Profiles are sets of configuration options that can be flexibly applied to a pipeline.
They are also commonly defined in the `nextflow.config` file in the root of the pipeline directory.

Profiles that come with nf-core pipelines can be broadly categorized into two groups:

- Software management profiles
  - Profiles for the management of software dependencies using container or environment management tools, for example, `docker`, `singularity`, and `conda`.
- Test profiles
  - Profiles to execute the pipeline with a standardized set of test data and parameters, for example, `test` and `test_full`.

nf-core pipelines are required to define software containers and environments that can be activated using profiles. Although it is possible to run the pipelines with software installed by other methods (e.g., environment modules or manual installation), using container technology is more sharable, convenient, and reproducible.

### Shared configuration files

nf-core pipelines can also load custom institutional profiles that have been submitted to the [nf-core config repository](https://github.com/nf-core/configs). At run time, nf-core pipelines will fetch these configuration profiles from the [nf-core config repository](https://github.com/nf-core/configs) and make them available.

For shared resources such as an HPC cluster, you may consider developing a shared institutional profile.

Follow [this tutorial](https://nf-co.re/docs/tutorials/config_institutional_profile) to set up your own institutional profile.

### Custom parameter and configuration files

Nextflow will look for files that are external to the pipeline project directory. These files include:

- The config file `$HOME/.nextflow/config`
- A config file named `nextflow.config` in your current directory
- Custom configuration files specified using the command line
  - A parameter file that is provided using the `-params-file` option
  - A config file that are provided using the `-c` option

**Parameter file format**

Parameter files are `.json` files that can contain an unlimited number of parameters:

```json title="nf-params.json"
{
  "<parameter1_name>": 1,
  "<parameter2_name>": "<string>",
  "<parameter3_name>": true
}
```

You can override default parameters by creating a `.json` file and passing it as a command-line argument using the `-param-file` option.

```bash
nextflow run nf-core/rnaseq -profile docker --input <path/to/input? --outdir <results> -param-file <path/to/nf-params.json>
```

**Configuration file format**

Configuration files are `.config` files that can contain various pipeline properties and can be passed to Nextflow using the `-c` option in your execution command:

```bash
nextflow run nf-core/rnaseq  -profile docker --input <path/to/input> --outdir <results> -c <path/to/custom.config>
```

Custom configuration files are the same format as the configuration file included in the pipeline directory.

Configuration properties are organized into [scopes](https://www.nextflow.io/docs/latest/config.html#config-scopes) by dot prefixing the property names with a scope identifier or grouping the properties in the same scope using the curly brackets notation. For example:

```groovy
alpha.x  = 1
alpha.y  = 'string value'
```

Is equivalent to:

```groovy
alpha {
    x = 1
    y = 'string value'
}
```

[Scopes](https://www.nextflow.io/docs/latest/config.html#config-scopes) allow you to quickly configure settings required to deploy a pipeline on different infrastructure using different software management.

A common scenario is for users to write a custom configuration file specific to running a pipeline on their infrastructure.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for tuning process resource specifications, other infrastructural tweaks (such as output directories), or module arguments (`args`).
:::

Multiple scopes can be included in the same `.config` file using a mix of dot prefixes and curly brackets.

```groovy
executor.name = "sge"

singularity {
    enabled    = true
    autoMounts = true
}
```

See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-scopes) for a full list of scopes.
