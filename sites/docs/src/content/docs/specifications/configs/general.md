---
title: General
subtitle: Follow general config guidelines
markdownPlugin: addNumbersToHeadings
weight: 1
---

The keywords "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Permissions to host on nf-core/configs

Configs hosted on nf-core/configs SHOULD have permission from the administrators of the given infrastructure to host the config publicly within the nf-core GitHub organisation, except for the exceptions described below.

- Configs for sensitive data clusters MUST have permission from the system administrators.

:::tip
A config MAY declare the config is 'unofficial' if the system administrators are OK with hosting publicly but do not maintain it themselves.
:::

## Alignment of config with local policies

Configurations SHOULD comply with and document administrative policies of the infrastructure as far as possible.
For example, if multiple possible partitions exist but there is a policy to use specific partitions for specific cases, this SHOULD be represented in the config.
Another example is if executing the main Nextflow run command on login/submit nodes is not allowed, include a sample job submission script (for example a SLURM `sbatch` script) in the documentation.

## Sensitive data and offline environments

### Sensitive data clusters

Configs for infrastructure handling sensitive or restricted data SHOULD document
any relevant data governance policies that affect pipeline execution.
A config for sensitive data infrastructure SHOULD use a local container registry
rather than pulling from public registries (e.g., Docker Hub, Quay.io).

### Offline/air-gapped clusters

A config for infrastructure without external internet access SHOULD describe how to set
`singularity.cacheDir` or `apptainer.cacheDir` to a directory pre-populated
using [`nf-core download`](https://nf-co.re/docs/nf-core-tools/cli/pipelines/download) rather than configuring a container registry.
It SHOULD set `params.igenomes_ignore = true` and provide paths to locally
available reference genomes instead.

## Size of configs

### Number of infrastructures in a config

A single configuration file SHOULD only be used to represent a single cluster or type of infrastructure.

For HPC infrastructure, a single config MAY represent multiple similar or linked-HPCs that are dynamically selected within the config.

### Structure when multiple configs

If multiple HPCs are supported in a single config, any sub-configs that are selected based on a condition in the main config MUST be placed in a subdirectory.
The subdirectory MUST have the same base name as the config file (i.e., for a config called `myinstitute.config`, a directory named `myinstitute`) and be loaded in the main config with `includeConfig()`.

Example:

```tree
conf/
├── <myinstitute>/
│   ├── <hpc1>.config
│   └── <hpc2>.config
└── <myinstitute>.config
```

## Scope of configs

Two config types are possible:

- **[Institutional config](#institutional-config)** that directs how any pipeline works on a given infrastructure.
- **[Pipeline specific config](#pipeline-specific-institutional-config)** that optimises how a specific pipeline works on a given infrastructure.

### Institutional config

An institutional config MUST be compatible with any nf-core pipeline or user.
It defines how the pipeline interacts with the infrastructure, such as scheduling options, software environment settings, and resource limits.

An institutional config SHOULD NOT define any resource defaults with `withName` or `withLabel`.
It SHOULD provide reasonable default settings for operating on the infrastructure (e.g. `resourceLimits`, `beforeScript`, `clusterOptions`, `runOptions`).

### Pipeline specific institutional config

A pipeline specific config MAY modify the default resource values (memory, CPUs, time).

Where possible, it's RECOMMENDED to provide defaults for using locally available references and similar resources.

## Naming

### Size of name

Institutional configs SHOULD use a short name or acronym as the config name.

### Formatting

Config names MUST be written (e.g. file names) or referred to (in documentation) in all lowercase letters or numbers.

Configs names MAY use an underscore.
Config names MUST NOT use any other symbols.

### Names when multiple infrastructure in a single institution

When multiple computational infrastructure exist for a single institution, an institutional prefix SHOULD be used.

For example, for the Max Planck Computing and Data Facilitation (MPCDF) institution that has two HPCs named `raven` and `viper`:

```groovy
mpcdf_raven.conf
mpcdf_viper.conf
```

## Required files

An institutional config MUST consist of two files:

- **`conf/<config name>.config`** for the config itself.
- **`docs/<config name>.md`** for documentation about the config.

Files for a pipeline-specific institutional config must be located:

- **[`conf/pipeline/<pipeline name>/<config name>.config`](https://github.com/nf-core/configs/tree/master/conf/pipeline)**
- **[`docs/pipeline/<pipeline name>/<config name>.md`](https://github.com/nf-core/configs/tree/master/docs/pipeline)**

Furthermore, the config MUST be referred to in two additional places

- **[nf-core_custom.config](https://github.com/nf-core/configs/blob/master/nfcore_custom.config)**: so the configs can be loaded by any nf-core pipeline
- **[.github/workflows/main.yml](https://github.com/nf-core/configs/blob/master/.github/workflows/main.yml)**: so the config is tested and linted using the GitHub action CI tests

## Required information

A config MUST have a current contact person responsible for maintaining the config.

## Parameters

## Required parameters

A config MUST include three descriptive parameters:

| Parameter                    | Purpose                                                                    |
| ---------------------------- | -------------------------------------------------------------------------- |
| `config_profile_description` | A short description which infrastructure the config is used for.           |
| `config_profile_contact`     | The name and GitHub handle of the person currently maintaining the config. |
| `config_profile_url`         | A URL to details about the infrastructure or the institution.              |

```groovy
params {
    config_profile_description = 'The <name of infrastructure> cluster profile'
    config_profile_contact     = '<maintainer name> (@<github handle>)'
    config_profile_url         = 'https://<url>.com'
}
```

## Optional parameters

### Parameters

A config MAY also define the `max_*` parameters with the same values as the `resourceLimit` directive.
This provides backwards compatibility of older pipelines with older Nextflow versions.

```groovy
params {
    igenomes_ignore = true
    max_memory      = 750.GB
    max_cpus        = 200
    max_time        = 30.d
}
```

### Custom parameters

Custom config- or infrastructure-specific parameters MAY be used, such as for cluster scheduler 'account' or 'project' parameters.

Custom config- or infrastructure-specific parameters MUST be documented in the configs `.md`.
Custom config- or infrastructure-specific parameters MUST be included in an nf-schema validation scope `ignoreParams` parameter.

For example:

```groovy
validation {
    ignoreParams = ['cluster_account']
}
```

## Resource limits

### Directive

A config MUST define the maximum resource limits of a computing infrastructure using the [`resourceLimits`](https://www.nextflow.io/docs/latest/reference/process.html#resourcelimits) process directive.

```groovy
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
    executor = 'slurm'
    queue    = { task.memory <= 250.GB ? (task.time <= 24.h ? 'fast' : 'long') : 'bigmem' }
}
```
