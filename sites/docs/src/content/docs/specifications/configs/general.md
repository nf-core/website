---
title: General
subtitle: Follow general config guidelines
markdownPlugin: addNumbersToHeadings
weight: 1
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Permissions

Configs hosted on nf-core/configs SHOULD have permission from the administrators of the given infrastructure to host the config publicly within the nf-core GitHub organisation.

:::tip
A config MAY declare the config is 'unofficial' if the system administrators are OK with hosting publicly, as long as they do not maintain it themselves.
:::

## Size of configs

A single config SHOULD represent a single cluster or type of infrastructure.

A single config MAY include multiple configs using sub `profiles`, when <!-- TODO -->

## Scope of configs

There are two config types:

- **[Institutional config](#institutional-config)** that directs how any pipeline works on a given infrastructure
- **[Pipeline specific config](#pipeline-specific-institutional-config)** that optimises how a specfic pipeline works on a given infrastructure

### Institutional config

An institutional config MUST be compatible with any nf-core pipeline or user.
It defines how the pipeline interacts with the infrastructure, such as scheduling options, software environment settings, and resource limits.

An institutional config SHOULD NOT define any resource defaults with `withName` or `withLabel`. It SHOULD provide reasonable default settings for operating on the infrastructure (e.g. `resourceLimits`, `beforeScript`, `clusterOptions`, `runOptions`).

### Pipeline specific institutional config

A pipeline specific config MAY modify the default resource values (memory, CPUs, time).

Where possible, it's RECOMMENDED to provide defaults for using locally available references and similar resources.
## Naming

Institutional configs SHOULD use a short name or acronym as the config name.

The name MUST be used

When multiple computational infrastructures exist for a single institution, an institutional prefix MAY be used.

Example, for the MPCDF institution that has two HPCs named `raven` and `viper`:

```groovy
mpcdf_raven.conf
mpcdf_viper.conf
```

## Required files

An institutional config MUST consist of two files:

- `conf/<config name>.conf`
- `docs/<config name>.md`

Files for a pipeline-specific institutional config must be located:

- **[`conf/pipeline/<pipeline name>/<config name>.conf`](https://github.com/nf-core/configs/tree/master/conf/pipeline)**
- **[`docs/pipeline/<pipeline name>/<config name>.md`](https://github.com/nf-core/configs/tree/master/conf/pipeline)**

Furthermore, the config MUST be referred to in two additional places

- **[nf-core_custom.config](https://github.com/nf-core/configs/blob/master/nfcore_custom.config)**: so the configs can be loaded by any nf-core pipeline
- **[.github/workflows/main.yml](https://github.com/nf-core/configs/blob/master/.github/workflows/main.yml)**: so the config is tested and linted using the GitHub action CI tests

## Required information

A config MUST have a current contact person responsible for maintaining the config.

## Required parameters

A config MUST include three descriptive parameters:

| Parameter                    | Purpose                                                                  |
| ---------------------------- | ------------------------------------------------------------------------ |
| `config_profile_description` | A short description which infrastructure the config is used for.         |
| `config_profile_contact`     | The name and GitHub handle of the person currently maintaing the config. |
| `config_profile_url`         | A URL to details about the infrastructure or the institution.            |

```groovy
params {
    config_profile_description = 'The <name of infrastructure> cluster profile'
    config_profile_contact     = '<maintainer name> (@<github handle>)'
    config_profile_url         = 'https://<url>.com'
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
