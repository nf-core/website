---
title: ext properties/keys
description: Learn about the ext properties or keys in nf-core components
shortTitle: <code>ext.args</code>
parentWeight: 30
---

Ext properties or keys are special process directives (See: [ext directive](https://www.nextflow.io/docs/latest/process.html#ext) ) that insert strings into the module scripts. For example, an nf-core module uses the string assigned to `ext.args` ( or `ext.args2`, `ext.args3`, ... ) to insert tool specific options in a module script:

Example:

The configuration

```groovy
process {
  withName: 'TOOL_SUBTOOL' {
    ext.args = '-T -K'
  }
}
```

inserts the string `-T -K` as options to the module script:

```groovy
process TOOL_SUBTOOL {
  input:
  tuple val(meta), path(bam)

  output:
  tuple val(meta), path("*.log"), emit: log
  path "versions.yml",            emit: versions

  script:
  def args   = task.ext.args ?: ''          // If ext.args is defined assign it to args
  def prefix = task.ext.prefix ?: meta.id   // If ext.prefix is defined assign it to prefix, otherwise assign meta.id value
  """
  tool subtool $args $bam > ${prefix}.log
  """
}
```

so the script becomes:

```bash
#! /usr/env/bin bash
tool subtool -T -K test.bam > test.log
```

## Permitted ext keys

The following table lists the permitted keys used in nf-core modules. Other keys must be discussed with
maintainers, and added here when permitted.

| Key         | Description                                            |
| ----------- | ------------------------------------------------------ |
| ext.args    | Additional arguments appended to command in module.    |
| ext.args2   | Second set of arguments appended to command in module. |
| ext.args3   | Third set of arguments appended to command in module.  |
| ext.argsN   | Nth set of arguments appended to command in module.    |
| ext.prefix  | File name prefix for output files.                     |
| ext.when    | Boolean expression to determine when a module runs.    |
| ext.use_gpu | Determines whether the module uses GPU settings.       |
| ext.singularity_pull_docker_container | Whether to use the Docker URI instead of the Singularity URI under the `singularity` profile |

:::note
The order of the numeric ID of `args` must match the order of the tools as used in the module.
:::

To see some more advanced examples of these keys in use see:

- [Set ext.args based on parameter settings](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/modules.config#L222-L226)
- [Set ext.prefix based on task inputs](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/modules.config#L297)
- [Set ext.args based on both parameters and task inputs](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/modules.config#L377-L381)

## Rational

The ext keys are properties that are intended to allow a pipeline user to change via configuration file how a module behaves. They are not a means for a developer to bypass passing values via a module's `input:`. The majority of optional command-line flags are passed through `ext.args`. If changing an ext key would lead to pipeline instability, it should be an `input:`. This allows module parameters to be fully described and documented.
