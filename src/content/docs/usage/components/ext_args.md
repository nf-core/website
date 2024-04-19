---
title: ext properties/keys
description: Learn about the ext properties or keys in nf-core components
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

The following table lists the available keys commonly used in nf-core modules.

| Key        | Description                                            |
| ---------- | ------------------------------------------------------ |
| ext.args   | Additional arguments appended to command in module.    |
| ext.args2  | Second set of arguments appended to command in module. |
| ext.args3  | Third set of arguments appended to command in module.  |
| ext.prefix | File name prefix for output files.                     |

:::note
The order of the numeric ID of `args` must match the order of the tools as used in the module.
:::

To see some more advanced examples of these keys in use see:

- [Set ext.args based on parameter settings](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/modules.config#L222-L226)
- [Set ext.prefix based on task inputs](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/modules.config#L297)
- [Set ext.args based on both parameters and task inputs](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/modules.config#L377-L381)
