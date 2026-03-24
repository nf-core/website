---
title: Formatting
subtitle: Configuration file formatting requirements
markdownPlugin: addNumbersToHeadings
weight: 9
---

The keywords "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Strict syntax

nf-core/configs MUST be [Nextflow strict syntax](https://www.nextflow.io/docs/latest/strict-syntax.html) compliant to ensure compatibility with all nf-core pipelines.

Compliance of a config to the strict syntax can be checked using the methods described in the strict syntax migration [guide](../../developing/migration-guides/config-strict-syntax).

### Variables

Variables MUST be defined as a parameter (`param`) for reusing information to different parts of a config.

```groovy
- def variable_name = <code>
+ params.variable_name = <code>
```

Variables that are defined using a condition MUST be wrapped in a closure and executed with the `.call()` function.

```groovy
params.variable_name = {
    if (large_data == true) {
        return 'small'
    } else {
        return 'big'
    }
}.call()
```

### Functions

Functions MUST NOT be used within a config.

### Conditional statements

Simple if-else statements SHOULD be written using a ternary operator.

```groovy
process.executor = params.slurm ? 'slurm' : 'local'
```

Larger conditions SHOULD be written using if-else statements.
Switch conditions MUST NOT be used.

```groovy
queue = {
  if (task.memory >= 216.GB) {
      if (task.time >= 7.d) {
          return 'longmem'
      } else {
          return 'mem'
      }
  } else {
      if (task.time >= 21.d) {
          return 'long60'
      } else if (task.time >= 7.d) {
          return 'long'
      } else if (task.time >= 48.h) {
          return 'medium'
      } else {
          return 'short'
      }
  }
}
```

### Environmental Variables

Environmental variables MUST be referenced with a `System.getenv` call.

```groovy
scratch      = "/scratch/${System.getenv('USER')}"
```
