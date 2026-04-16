---
title: Migrating to strict syntax
subtitle: Update nf-core/configs for Nextflow strict syntax compliance
description: How to update Nextflow configs for strict syntax compliance
shortTitle: Migrating to strict syntax
---

From version 26.04, Nextflow enables [strict syntax](https://www.nextflow.io/docs/latest/strict-syntax.html#preparing-for-strict-syntax) mode by default, a more restrictive form of Nextflow that improves error messages and makes code more consistent.

This guide describes how to update your config for strict syntax compliance, for example in [nf-core/configs](https://github.com/nf-core/configs/).

:::note{title="Prerequisites"}

You will need the following to get started:

- A local copy of the `nf-core/configs` repository (ideally a fork on a separate branch if you need to make changes)
- VS Code with the [Nextflow extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow), or Nextflow version `26.02.0-edge` or later for command-line testing

:::

## Check for strict syntax compliance

You can check compliance in two ways:

- Using the Nextflow VS Code extension
- Using Nextflow on the command line

### VS Code

VS Code with the [Nextflow extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow) is the quickest way to check for problems with your config.

To check your config with the Nextflow VS Code extension:

1. Open the `nf-core/configs` repo as a VS Code project.

   The extension language server highlights any issues using [hover hints or the diagnostics window](https://www.nextflow.io/docs/latest/vscode.html#diagnostics).

### Command line

If you don't use VS Code, you can use Nextflow on the command line.

:::warning
This method checks syntax only.
See [Final testing](#final-testing) to verify functionality.
:::

To test whether your config is strict syntax compliant, create an empty workflow and run it using your config as a custom config file.

1. In your local copy or fork of `nf-core/configs`, create a file called `main.nf`:

   ```bash
   echo "workflow{}" > main.nf
   ```

1. Run the minimal workflow with the latest edge version of Nextflow and point to the config file you want to test:

   ```bash
   NXF_VER=26.02.0-edge nextflow run main.nf -c conf/<your_config>.config
   ```

If you see no warning or error messages, your config is already compliant.
If you see messages, see the next section for common fixes.

## Common fixes

This section describes common fixes and solutions for issues in existing `nf-core/configs`.

### Simple variables

The strict syntax does not allow variables in configs, so you will see an error like this:

```console
Error <config_name>.config:1:1: Variable declarations cannot be mixed with config statements
```

If your config uses 'simple' variables, for example static values like `scratch_dir = '/ptmp'`, follow the steps below.
For variables with conditions inside, see [Dynamic variables](#dynamic-variables).

First, check whether the variable is used multiple times.
If it isn't, remove the variable entirely and use its contents directly where the variable was being used.
If the variable is used multiple times, convert it to a parameter:

```diff groovy
- def variable_name = <code>
+ params.variable_name = <code>
```

:::warning
Make sure the parameter names are unique and isolated to the config so they don't overwrite anything in the pipelines themselves.
We recommend the format `<config_name>_<variable_name>`, but you can make it more unique if needed.
:::

To prevent [nf-schema](https://github.com/nextflow-io/nf-schema) warnings during pipeline initialisation, add custom parameters to the `ignoreParams` list in your config:

```groovy
validation {
    ignoreParams = [
        <list_all_used_parameters>
    ]
}
```

For a real-world example, see [nf-core/configs PR #1013](https://github.com/nf-core/configs/pull/1013/changes#diff-4adf4ac5644b18a57f91cb4b39a6d52d8f6618253441c4a61f4efa9a42a25956R107).

### Dynamic variables

The strict syntax does not allow variables in configs, so you will see an error like this:

```console
Error <config_name>.config:1:1: Variable declarations cannot be mixed with config statements
```

If your config uses dynamic variables, for example variables with a condition inside, follow the steps below.
For static variables, see [Simple variables](#simple-variables).

Convert these variables to parameters.
For simple cases, the code can often be converted to a ternary one-liner:

```groovy
params.random_var = something ? 'not_random' : 'random'
```

For more complicated cases, embed the conditions inside a parameter containing a closure and use `.call()` to evaluate it on config resolution.
For example, the following assignment:

```groovy
def random_var = ''
if (something == true) {
    random_var = 'not_random'
} else {
    random_var = 'random'
}
```

Becomes:

```groovy
params.random_var = {
    if (something == true) {
        return 'not_random'
    } else {
        return 'random'
    }
}.call()
```

:::warning
Don't forget to add `.call()` at the end.
This ensures the code is evaluated during config resolution.
Without `.call()`, the parameter will be a closure instead of the expected value.
:::

To prevent [nf-schema](https://github.com/nextflow-io/nf-schema) warnings during pipeline initialisation, add custom parameters to the `ignoreParams` list in your config:

```groovy
validation {
    ignoreParams = [
        <list_all_used_parameters>
    ]
}
```

For a real-world example, see [nf-core/configs PR #1013](https://github.com/nf-core/configs/pull/1013/changes#diff-4adf4ac5644b18a57f91cb4b39a6d52d8f6618253441c4a61f4efa9a42a25956R95-R101).

### Functions

The strict syntax no longer allows functions in configs, so you will see an error like this:

```console
Error <config name>.config:630:14: Unexpected input: '('
│ 630 | def check_max(obj, type) {
```

Refactor your code to avoid using functions.

:::warning
There is currently no clean solution for this case because of a conflict with nf-schema, which is used for pipeline input validation.
For now, reimplement the code at each use of the function within a closure.
:::

~~All functions should be converted to callable closures that are assigned to a parameter.
For example, the following function:~~

```groovy
def calculate_something(memory, time) {
    def output = null
    // function code
    return output
}
```

~Becomes:~

```groovy
params.calculate_something = { memory, time ->
    def output = null
    // function code
    return output
}
```

~Calling the function can then be done via `params.calculate_something(memory, time)` instead of `calculate_something(memory, time)`.~

~For a real-world example, see [nf-core/configs PR #1015](https://github.com/nf-core/configs/pull/1015/changes#diff-c60bd9c6097498d07b2f2eb3937b7d4ab3cb15e9167bacf80cb49c9848806e6fR117-R119).~

### Basic if-statements

The strict syntax no longer allows full if-else statements, so you will see an error like this:

```console
Error <config>.config:48:5: If statements cannot be mixed with config statements
```

Shorten basic if statements (if-statements with one line per condition) using the Groovy `?:` ternary syntax.
For example, the following if-else block:

```groovy
if(params.slurm) {
    process.executor = 'slurm'
} else {
    process.executor = 'local'
}
```

Becomes:

```groovy
process.executor = params.slurm ? 'slurm' : 'local'
```

For long lines, spread this over multiple lines for readability:

```groovy
process.executor = params.slurm ?
    'slurm' :
    'local'
```

:::note
You can also use a closure-wrapped approach:

```groovy
process.executor = {
    if (params.slurm) {
        return 'slurm'
    }
    return 'local'
}.call()
```

This is similar to a ternary operator but less readable.
This approach is discouraged for simple if-statements, but it can be useful for more complex conditions.
It's up to the config developer to decide when to use each method.
:::

For a real-world example, see [nf-core/configs PR #1013](https://github.com/nf-core/configs/pull/1013/changes#diff-4adf4ac5644b18a57f91cb4b39a6d52d8f6618253441c4a61f4efa9a42a25956R119-R157).

### Environment variables

The strict syntax no longer allows calling execution (shell) environment variables directly, so you will see an error like this:

```console
Error <config name>.config:21:32: `USER` is not defined (hint: use `env('...')` to access environment variable)
```

Wrap any environment variables in the Groovy `System.getenv()` function.
For example:

```groovy
scratch      = "/scratch/${USER}"
```

Becomes:

```groovy
scratch      = "/scratch/${System.getenv('USER')}"
```

For a real-world example, see [nf-core/configs PR #1019](https://github.com/nf-core/configs/pull/1019/changes#diff-3a625144f9f4ce62e80059132891092abbf804c666dd1dd983ad1068e91f46d0R21).

### Switch statements

The strict syntax no longer allows switch statements, so you will see an error like this:

```console
Error <config name>.config:27:64: Unexpected input: ':'
```

Change the switch statement to a closure-wrapped if-else statement.
For example, the following switch block:

```groovy
queue = {
  switch (task.memory) {
      case { it >= 216.GB }:
          switch (task.time) {
              case { it >= 7.d }:
                  return 'longmem'
              default:
                  return 'mem'
          }
      default:
          switch (task.time) {
              case { it >= 21.d }:
                  return 'long60'
              case { it >= 7.d }:
                  return 'long'
              case { it >= 48.h }:
                  return 'medium'
              default:
                  return 'short'
          }
  }
}
```

Becomes:

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

For a real-world example, see [nf-core/configs PR #1025](https://github.com/nf-core/configs/pull/1025/changes#diff-ddd72bdeaafae60efe36de07b84d991665e791811765d7987a7bf90cc1cd8584L24).

## Complete examples

For examples of entire configs that were made strict syntax compliant, see:

- [nf-core/configs PR #1013](https://github.com/nf-core/configs/pull/1013)
- [nf-core/configs PR #1015](https://github.com/nf-core/configs/pull/1015)

## Final testing

Once you have made your changes and tested them with the [basic workflow](#command-line), we recommend testing the updated config on a real nf-core pipeline.

If you are updating the config on a fork or branch, you can use these two parameters:

- `--custom_config_version`: Specify a different branch (for example, a branch within `nf-core/configs`)
- `--custom_config_base`: Specify a different fork

For example, for a fork:

```bash
nextflow pull nf-core/demo
NXF_VER=26.02.0-edge nextflow run nf-core/demo -profile test,<config name> --custom_config_base 'https://github.com/<your user name>/nf-core-configs/raw/refs/heads/<your branch name>/nfcore_custom.config'
```

For a branch on `nf-core/configs`:

```bash
nextflow pull nf-core/demo
NXF_VER=26.02.0-edge nextflow run nf-core/demo -profile test,<config name> --custom_config_version '<your-fixes-branch>'
```
