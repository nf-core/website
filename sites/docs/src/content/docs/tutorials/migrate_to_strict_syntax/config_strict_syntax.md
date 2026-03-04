---
title: Making nf-core/configs strict syntax compliant
subtitle: Announcement and recommendations for updating nf-core/configs Nextflow strict syntax compliant
---

## Introduction

From 26.04, Nextflow will make its new '[Strict Syntax](https://www.nextflow.io/docs/latest/strict-syntax.html#preparing-for-strict-syntax)' mode default.
It is a more restrictive way of writing Nextflow that improves error messages and makes the code more consistent.

This page describes ways to make your config syntax compliant. For example, for [nf-core/configs](https://github.com/nf-core/configs/).

## Check for strict syntax compliance

There are two ways to check if your existing config is compliant with the new strict syntax.

In both cases, you need a local copy of the nf-core/configs repository (ideally as a fork on a separate branch if you need to make changes).

### VS Code

VS Code with the [Nextflow extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow) is the easiest way to check for problems with your config.

To check your config with the Nextflow VS Code extension:

1. Open the nf-core/configs repo as a VS Code project
   - The extension language server will highlight any issues using [hover hints or the diagnostics window](https://www.nextflow.io/docs/latest/vscode.html#diagnostics).

### Command line

If you don't use VSCode, you can instead use Nextflow itself on the command line.

> [!WARNING]
> This does not check the config functionality works!
> See end of instructions

The fastest way to test whether your config is Nextflow strict syntax compliant is to make an empty workflow, and run it using your config as a custom config file.

For example:

1. In your local copy of fork of nf-core/configs, make a file called `main.nf`

   ```bash
   echo "workflow{}" > main.nf
   ```

2. Run the minimal workflow with the latest edge version of Nextflow (`26.02.0-edge`) and pointing to the config file you want to test

   ```bash
   NXF_VER=26.02.0-edge nextflow run main.nf -c conf/<your_config>.config
   ```

If you have no warning or errors messages, you are good to go - your config is already compliant!

If you have messages, see the next section on how to resolve the most common errors.

## Common fixes

This section describes common fixes and solutions seen across existing nf-core/configs.

### Simple variables

#### Problem

Variables are not allowed in the configs under the strict syntax.

If you get an error like:

```groovy
Error <config_name>.config:1:1: Variable declarations cannot be mixed with config statements
```

If your config use variables that are 'simple', i.e. static values `scratch_dir = ''/ptmp'`, continue reading this section.
If you have more complicated variables (e.g. variables with a condition inside), see [Dynamic variables](#dynamic-variables).

#### Solution

You will need to convert the variable to a parameter.

Check if the variable is actually used multiple times.

If it isn't, then remove the variable entirely, and just directly use the contents of the variable where the variable was being used.

If the variable is used multiple times, you can convert the variable to a parameter. For example:

```diff groovy
- def variable_name = <code>
+ params.variable_name = <code>
```

> [!WARNING]
> Make sure the parameter names are unique and isolated to the config so they don't overwrite anything in any pipelines themselves!
> We recommend: `<config_name>_<variable_name>`, but feel free to make it more unique.

To prevent [nf-schema](https://github.com/nextflow-io/nf-schema) warnings during pipeline initialisation, you should also add the following to your config:

```groovy
validation {
    ignoreParams = [
        <list_all_used_parameters>
    ]
}
```

#### Example

Real world example: [https://github.com/nf-core/configs/pull/1013/changes#diff-4adf4ac5644b18a57f91cb4b39a6d52d8f6618253441c4a61f4efa9a42a25956R107](https://github.com/nf-core/configs/pull/1013/changes#diff-4adf4ac5644b18a57f91cb4b39a6d52d8f6618253441c4a61f4efa9a42a25956R107)

### Dynamic variables

#### Problem

Variables are not allowed in the configs under the strict syntax.

If you get an error like:

```groovy
Error <config_name>.config:1:1: Variable declarations cannot be mixed with config statements
```

If your config uses variables that are dynamic (e.g. variables with a condition inside), continue reading this section.
If your config uses variables that are 'simple', i.e. static values `scratch_dir = ''/ptmp'`, see [Simple variables](#simple-variables).

#### Solution

You will need to convert these variables to a parameter.

Check if the code could not be converted to a ternary one-liner, e.g.

```groovy
params.random_var = something ? 'not_random' : 'random'
```

This is however only possible in very simple cases.

If you have more complicated cases, you will need to embed the conditions inside a parameter containing a closure and `call()` it.

The following example shows how to convert the assignment to a closure and call it on config resolution.

Example:

```groovy
def random_var = ''
if (something == true) {
    random_var = 'not_random'
} else {
    random_var = 'random'
}
```

Should be converted to the following:

```groovy
params.random_var = {
    if (something == true) {
        return 'not_random'
    } else {
        return 'random'
    }
}.call()
```

> [!WARNING]
> Don't forget to set the `.call()` here at the end.
> This makes sure the code is evaluated during config resolution.
> If you don't add `.call()`, the parameter will be a closure instead of the expected value.

To prevent [nf-schema](https://github.com/nextflow-io/nf-schema) warnings during pipeline initialisation, you should also add the following to your config:

```groovy
validation {
    ignoreParams = [
        <list_all_used_parameters>
    ]
}
```

#### Example

Real world example: [https://github.com/nf-core/configs/pull/1013/changes#diff-4adf4ac5644b18a57f91cb4b39a6d52d8f6618253441c4a61f4efa9a42a25956R95-R101](https://github.com/nf-core/configs/pull/1013/changes#diff-4adf4ac5644b18a57f91cb4b39a6d52d8f6618253441c4a61f4efa9a42a25956R95-R101)

### Functions

#### Problem

The use of functions in configs is no longer allowed by the strict syntax.

```groovy
Error <config name>.config:630:14: Unexpected input: '('
│ 630 | def check_max(obj, type) {
```

#### Solution

Try to refactor your code to not use any functions.

> [!WARNING]
> We do not currently have a good solution to this, due to a conflict with nf-schema used for pipeline input validation!
> Our only solution is to repeatedly implement the code at each use of the function within a closure.

~~All functions should be converted to callable closures that are assigned to a parameter.~~

~~Example:~~

```groovy
def calculate_something(memory, time) {
    def output = null
    // function code
    return output
}
```

~~Should become:~~

```groovy
params.calculate_something = { memory, time ->
    def output = null
    // function code
    return output
}
```

~~Calling the function can then be done via `params.calculate_something(memory, time)` instead of `calculate_something(memory, time)~~

#### Example

~~Real world example: [https://github.com/nf-core/configs/pull/1015/changes#diff-c60bd9c6097498d07b2f2eb3937b7d4ab3cb15e9167bacf80cb49c9848806e6fR117-R119](https://github.com/nf-core/configs/pull/1015/changes#diff-c60bd9c6097498d07b2f2eb3937b7d4ab3cb15e9167bacf80cb49c9848806e6fR117-R119)~~

### Basic if-statements

#### Problem

You can no longer use 'full' if-else statements with the strict syntax.

```groovy
Error <config>.config:48:5: If statements cannot be mixed with config statements
```

### Solution

Basic if statements (if-statements that usually have one line per condition) can be shortened using the `?:` groovy ternary syntax.

Example:

```groovy
if(params.slurm) {
    process.executor = 'slurm'
} else {
    process.executor = 'local'
}
```

Can become:

```groovy
process.executor = params.slurm ? 'slurm' : 'local'
```

This can also be spread out over multiple lines for readability of long lines:

```groovy
process.executor = params.slurm ?
    'slurm' :
    'local'
```

> [!NOTE]
> This could also be done using the following:
>
> ```groovy
> process.executor = {
>     if (params.slurm) {
>         return 'slurm'
>     }
>     return 'local'
> }.call()
> ```
>
> This is similar to using ternary operators but less readable.
> This is discouraged for simple if-statements, however you can use this in cases of more complex conditions.
> It's up to the config developer to decide when to use what method.

#### Example

Real world example: [https://github.com/nf-core/configs/pull/1013/changes#diff-4adf4ac5644b18a57f91cb4b39a6d52d8f6618253441c4a61f4efa9a42a25956R119-R157](https://github.com/nf-core/configs/pull/1013/changes#diff-4adf4ac5644b18a57f91cb4b39a6d52d8f6618253441c4a61f4efa9a42a25956R119-R157)

### Environmental variables

#### Problem

Calling execution (shell) environment variables directly are no longer allowed under the strict syntax.

```groovy
Error <config name>.config:21:32: `USER` is not defined (hint: use `env('...')` to access environment variable)
```

#### Solution

Wrap any environment variables in the `System.getenv('USER')` Groovy function.

For example:

```groovy
scratch      = "/scratch/${USER}"
```

Becomes

```groovy
scratch      = "/scratch/${System.getenv('USER')}"
```

#### Example

Real world example: [https://github.com/nf-core/configs/pull/1019/changes#diff-3a625144f9f4ce62e80059132891092abbf804c666dd1dd983ad1068e91f46d0R21](https://github.com/nf-core/configs/pull/1019/changes#diff-3a625144f9f4ce62e80059132891092abbf804c666dd1dd983ad1068e91f46d0R21)

### Switch statements

#### Problem

Switch statements are no longer allowed by the strict syntax.

```groovy
Error <config name>.config:27:64: Unexpected input: ':'
```

#### Solution

Change the switch to a closure-wrapped if-else statement.

For example:

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

Becomes

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

#### Example

Real world example: [https://github.com/nf-core/configs/pull/1025/changes#diff-ddd72bdeaafae60efe36de07b84d991665e791811765d7987a7bf90cc1cd8584L24](https://github.com/nf-core/configs/pull/1025/changes#diff-ddd72bdeaafae60efe36de07b84d991665e791811765d7987a7bf90cc1cd8584L24)

## Complete examples

Here is a list of example PRs of entire configs that were made strict syntax compliant.

- [https://github.com/nf-core/configs/pull/1013](https://github.com/nf-core/configs/pull/1013)
- [https://github.com/nf-core/configs/pull/1015](https://github.com/nf-core/configs/pull/1015)

## Final Testing

Once you have made your changes and you have tested with the [basic workflow](#testing-your-config-locally), it is HIGHLY recommended to test the updated config on a real nf-core pipeline.

If you are updating the config on a fork or branch, you can use the two parameters:

- `--custom_config_version`: to specify a different branch (e.g. on a branch within nf-core/configs)
- `--custom_config_base`: to specify a different fork

For example, for a fork:

```bash
nextflow pull nf-core/demo
NXF_VER=26.02.0-edge nextflow run nf-core/demo -profile test,<config name> --custom_config_base 'https://github.com/<your user name>/nf-core-configs/raw/refs/heads/<your branch name>/nfcore_custom.config'
```

and for a branch on nf-core configs:

```bash
nextflow pull nf-core/demo
NXF_VER=26.02.0-edge nextflow run nf-core/demo -profile test,<config name> --custom_config_version '<your-fixes-branch>'
```
