---
title: Migrating modules to static types and records
subtitle: Migrate modules to static types and records
description: How to migrate nf-core modules to Nextflow static types and records
shortTitle: Migrating to static types
---

[Static types](https://docs.seqera.io/nextflow/process-typed) let a module declare one typed input record and one typed output record, instead of several tuple channels.
Version reporting, `ext` variables, and staged scripts also change format.

This guide describes how to migrate an nf-core module to static types and records.
For the rules that a migrated module must follow, see the [static types and records specifications](/docs/specifications/components/modules/static-types).

:::note{title="Prerequisites"}

You will need the following to get started:

- Nextflow 26.10.0 or later
- A clone of the `nf-core/modules` repository
- nf-core tools
- nf-test

:::

## Migrate the module `main.nf`

To migrate the module `main.nf` file:

1. Add `nextflow.enable.types = true` as the first line of the file.

1. Add an `ext` directive that defines the prefix and arguments as closures:

   ```groovy title="main.nf"
   ext [
       prefix: { output_prefix ?: "${id}" },
       args: { <tool>_args ?: '' }
   ]
   ```

   Replace `<tool>` with the tool name, for example `samtools_stats`.

1. Replace the input tuples with a single input record:

   ```groovy title="main.nf"
   input:
   record(
       id: String,
       meta: Record?,
       <input_file>: Path,
       <tool>_args: String?,
       output_prefix: String?
   )
   ```

   Add a field for each required and optional input file of the process.
   Mark optional fields with `?` after the type.

1. Remove the `def args` and `def prefix` lines from the `script:` and `stub:` blocks.
   Use `${task.ext.args}` and `${task.ext.prefix}` instead.

1. Remove the `when:` block:

   ```groovy title="main.nf"
   when:
   task.ext.when == null || task.ext.when
   ```

1. Replace the output tuples with a single output record and remove the `emit:` names:

   ```groovy title="main.nf"
   output:
   record(
       id: id,
       meta: meta,
       bed: file("*.bed")
   )
   ```

   For optional outputs, use `file("*.bed", optional: true)`.

1. Replace the version output with a `topic:` block:

   ```groovy title="main.nf"
   topic:
   record(process: "${task.process}", tool: '<tool>', version: eval('<version-command>')) >> 'versions'
   ```

   Replace `<tool>` with the tool name and `<version-command>` with the command that prints its version.
   If the module writes a `versions.yml` file, remove the `versions.yml` HEREDOC from the `script:` and `stub:` blocks.

## Replace a template with a staged script

Modules with static types stage a script file and run it directly, instead of using the `template` function.

To replace a template with a staged script:

1. Add a `stage:` block that stages the script into the task directory:

   ```groovy title="main.nf"
   stage:
   stageAs file("${moduleDir}/templates/<script>"), '*'
   ```

1. Set the executable bit on the script:

   ```bash
   chmod +x templates/<script>
   ```

1. Run the script in the `script:` block:

   ```groovy title="main.nf"
   script:
   """
   ./<script> ${task.ext.args} ${input} ${task.ext.prefix}
   """
   ```

   Replace `input` with the name of the input file field.

## Convert the template file to a standalone script

Nextflow doesn't process a staged script as a template, so the script must read its values from command-line arguments.

To convert the template file:

1. Unescape the special characters that don't need escaping outside a Nextflow string.
   For example, replace `\\` with `\` and `\$` with `$`.

1. Replace each Nextflow variable, such as `${task.ext.prefix}`, with a command-line argument.
   Parse options from `${task.ext.args}` with the script's argument parser.

1. Remove the version reporting block from the script.
   The `topic:` block in `main.nf` now reports the version.

1. Check the script syntax.
   For example, for a Perl script:

   ```bash
   perl -c templates/<script>
   ```

## Update the tests

Tests must pass inputs as records.
Nextflow rejects a plain map as input to a module with static types.

To update the tests:

1. In `tests/main.nf.test`, pass each input as a record:

   ```groovy title="main.nf.test"
   record(id: '<sample-id>', <input_file>: file('<path-to-test-file>'), <tool>_args: null, output_prefix: null)
   ```

1. Move any per-test `ext.args` from `tests/nextflow.config` to the `<tool>_args` field of the record.

1. Run the tests and update the snapshot:

   ```bash
   nf-test test <path> --profile=+docker --update-snapshot
   ```

   :::warning
   `--update-snapshot` overwrites the `.snap` file even when the test fails.
   If a test run fails, restore the snapshot with `git checkout <path>/tests/main.nf.test.snap` before you run the tests again.
   :::

1. Check that the output md5 checksums in the new snapshot match the old snapshot.

## Update documentation and callers

To finish the migration:

1. Update the module `meta.yml` file to match the new input and output records.

1. Lint the module:

   ```bash
   nf-core modules lint <module-name>
   ```

1. Update the subworkflows and pipelines that call the module, and update their snapshots.
