---
title: Migrating to topic channels
description: Learn how to migrate nf-core modules and pipelines to use topic channels
shortTitle: Migrating to topic channels
parentWeight: 100
---

[Topic channels](https://www.nextflow.io/docs/latest/process-typed.html#topics) are a new feature in Nextflow that allow for more flexible and efficient handling of version outputs across modules and pipelines. Instead of collecting versions through YAML files, topic channels enable direct version tracking through structured channel outputs.

This migration guide provides step-by-step instructions for three different scenarios:

- **Updating modules**: For standard modules that generate version outputs in scripts
- **Updating modules with templates**: For modules that use template scripts to generate version information
- **Updating pipelines**: For pipelines that consume modules with topic channels

:::note{title="Prerequisites"}
You will need the following to get started:

- nf-core tools version 3.5.0 or later
- A clone of the `nf-core/modules` repository (for module updates)
- Nextflow version 25.04.0 or later
  :::

## Migrate modules

To migrate a module to use topic channels for version outputs:

1. Open the module `main.nf` file.

1. Identify the code that looks similar to the following:

   ```groovy title="main.nf"
   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       <tool1>: \$(tool1 --version)
       <tool2>: \$(tool2 --version)
   END_VERSIONS
   ```

1. Remove the `versions.yml` file from the `output` block.

1. Add outputs for each tool in the module:

   ```groovy
   tuple val("${task.process}"), val('<tool1>'), eval('tool1 --version'), emit: versions_tool1, topic: versions
   ```

   Replace `<tool1>` and `tool1 --version` with the tool name and version command.

   - Repeat this for each tool used in the module.

1. Run `nf-core modules lint <module_name> --fix` to migrate the `meta.yml` file with the new topic outputs.

1. Add a `type` and `description` for each field in the versions output:

   ```yaml title="meta.yml"
   versions_<tool>:
     - - ${task.process}:
           type: string
           description: The process the versions were collected from
       - <tool>:
           type: string
           description: The tool name
       - <versions_command>:
           type: string
           description: The version of the tool
   ```

   Update `<tool>` and `<versions_command>` to the output channel values.

1. Add the topics block to the `meta.yml`file below the outputs section:

   ```yaml title="meta.yml"
   topics:
     - versions:
         - - process:
               type: string
               description: The process the versions were collected from
           - tool:
               type: string
               description: The tool name
           - version:
               type: string
               description: The version of the tool
   ```

1. Update the `main.nf.test` file to check for the new version outputs.

   - If the test runs on all process output (`snapshot(process.out).match()`), do nothing.
   - If the test checks for specific outputs, update it to check for the new version outputs. `process.out.versions` should be changed to `process.out.findAll { key, val -> key.startsWith('versions') }`.

1. Run tests to regenerate the snapshots:

   ```bash
   nf-core modules test <module_name> --update
   ```

1. Check that the snapshot is correct and that all versions are being captured correctly.

   - If not, return to step 3 and migrate any incorrect information in the module.

1. Commit and push changes to your fork.

1. Open a pull request to the main `nf-core/modules` repository.

## Migrate modules with template scripts

To migrate modules that use template scripts to use topic channels for version outputs:

1. Open the modules `main.nf` file.

1. Update the `path "versions.yml", emit: versions` line:

   ```groovy title="main.nf"
   path "versions.yml", emit: versions, topic: versions
   ```

1. Add the following to the `meta.yml` file below the outputs section:

   ```yaml title="meta.yml"
   versions:
     - versions.yml:
         type: file
         description: YAML file containing versions of tools used in the module
   ```

1. Commit and push changes to your fork.

1. Open a pull request to the main `nf-core/modules` repository.

## Migrate pipelines

To migrate a pipelines to use topic channels for version outputs:

1. Update the pipeline template a version that includes support for topic channels (version 3.5.0 or later).

1. Run `nf-core modules update` to pull the latest changes.

1. Find and remove instances in the pipeline code where `versions` is referenced as an invalid process output. These outputs are now handled by topic channels and should be removed.

   You may encounter errors like this:

   ```console
   ERROR ~ No such variable: Exception evaluating property 'versions' for nextflow.script.ChannelOut, Reason: groovy.lang.MissingPropertyException: No such property: versions for class: groovyx.gpars.dataflow.DataflowBroadcast
   ```

   **Example**: A workflow using the `samtools/sort` module might have code like this:

   ```nextflow
   SAMTOOLS_SORT(input)
   ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
   ```

   **Fix**: Remove the line referencing `SAMTOOLS_SORT.out.versions`:

   ```diff
   SAMTOOLS_SORT(input)
   -ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
   ```
