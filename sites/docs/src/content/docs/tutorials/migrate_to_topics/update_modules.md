---
title: Update modules to use topic channels
description: How to update nf-core modules to use the new topic channels feature
---

In this tutorial, we will go through the steps needed to update existing nf-core modules to use the new topic channels feature for version outputs. This will be the main way to update most modules, the only exceptions being modules that use template scripts. Take a look at the [template migration](#template-migration) chapter for more information on how to update those modules.

## Prerequisites

- nf-core tools version 3.5.0 or higher
- A clone of the `nf-core/modules` repository
- Nextflow version 25.04.0 or higher

## Updating standard modules

To update a module to use topic channels for version outputs:

1. Open the modules `main.nf` file.

1. Identify the code that looks similar to the following:

   ```groovy title="main.nf"
   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       <tool1>: \$(tool1 --version)
       <tool2>: \$(tool2 --version)
   END_VERSIONS
   ```

1. Remove the `versions.yml` file from the `output` block.

1. Add the following outputs for each tool in the module:

   ```groovy title="main.nf"
   tuple val("${task.process}"), val('<tool1>'), eval('tool1 --version'), emit: versions_tool1, topic: versions
   ```

   Update `<tool1>` and `tool1 --version` with the tool name and version command. Repeat this for each tool used in the module.

1. Run

   ```bash
   nf-core modules lint <module_name> --fix
   ```

   to automatically update the `meta.yml` file with the new topic outputs.

1. Add a `type` and `description` field for each field in the versions output. Use the following template:

   ```yaml title="meta.yml"
   - - "${task.process}":
          type: string
          description: The name of the process
      - "{{ component }}":
          type: string
          description: The name of the tool
      - "{{ component }} --version":
          type: eval
          description: The expression to obtain the version of the tool
   ```

   Update `<tool>` and `<versions_command>` to the output channel values.

1. Add the topics block to the `meta.yml`file underneath the outputs section:

   ```yaml title="meta.yml"
   topics:
     versions:
       - - "${task.process}":
           type: string
           description: The name of the process
       - "{{ component }}":
           type: string
           description: The name of the tool
       - "{{ component }} --version":
           type: eval
           description: The expression to obtain the version of the tool
   ```

1. Update the `main.nf.test` file to check for the new version outputs.
   - If the test runs on all process output (`snapshot(process.out).match(){:groovy}`), do nothing.
   - If the test checks for specific outputs, update it to check for the new version outputs. `process.out.versions{:groovy}` should be changed to `process.out.findAll { key, val -> key.startsWith('versions') }{:groovy}`.

1. Run the tests to regenerate the snapshots:

   ```bash
   nf-core modules test <module_name> --update
   ```

1. Check that the snapshot is correct and that all versions are being captured correctly. If not, return to step 3 and update any incorrect information in the module.

1. Commit and push your changes to your fork.

1. Open a pull request to the main `nf-core/modules` repository.

## Template migration

Modules that use template scripts for version outputs will need to be updated slightly differently:

1. Open the modules `main.nf` file.

1. Update the `path "versions.yml", emit: versions{:groovy}` line to the following:

   ```groovy title="main.nf"
   path "versions.yml", emit: versions, topic: versions
   ```

1. Add the following lines to the `meta.yml` file underneath the outputs section:

   ```yaml title="meta.yml"
   versions:
     - versions.yml:
         type: file
         description: YAML file containing versions of tools used in the module
   ```

1. Commit and push your changes to your fork.

1. Open a pull request to the `nf-core/modules` repository.
