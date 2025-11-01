---
title: Update modules to use topic channels
description: How to update nf-core modules to use the new topic channels feature
---

In this tutorial, we will go through the steps needed to update existing nf-core modules to use the new topic channels feature for version outputs. This will be the main way to update most modules, the only exceptions being modules that use template scripts. Take a look at the [template migration](#template-migration) chapter for more information on how to update those modules.

## Prerequisites

- The latest version of `nf-core/tools` (version 3.5.0 or higher)
- The modules repository is cloned locally or open in some development environment.
- Nextflow version 25.04.0 or higher is installed.

## Updating standard modules

To update a module to use topic channels for version outputs, follow these steps:

1. Go the the `main.nf` file of the module you want to update.

   Each one of these files will have a piece of code that looks like this:

   ```groovy
   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       <tool1>: \$(tool1 --version)
       <tool2>: \$(tool2 --version)
   END_VERSIONS
   ```

2. Remove the `versions.yml` file from the `output` block.

3. Add the following new outputs for each tool used in the module:

   ```groovy
   tuple val("${task.process}"), val('<tool1>'), eval('tool1 --version'), emit: versions_tool1, topic: versions
   ```

   Update `<tool1>` and `tool1 --version` with the actual tool name and version command. Repeat this for each tool used in the module.

4. Run `nf-core modules lint <module_name> --fix` to automatically update the `meta.yml` file with the new topic outputs.

5. Add a `type` and `description` field for each field in the versions output. You can use the following template for this:

   ```yaml
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

   Update `<tool>`, `<versions_command>` to the actual values in the output channel.

6. Add the topics block to the `meta.yml`file. Ideally it would be located underneath the outputs section:

```yaml
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

7. Update the `main.nf.test` file to check for the new versions outputs.

- If the test runs on all process output (`snapshot(process.out).match()`), you don't have to do anything.
- If the test checks for specific outputs, you will need to update it to check for the new version outputs. `process.out.versions` should be changed to `process.out.findAll { key, val -> key.startsWith('versions') }`.

8. Run the tests to regenerate the snapshots:

   ```bash
   nf-core modules test <module_name> --update
   ```

9. Check that the snapshot is correct and that all versions are being captured correctly. If not, please return to step 3 and update anything wrong with the module.

10. Commit and push your changes to your fork and open a pull request to the main nf-core/modules repository.

## Template migration

Modules that use template scripts for version outputs will need to be updated slightly differently. Follow these steps:

1. Go to the `main.nf` file of the module you want to update.

2. Update the `path "versions.yml", emit: versions` line to this:

   ```groovy
   path "versions.yml", emit: versions, topic: versions
   ```

3. Add the following lines to the `meta.yml` file underneath the outputs section:

   ```yaml
   versions:
     - versions.yml:
         type: file
         description: YAML file containing versions of tools used in the module
   ```

4. Commit and push your changes to your fork and open a pull request to the main nf-core/modules repository.
