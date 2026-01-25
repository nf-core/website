---
title: Adding modules to pipelines
subtitle: How to add nf-core modules to pipelines
shortTitle: Adding modules
weight: 1
---

nf-core modules are standardized, reusable components that wrap individual bioinformatics tools.
Each module includes the tool's process definition, container specifications, and metadata.
By using modules from the central nf-core/modules repository, you benefit from community testing, standardized structure, and easier maintenance.

You can use nf-core tools to add these modules directly to your pipeline.
Follow these steps to integrate a new nf-core module.

:::note
Pipeline workflows may vary.
This guide covers the most common integration steps.
:::

## Add an nf-core module

To add an nf-core module to your pipeline with nf-core tools:

1. Install the module to your pipeline using the nf-core command-line tool:

   ```bash
   nf-core modules install <tool>/<subtool>
   ```

1. Add the `include` statement at the workflow level to import the module.
   - The installation command from the previous step will provide the exact syntax you need.

1. Insert the module call in the relevant workflow script.
   - Ensure version information is mixed into the versions channel and route output files to MultiQC if applicable.

1. Create a dedicated section in `conf/modules.conf` that specifies the default output directory and file patterns for the module.

1. Add module-specific parameters in `nextflow.config` with sensible defaults.
   - For pipeline-level parameters that need to be passed to the module, insert them into `ext.args` using proper closure syntax:

     ```groovy
     ext.args = { "--option $params.option" }
     ```

1. Update `nextflow_schema.json` to include the new parameters:

   ```bash
   nf-core pipelines schema build
   ```

1. Update the following documentation files:
   - Add tool citations to `CITATIONS.md`
   - Document usage information in `docs/usage.md`
   - Describe output directories in `docs/output.md`
   - Update the pipeline README and workflow diagrams
   - Update the MultiQC configuration file

1. Run the nf-core linting tool to ensure your code meets nf-core standards:

   ```bash
   nf-core pipelines lint
   ```

1. Run Prettier to format your code.

   ```bash
   prettier -w
   ```

1. Run the pipeline locally with test data to verify:
   - The module executes successfully
   - Outputs are generated correctly
   - No errors or warnings occur
