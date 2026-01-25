---
title: Create components
subtitle: Writing modules and subworkflows for nf-core
shortTitle: Create components
---

nf-core modules and subworkflows are reusable Nextflow components shared across pipelines.
Modules are individual tool wrappers that represent single processes, while subworkflows combine multiple modules into cohesive analysis units.
Shared components enable standardized, reproducible analyses across research groups worldwide, reduce duplicated development effort, and accelerate scientific discovery.

This guide explains how to write and develop modules and subworkflows for nf-core.
Once your component is ready, see [Contributing Components](../../contributors/contribute-components.md) for instructions on contributing it to nf-core.

:::note{title="Prerequisites"}
You will need the following to get started:

- [Nextflow version 21.04.0 or later](../../get_started/environment_setup/nextflow.md)
- [nf-core/tools version 2.7 or later](../../get_started/environment_setup/nf-core-tools.md)
- [nf-test](https://www.nf-test.com/) for testing
- [pre-commit](https://pre-commit.com/) for code quality checks
- A container engine ([Docker](../../get_started/environment_setup/software-dependencies.md#docker), [Singularity](../../get_started/environment_setup/software-dependencies.md#singularity), or [Conda](../../get_started/environment_setup/software-dependencies.md#condamamba))

:::

## Check if the component exists

Before creating a new component, verify it doesn't already exist by:

1. Using the `nf-core modules list` command:

   ```bash
   nf-core modules list
   ```

1. Checking the [nf-core/modules repository](https://github.com/nf-core/modules).

1. Reviewing open [pull requests](https://github.com/nf-core/modules/pulls) and [issues](https://github.com/nf-core/modules/issues).

If the component is new, create an [issue](https://github.com/nf-core/modules/issues) with an appropriate title (e.g., "new module: fastqc") and assign yourself to let others know you're working on it.

## Create the component structure

To create a new module or subworkflow:

1. Create a new branch in your local repository:

   ```bash
   git checkout -b <component-name>
   ```

1. Use nf-core/tools to generate the component structure:
   - For modules:

     ```bash
     nf-core modules create
     ```

   - For subworkflows:

     ```bash
     nf-core subworkflows create
     ```

   The create command generates three core files with a standard structure:
   - **`main.nf`**: The main script containing the process or workflow definition with TODO statements to guide you through implementation.
   - **`meta.yml`**: Metadata file storing general information, author details, and input/output descriptions. This file is validated against a JSON schema to ensure consistency.
   - **`tests/main.nf.test`**: Test workflow to unit test the component outputs. All components must include tests.

## Write your component

To implement your component:

1. Follow the TODO statements in the generated `main.nf` file.
   - Implement the tool command or workflow logic.
   - Populate version channels using eval output qualifiers.
   - Ensure adherence to [component specifications](LINK).

1. Complete the `meta.yml` file with comprehensive metadata.
   - Fill in all required metadata fields.
   - Document all inputs and outputs with descriptions.
   - Include tool references and version information.

1. Write comprehensive tests in `tests/main.nf.test`.
   - Use minimal test data from `tests/config/test_data.config` when possible.
   - Follow test-data specifications to keep test datasets small and fast.
   - Create test snapshots to verify outputs match expected results.
   - Add multiple test scenarios if applicable (e.g., single-end and paired-end data).
   - For large datasets requiring stub tests, document alternative testing procedures.
   - See the [Testing documentation](../testing/overview.md) for comprehensive guides on writing assertions and handling complex testing scenarios.

:::note
All components require test workflows.
Tests validate that your component works correctly and prevent regressions when changes are made.
:::

## Test your component

Before contributing your component, thoroughly test it to ensure it meets nf-core standards:

1. Run the linting tool to check code quality and adherence to nf-core standards:
   - For modules:

     ```bash
     nf-core modules lint <module_name>
     ```

   - For subworkflows:

     ```bash
     nf-core subworkflows lint <subworkflow_name>
     ```

   The linting tool checks for common issues such as missing files, incorrect formatting, and invalid metadata.

1. Run the component tests:
   - For modules:

     ```bash
     nf-core modules test <module_name>
     ```

   - For subworkflows:

     ```bash
     nf-core subworkflows test <subworkflow_name>
     ```

   All tests must pass before you can contribute your component.

:::tip
GitHub Actions automatically runs these tests when you submit a pull request, but running them locally first helps catch issues early.
:::

## Next steps

Once your component is written and all tests pass, you can contribute it to nf-core:

- Follow the [Contributing Components](../../contributors/contribute-components.md) guide to submit your component to the nf-core/modules repository.
- Your component will be reviewed by the nf-core community and, once approved, will be available to all nf-core pipelines and the Nextflow community.

## Additional resources

If you need assistance while developing your component:

- Ask in the `#modules` channel on [nf-core Slack](https://nf-co.re/join).
- Check existing components in the [nf-core/modules repository](https://github.com/nf-core/modules) for examples.

The nf-core community is here to help.
Don't hesitate to ask questions.
