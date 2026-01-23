---
title: Contributing components
subtitle: How to contribute modules and subworkflows to nf-core
shortTitle: Components
---

nf-core modules and subworkflows are reusable Nextflow components shared across pipelines.
Modules are individual tool wrappers that represent single processes, while subworkflows combine multiple modules into cohesive analysis units.
Shared components enable standardized, reproducible analyses across research groups worldwide, reduce duplicated development effort, and accelerate scientific discovery.

Follow these steps to contribute modules and subworkflows to the [nf-core/modules repository](https://github.com/nf-core/modules), making them available to all nf-core pipelines and the Nextflow community.
The nf-core community strongly encourages all members to contribute their expertise and welcomes contributions of any scale.

:::note{title="Prerequisites"}
You will need the following to get started:

- [Nextflow version 21.04.0 or later](../get_started/environment_setup/nextflow.md)
- [nf-core/tools version 2.7 or later](../get_started/environment_setup/nf-core-tools.md)
- [nf-test](https://www.nf-test.com/) for testing
- [pre-commit](https://pre-commit.com/) for code quality checks
- A container engine ([Docker](../get_started/environment_setup/software-dependencies.md#docker), [Singularity](../get_started/environment_setup/software-dependencies.md#singularity), or [Conda](../get_started/environment_setup/software-dependencies.md#condamamba))
- A [GitHub account](https://github.com/signup)

:::

:::tip
For detailed information about writing components, see [Create Components](../developers/components/create-components.md).
:::

## Check if the component exists

Before starting, verify the component doesn't already exist:

1. Use the `nf-core modules list` command to search for existing modules:

   ```bash
   nf-core modules list
   ```

1. Check the [nf-core/modules repository](https://github.com/nf-core/modules) for existing components.

1. Review open [pull requests](https://github.com/nf-core/modules/pulls) to see if someone is already working on the component.

1. Search open [issues](https://github.com/nf-core/modules/issues) for related discussions.

## Create an issue

To let others know you're working on a new component:

1. Create an [issue](https://github.com/nf-core/modules/issues) in the nf-core/modules repository.
   - Use a clear title, e.g., "New module: fastqc".
   - Describe what the component will do and how it will be used.
   - Add yourself to the **Assignees** to indicate you're working on it.

This helps prevent duplicate work and allows the community to provide early feedback.

## Fork and clone the repository

To set up your development environment:

1. Fork the [nf-core/modules repository](https://github.com/nf-core/modules) to your GitHub account.
   - Go to the nf-core/modules repository and select **Fork**.

1. Clone your fork locally and set up the upstream remote:

   ```bash
   git clone https://github.com/<github_username>/modules.git nf-core-modules
   cd nf-core-modules
   git remote add upstream https://github.com/nf-core/modules.git
   git checkout -b <component_name>
   ```

1. Install pre-commit hooks:

   ```bash
   pre-commit install
   ```

:::note
The `pre-commit install` command sets up automatic code quality checks that will run when you commit changes.
These checks help ensure your code meets nf-core standards.
:::

## Create your component

To create the component structure:

1. Use nf-core/tools to generate the component files:
   - For modules:

     ```bash
     nf-core modules create
     ```

   - For subworkflows:

     ```bash
     nf-core subworkflows create
     ```

   This command creates the necessary files and directory structure with template code to guide you.

1. Write your component following the [Create Components guide](../developers/components/create-components.md).
   - Implement the tool command or workflow logic in `main.nf`.
   - Complete the metadata in `meta.yml`.
   - Write comprehensive tests in `tests/main.nf.test`.

## Test your component

Before submitting, thoroughly test your component:

1. Run the linting tool to check code quality:
   - For modules:

     ```bash
     nf-core modules lint <module_name>
     ```

   - For subworkflows:

     ```bash
     nf-core subworkflows lint <subworkflow_name>
     ```

   The linting checks ensure your component follows nf-core standards and includes all required files.

1. Run the component tests:
   - For modules:

     ```bash
     nf-core modules test <module_name>
     ```

   - For subworkflows:

     ```bash
     nf-core subworkflows test <subworkflow_name>
     ```

   All tests must pass before submission.

:::note
GitHub Actions automatically runs these tests when you submit a pull request, but running them locally first helps catch issues early and speeds up the review process.
:::

## Submit your component

To submit your component for review:

1. Commit your changes with a clear message:

   ```bash
   git add .
   git commit -m "Add new component: <component_name>"
   ```

1. Sync your branch with upstream before pushing:

   ```bash
   git pull --rebase upstream master
   git push origin <component_name>
   ```

:::tip
Rebasing with upstream before creating a pull request applies your changes on top of the latest updates and helps prevent merge conflicts.
:::

## Create a pull request

To request review of your component:

1. Go to the [nf-core/modules repository](https://github.com/nf-core/modules) on GitHub.

1. Create a pull request from your branch.
   - Use the pull request template provided by the repository.
   - Reference the issue you created.
   - Describe what the component does and how you tested it.
   - Include example commands or screenshots if helpful.

1. Add the "Ready for Review" label to indicate your component is ready for maintainer review.

1. Request reviews from `nf-core/modules-team`.

:::note
Components are tested via GitHub Actions CI using Docker, Singularity, and Conda to ensure they work across different environments.
The automated tests must pass before your component can be merged.
:::

## Address review feedback

The nf-core maintainers will review your contribution:

1. Wait for maintainer feedback on your pull request.

   Common feedback includes:
   - Code style improvements.
   - Additional tests or test scenarios.
   - Documentation clarifications.
   - Container or dependency updates.

1. Make requested changes and push them to your branch:

   ```bash
   git add .
   git commit -m "Address review feedback"
   git push origin <component_name>
   ```

   The pull request will update automatically with your new changes.

1. Respond to reviewer comments to discuss any questions or clarifications.

## Complete the contribution

Once your pull request is approved:

1. A maintainer will merge your pull request.

1. Your component will be available in the nf-core/modules repository.

1. All nf-core pipelines and Nextflow users can now use your component.

Congratulations on contributing to nf-core! Your work helps researchers worldwide conduct reproducible analyses.

## Additional resources

If you need assistance during the contribution process:

- Ask in the `#modules` channel on [nf-core Slack](https://nf-co.re/join).
- Check existing components in the [nf-core/modules repository](https://github.com/nf-core/modules) for examples.
- Review the [Create Components guide](../developers/components/create-components.md) for technical details.

The nf-core community is here to help.
Don't hesitate to ask questions.
