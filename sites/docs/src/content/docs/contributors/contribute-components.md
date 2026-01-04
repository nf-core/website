---
title: Contributing components
subtitle: How to contribute modules and subworkflows to nf-core
shortTitle: Contributing components
---

<!-- TODO: Add links -->
<!-- TODO: Revise is test data and testing should be expanded here, my 2 cents is is it should be another page -->

nf-core modules and subworkflows are reusable Nextflow components shared across pipelines:

- **Modules**: Individual tool wrappers that represent single processes
- **Subworkflows**: Combinations of multiple modules into cohesive analysis units

Shared components enable standardized, reproducible analyses across research groups worldwide, reduce duplicated development effort, and accelerate scientific discovery. The nf-core community strongly encourages all members to contribute their expertise and welcomes contributions of any scale.

This guide explains how to contribute modules and subworkflows to the [nf-core/modules repository](https://github.com/nf-core/modules). For detailed information about writing components, see[Writing Components](LINK).

:::note{title="Prerequisites"}
You will need the following to get started:

- [Nextflow version 21.04.0 or later](../../get_started/environment_setup/nextflow.md)
- [nf-core/tools version 2.7 or later](../../get_started/environment_setup/nf-core-tools.md)
- [nf-test](https://www.nf-test.com/) for testing
- [pre-commit](https://pre-commit.com/) for code quality checks
- A container engine ([Docker](../../get_started/environment_setup/software-dependencies.md#docker), [Singularity](../../get_started/environment_setup/software-dependencies.md#singularity), or [Conda](../../get_started/environment_setup/software-dependencies.md#condamamba))
- A [GitHub account](https://github.com/signup)

:::

## Contribute a component

1. Check if the component already exists:

    ```bash
    nf-core modules list
    ```

    Also check:

    - The [nf-core/modules repository](https://github.com/nf-core/modules)
    - Open [pull requests](https://github.com/nf-core/modules/pulls)
    - Open [issues](https://github.com/nf-core/modules/issues)

1. Create an [issue](https://github.com/nf-core/modules/issues) to let others know you're working on it:

    - Use a clear title, e.g., "New module: fastqc"
    - Describe what the component will do
    - Add yourself to the **Assignees**

1. Fork the repository to your GitHub account and clone it locally:

    ```bash
    git clone https://github.com/<github_username>/modules.git nf-core-modules
    cd nf-core-modules
    git remote add upstream https://github.com/nf-core/modules.git
    git checkout -b <component_name>
    pre-commit install
    ```

    :::note
    The `pre-commit install` command sets up automatic code quality checks.
    :::

1. Write your component:

    - Use nf-core/tools to generate the component structure:

        - **Modules:**

          ```bash
          nf-core modules create
          ```

        - **Subworkflows:**

          ```bash
          nf-core subworkflows create
          ```

    :::note
    The `create` command creates the necessary files and directory structure:

    - **`main.nf`**: The main script containing the process definition with TODO statements to guide you
    - **`meta.yml`**: Metadata file storing general information about the module and input/output descriptions
    - **`tests/main.nf.test`**: Test workflow to unit test the module outputs (required for all modules)

    For detailed information about writing modules and subworkflows, see the [writing components guide](LINK).
    :::

1. Before submitting, thoroughly test your component:

    ```bash
    nf-core modules lint <component_name>
    ```

    ```bash
    nf-core modules test <component_name>
    ```

    :::note
    The linting checks ensure your component follows nf-core standards. All tests must pass before submission. See LINK for more information about lint and test.
    :::

1. Commit and push your changes with a clear message:

    ```bash
    git add .
    git commit -m "Add new module: <component_name>"
    git push origin <component_name>
    ```

    :::tip
    Before creating a pull request, sync your branch with upstream to apply your changes on top of the latest updates:

    ```bash
    git pull --rebase upstream master
    ```
    :::
1. Go to the [nf-core/modules repository](https://github.com/nf-core/modules) and create a pull request from your branch:

    - Use the pull request template provided by the repository:

        - Reference the issue you created
        - Describe what the component does
        - Explain how you tested it
        - Include example commands or screenshots if helpful

1. Address review feedback

    :::note
    The nf-core maintainers will review your contribution and may request changes. Common feedback includes:

    - Code style improvements
    - Additional tests
    - Documentation clarifications
    - Container or dependency updates
    :::

1. When approved, merge your pull request

## Additional resources

If you need assistance:

- Ask in the `#modules` channel on [nf-core Slack](https://nf-co.re/join)
- Check existing components in the [nf-core/modules repository](https://github.com/nf-core/modules) for examples
- Review the [writing components guide](LINK) for technical details
- Consult the [nf-core contributing guidelines](https://nf-co.re/docs/contributing/guidelines)

The nf-core community is here to help. Don't hesitate to ask questions.

<!-- TODO: Add links to -->
