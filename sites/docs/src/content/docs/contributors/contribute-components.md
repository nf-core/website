---
title: Contributing components
subtitle: How to contribute modules and subworkflows to nf-core
shortTitle: Contributing components
---

<!-- TODO: Add links -->
<!-- TODO: Revise is test data and testing should be expanded here, my 2 cents is is it should be another page -->

This guide explains how to contribute modules and subworkflows to the [nf-core/modules repository](https://github.com/nf-core/modules). For detailed information about writing components, see the [writing components guide](LINK).

## What are nf-core components?

nf-core modules and subworkflows are reusable Nextflow components shared across pipelines:

- **Modules**: Individual tool wrappers that represent single processes
- **Subworkflows**: Combinations of multiple modules into cohesive analysis units

All components in nf-core/modules are reviewed, tested, and maintained by the community.

## Prerequisites

Before contributing, ensure you have:

- [Nextflow version 21.04.0 or later](../../get_started/environment_setup/nextflow.md)
- [nf-core/tools version 2.7 or later](../../get_started/environment_setup/nf-core-tools.md)
- [nf-test](https://www.nf-test.com/) for testing
- [pre-commit](https://pre-commit.com/) for code quality checks
- A container engine ([Docker](../../get_started/environment_setup/software-dependencies.md#docker), [Singularity](../../get_started/environment_setup/software-dependencies.md#singularity), or [Conda](../../get_started/environment_setup/software-dependencies.md#condamamba))
- A [GitHub account](https://github.com/signup)

## Contribution workflow

### 1. Check if the component already exists

Before starting, verify the component doesn't already exist:

```bash
nf-core modules list
```

Also check:

- The [nf-core/modules repository](https://github.com/nf-core/modules)
- Open [pull requests](https://github.com/nf-core/modules/pulls)
- Open [issues](https://github.com/nf-core/modules/issues)

### 2. Create an issue

If the component doesn't exist, create an [issue](https://github.com/nf-core/modules/issues) to let others know you're working on it:

- Use a clear title, e.g., "New module: fastqc"
- Describe what the component will do
- Add yourself to the **Assignees**

This prevents duplicate work and allows the community to provide feedback early.

### 3. Fork and clone the repository

Fork the repository to your GitHub account and clone it locally:

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

### 4. Write your component

Use nf-core/tools to generate the component structure:

- **For modules:**

  ```bash
  nf-core modules create
  ```

- **For subworkflows:**

  ```bash
  nf-core subworkflows create
  ```

The `create` command creates the necessary files and directory structure:

- **`main.nf`**: The main script containing the process definition with TODO statements to guide you
- **`meta.yml`**: Metadata file storing general information about the module and input/output descriptions
- **`tests/main.nf.test`**: Test workflow to unit test the module outputs (required for all modules)

For detailed information about writing modules and subworkflows, see the [writing components guide](LINK).

### 5. Test your component

Before submitting, thoroughly test your component:

```bash
nf-core modules lint <component_name>
nf-core modules test <component_name>
```

The linting checks ensure your component follows nf-core standards. All tests must pass before submission.

See LINK for more information about lint and test.

### 6. Commit and push your changes

Commit your changes with a clear message:

```bash
git add .
git commit -m "Add new module: <component_name>"
git push origin <component_name>
```

The pre-commit hooks will automatically check your code before committing.

:::tip
Before creating a pull request, sync your branch with upstream to apply your changes on top of the latest updates:

```bash
git pull --rebase upstream master
```
:::

### 7. Create a pull request

Go to the [nf-core/modules repository](https://github.com/nf-core/modules) and create a pull request from your branch. In the pull request:

- Reference the issue you created
- Describe what the component does
- Explain how you tested it
- Include example commands or screenshots if helpful

Use the pull request template provided by the repository.

### 8. Respond to review feedback

The nf-core maintainers will review your contribution and may request changes. Common feedback includes:

- Code style improvements
- Additional tests
- Documentation clarifications
- Container or dependency updates

Address feedback by pushing new commits to your branch. The pull request will automatically update.

### 9. After merge

Once your pull request is approved and merged:

- Your component becomes available to all nf-core users
- It will be automatically tested with every nf-core/modules update
- You will be listed as a contributor to nf-core/modules

## Additional resources

If you need assistance:

- Ask in the `#modules` channel on [nf-core Slack](https://nf-co.re/join)
- Check existing components in the [nf-core/modules repository](https://github.com/nf-core/modules) for examples
- Review the [writing components guide](LINK) for technical details
- Consult the [nf-core contributing guidelines](https://nf-co.re/docs/contributing/guidelines)

The nf-core community is here to help. Don't hesitate to ask questions.

<!-- TODO: Add links to resources -->
