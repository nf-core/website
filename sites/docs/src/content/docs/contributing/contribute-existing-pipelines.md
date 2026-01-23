---
title: Contributing to existing pipelines
subtitle: How to contribute to existing nf-core pipelines
shortTitle: Existing pipelines
---

nf-core pipelines are community-driven bioinformatics workflows that enable reproducible analyses across research groups.
Contributing to existing pipelines helps improve functionality, fix bugs, and add new features that benefit the entire community.

This guide explains how to contribute code to existing nf-core pipelines.
Whether you're adding a new feature, fixing a bug, or improving documentation, the nf-core community welcomes your contributions.

:::note{title="Prerequisites"}
You will need the following to get started:

- [Nextflow version 21.04.0 or later](../../get_started/environment_setup/nextflow.md)
- [nf-core/tools version 1.10 or later](../../get_started/environment_setup/nf-core-tools.md)
- A container engine ([Docker](../../get_started/environment_setup/software-dependencies.md#docker), [Singularity](../../get_started/environment_setup/software-dependencies.md#singularity), or [Conda](../../get_started/environment_setup/software-dependencies.md#condamamba))
- [Git](https://git-scm.com/) for version control
- A [GitHub account](https://github.com/signup)

:::

## Contribute to an existing pipeline

1. Check if someone is already working on your idea:

   Check the pipeline repository for:
   - Open [issues](https://github.com/nf-core) to see if your idea has been discussed
   - Open [pull requests](https://github.com/nf-core) to avoid duplicating work

1. Create an [issue](https://github.com/nf-core) to let others know you're working on it:
   - Use a clear title that describes your contribution
   - Describe what you plan to add or change
   - Add yourself to the **Assignees**

1. Fork the pipeline repository to your GitHub account and clone it locally:

   ```bash
   git clone https://github.com/<github_username>/<pipeline_name>.git
   cd <pipeline_name>
   git remote add upstream https://github.com/nf-core/<pipeline_name>.git
   git checkout -b <branch_name>
   ```

   :::note
   Replace `<github_username>`, `<pipeline_name>`, and `<branch_name>` with appropriate values.
Use descriptive branch names like `add-fastqc-module` or `fix-memory-issue`.
   :::

1. Make your changes following the [pipeline conventions](#pipeline-contribution-conventions):
   - Make changes on the `dev` branch
   - Follow the pipeline's existing code style and structure
   - Add or update tests as needed
   - Update documentation to reflect your changes

1. Add new parameters to the pipeline schema:

   If you added new parameters, update the JSON schema:

   ```bash
   nf-core pipelines schema build
   ```

   :::note
   This command updates `nextflow_schema.json` with your new parameters and provides an interactive interface to add descriptions and validation rules.
   :::

1. Test your changes locally:

   ```bash
   nextflow run . -profile debug,test,docker --outdir <OUTDIR>
   ```

   :::tip
   The `debug` profile provides warnings about process selectors and other debug information that can help identify issues.
   :::

1. Run linting checks to ensure code quality:

   ```bash
   nf-core pipelines lint .
   ```

   :::note
   Address any failures or warnings.
Follow the listed URLs for documentation on how to fix specific issues.
   :::

1. Commit and push your changes with a clear message:

   ```bash
   git add .
   git commit -m "Add feature: <description>"
   git push origin <branch_name>
   ```

   :::tip
   Before creating a pull request, sync your branch with upstream to apply your changes on top of the latest updates:

   ```bash
   git pull --rebase upstream dev
   ```

   :::

1. Submit a pull request against the `dev` branch:

   Go to the pipeline repository on GitHub and create a pull request:
   - Use the pull request template provided
   - Reference the issue you created
   - Describe what you changed and why
   - Explain how you tested your changes
   - Include example commands or screenshots if helpful

1. Address review feedback

   :::note
   The nf-core maintainers and automated tests will review your contribution.
Common feedback includes:
   - Code style improvements
   - Additional or updated tests
   - Documentation clarifications
   - Linting fixes
   - Pipeline-specific conventions

   Pull requests are typically fully reviewed when automated tests pass.
   :::

1. When approved, a maintainer will merge your pull request

## Automated testing

When you create a pull request, [GitHub Actions](https://github.com/features/actions) automatically runs tests to ensure code quality and functionality.
Pull requests are typically fully reviewed when these tests pass.

### Lint tests

All nf-core pipelines must adhere to a [set of guidelines](https://nf-co.re/developers/guidelines).
The linting tool checks that your code follows these standards:

```bash
nf-core pipelines lint <pipeline-directory>
```

The tool verifies:

- Code structure and organization
- Documentation completeness
- Configuration file formatting
- Pipeline metadata

:::note
If any failures or warnings occur, follow the provided URLs for detailed documentation on how to fix each issue.
:::

### Pipeline tests

Each nf-core pipeline includes a minimal test dataset.
GitHub Actions runs the pipeline on this data to verify successful execution:

- Tests run with the latest version of Nextflow
- Tests run with the minimum required Nextflow version specified in the pipeline

If tests fail, review the GitHub Actions logs to identify and fix issues before requesting review.

## Patching bugs in releases

:::warning
Only use this workflow if a critical bug is discovered in a released version and needs an immediate patch.
:::

If a release contains a critical bug that requires immediate fixing:

1. Create a new `patch` branch based on `upstream/master`
1. Fix the bug and bump the patch version (X.Y.Z+1)
1. Submit a pull request to `master` (not `dev`) to directly resolve the bug

This workflow bypasses the normal development cycle for urgent production fixes.

## Pipeline contribution conventions

Follow these conventions to ensure code quality and consistency across nf-core pipelines.

### Adding a new process

When adding a new process to a pipeline:

1. Define the input channel from the expected previous process
1. Write the process block following pipeline conventions
1. Define output channels as needed
1. Add new parameters to `nextflow.config` with sensible defaults
1. Update `nextflow_schema.json` using `nf-core pipelines schema build`
1. Add validation for all relevant parameters
1. Test locally to verify the new process works correctly
1. Add test commands to `.github/workflow/ci.yml` if applicable
1. Update MultiQC configuration in `assets/multiqc_config.yml`:
   - Add relevant file suffixes and cleanup rules
   - Ensure module plots appear in the correct order
1. Document output files in `docs/output.md` with MultiQC screenshots if relevant

### Parameters

Define parameters with default values in `nextflow.config` under the `params` scope:

```groovy
params {
    my_parameter = 'default_value'
}
```

After adding parameters, update the schema:

```bash
nf-core pipelines schema build
```

This adds your parameters to `nextflow_schema.json` with descriptions and validation rules.

### Resource requirements

Define process resource requirements (CPUs, memory, time) in `conf/base.config` using `withLabel:` selectors:

```groovy
withLabel: process_low {
    cpus = 2
    memory = 4.GB
}
```

Use standardized labels from the [nf-core pipeline template](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/conf/base.config).
Reference resources dynamically in process blocks:

```groovy
script:
"""
tool --threads ${task.cpus} --memory ${task.memory}
"""
```

### Channel naming

Use consistent naming schemes for channels:

- Initial process output: `ch_output_from_<process>`
- Intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

Example:

```groovy
ch_output_from_fastqc
ch_fastqc_for_multiqc
```

### Nextflow version requirements

If you use new Nextflow features, update the minimum required version:

```bash
nf-core pipelines bump-version --nextflow . [min_nf_version]
```

### Images and figures

Follow the nf-core [style guidelines](https://nf-co.re/developers/design_guidelines) for pipeline diagrams, workflow images, and other visual documentation.

## Additional resources

If you need assistance:

- Ask in the pipeline-specific channel on [nf-core Slack](https://nf-co.re/join) (e.g., `#<pipeline-name>`)
- Consult the pipeline's documentation and usage guides
- Review existing code in the pipeline for examples
- Check the [nf-core contributing guidelines](https://nf-co.re/docs/contributing/guidelines)

The nf-core community is here to help.
Don't hesitate to ask questions.
