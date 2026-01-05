---
title: Contributing a new pipeline
subtitle: How to contribute a new nf-core pipeline
shortTitle: New pipelines
---

nf-core pipelines use a standardized template that provides file structure, boilerplate code, and git infrastructure for template synchronization. By starting with the nf-core template, you benefit from automated testing, consistent structure, and easier integration with the nf-core community.

Follow these steps to create a new nf-core pipeline, from initial proposal through to development. You will start by creating a new pipeline locally and working with it on your own GitHub account. Once you have a working version, you can request to move the repository to the nf-core GitHub organisation for final development, review, and release.

:::note{title="Prerequisites"}
You will need the following to get started:

- [Nextflow version 21.04.0 or later](../../get_started/environment_setup/nextflow.md)
- [nf-core/tools version 2.7 or later](../../get_started/environment_setup/nf-core-tools.md)
- [nf-test](https://www.nf-test.com/) for testing
- [pre-commit](https://pre-commit.com/) for code quality checks
- A container engine ([Docker](../../get_started/environment_setup/software-dependencies.md#docker), [Singularity](../../get_started/environment_setup/software-dependencies.md#singularity), or [Conda](../../get_started/environment_setup/software-dependencies.md#condamamba))
- A [GitHub account](https://github.com/signup)

:::

:::tip
To contribute features, bug fixes, or improvements to existing nf-core pipelines, see [Contributing to existing pipelines](/docs/contributors/contribute-existing-pipelines).
:::

## Request a new pipeline

To propose a new nf-core pipeline:

1. Create a proposal on the [nf-core/proposals](https://github.com/nf-core/proposals) repository.

    - Use the dedicated issue template for new pipeline proposals.
    - Fill out the form with key information about the pipeline you want to write.

1. Wait for the proposal to be reviewed.

    - The proposal will be discussed and checked for uniqueness to ensure minimal overlap with existing pipelines.
    - Once accepted, the proposal will be added to the [new pipeline tracking board](https://github.com/orgs/nf-core/projects/104) on GitHub.
    - The core team will create a Slack channel for your pipeline.

1. Join the following Slack channels:

    - [#pipeline-maintainers](https://nfcore.slack.com/channels/pipeline-maintainers) for major change announcements and general pipeline development discussion.
    - [#release-review-trading](https://nfcore.slack.com/channels/release-review-trading) to coordinate the two reviews required for your first release.

## Create the pipeline from the template

To create a new pipeline from the nf-core template:

1. Create the pipeline using the nf-core command-line tool:

    ```bash
    nf-core pipelines create
    ```

    - See the [nf-core tools documentation](/docs/nf-core-tools/pipelines/create) for detailed instructions.
    - This command creates the correct file structure, boilerplate code, and git infrastructure for template synchronization.

1. Switch to the `dev` branch:

    ```bash
    git checkout dev
    ```

    - All development should happen on `dev` or on other branches that get merged into `dev`.

:::tip
If you have an existing Nextflow pipeline, start fresh with the nf-core template and copy your code into the relevant places. This approach is usually easier than converting an existing pipeline. Contact the core team on Slack if you need guidance.
:::

## Push to GitHub

To push your new pipeline to GitHub:

1. Create an empty repository on GitHub under your personal account.

    - Go to GitHub and select **+** then **New Repository**.
    - Do not initialise the repository with any files (`README`, `LICENSE`, or `.gitignore`). The nf-core template already includes these files.

1. Add the GitHub repository as a remote to your local git repository:

    ```bash
    git remote add origin https://github.com/<YOUR_USERNAME>/<YOUR_REPOSITORY>.git
    ```

1. Push all branches to the remote GitHub repository:

    ```bash
    git push --all origin
    ```

    - The `nf-core pipelines create` command generates three standard nf-core branches (`master`, `dev`, and `TEMPLATE`) with an initial commit shared between them.
    - This git structure is required for automatic template synchronisation.

You should now see the nf-core template and branches in the GitHub web interface.

:::note
Repositories are generally not created within the nf-core organisation initially to avoid having empty pipelines listed if development takes longer than expected.
:::

## Develop your pipeline

To develop your pipeline:

1. Work on the `dev` branch following standard git development practices.

1. Run the nf-core linting tool regularly to ensure your pipeline meets nf-core standards:

    ```bash
    nf-core pipelines lint
    ```

    - See the [nf-core tools documentation](/docs/nf-core-tools/pipelines/lint) for more information.
    - GitHub Actions also runs this automatically, and you will receive notifications if tests fail.

1. Test the pipeline locally with test data to verify:

    - The pipeline executes successfully.
    - Outputs are generated correctly.
    - No errors or warnings occur.

1. Use the `debug` profile when testing to get additional information:

    ```bash
    nextflow run <pipeline_name> -profile debug
    ```

    - This enables warnings about process selectors, shows additional debug output, and disables cleanup.

## Set up test data

Test data is essential for validating your pipeline. nf-core uses the dedicated [nf-core/test-datasets](https://github.com/nf-core/test-datasets/) repository to store test data separately from pipeline code. This prevents repository bloat and ensures test data remains small and fast to execute.

:::note
GitHub Actions automatically run tests whenever someone updates the pipeline code. nf-core uses nf-test to check pipeline results, and regular test runs catch syntax errors and other problems that cause Nextflow to exit with an error.
:::

### Add test data to nf-core/test-datasets

To add your test data to the nf-core/test-datasets repository:

1. Fork the [nf-core/test-datasets](https://github.com/nf-core/test-datasets/) repository to your personal account.

1. Clone the repository and create a new branch for your pipeline:

    ```bash
    git clone https://github.com/<YOUR_USERNAME>/test-datasets.git
    cd test-datasets
    git checkout -b <pipeline_name>
    ```

1. Add your test data files.

    - Test data files must be very small due to GitHub file size limits and GitHub Actions timeout constraints.
    - Use minimal reference data such as PhiX, yeast, or part of a chromosome.
    - Aggressively subsample input data to be as small as possible while still testing pipeline functionality.

1. Ask for advice in the [#test-data](https://nfcore.slack.com/channels/test-data) Slack channel before adding test data.

1. Push the new files to GitHub:

    ```bash
    git add .
    git commit -m "Add test data for <pipeline_name>"
    git push --set-upstream origin <pipeline_name>
    ```

1. Create a pull request against the main [nf-core/test-datasets](https://github.com/nf-core/test-datasets/) repository.

    - First, create a new branch with your pipeline name using the GitHub web interface on the [nf-core/test-datasets repository](https://github.com/nf-core/test-datasets/).
    - Then open a pull request and select this branch as the target.

### Configure the test profile

To configure the test profile in your pipeline:

1. Edit the `conf/test.config` file.

    - Replace the example URLs with URLs to your test data files in the nf-core/test-datasets repository.
    - View files on GitHub and select **Raw** to get the direct URL.

1. Add required parameters to ensure the pipeline runs with minimal additional flags.

    - The `test` profile can be combined with other profiles such as `docker` or `conda`.
    - Do not specify a hardware environment in the test configuration.

1. Run the pipeline with the test profile to verify it works:

    ```bash
    nextflow run <pipeline_name> -profile test,docker --outdir <OUTDIR>
    ```

1. Update the `.github/workflows/` YAML files if you modify the test command.

### Test data requirements

All nf-core pipelines must meet these test data requirements:

- The `test` profile must exist and work.
- The test should be as comprehensive as possible.
- The test should run as much of the pipeline as possible.
- The test dataset should be as small as possible while still testing pipeline functionality.

See the [test data specifications](/docs/specifications/test-data/general) for complete requirements.

### Add full test data

You will need a full-size test dataset for the `test_full` profile:

1. Request upload to nf-core AWS S3 test buckets in the [#request-core](https://nfcore.slack.com/archives/C09H6NYHR9T) Slack channel.

    - This dataset produces realistic output when executed on each release.

## Move your pipeline to nf-core

Once your pipeline is written and tests pass, you can request to add your pipeline to the nf-core organisation.

To move your pipeline to the nf-core organisation:

1. Announce on the [#request-core](https://nfcore.slack.com/archives/C09H6NYHR9T) Slack channel that you need a core team member to move your repository.

1. Transfer repository ownership to nf-core:

    - Go to your repository settings.
    - Under the **General** page, find the **Danger Zone** section.
    - Select **Transfer Ownership** and transfer to the nf-core organisation.

1. Fork the transferred repository back to your account to continue development.

:::tip
You must be part of the nf-core organisation before transferring. If you are not a member, add a core team member to your repository and ask them to perform the transfer.
:::

:::note
Repositories should be transferred instead of forked to nf-core. Forking causes pull requests to default to the original user's fork rather than the nf-core repository, which can only be fixed by requesting manual detachment from GitHub support.
:::

### Configure branches

To set up the branch structure:

1. Ensure your repository has `dev` and `master` branches.

    - The `master` branch contains the latest stable release code.
    - The `dev` branch contains the latest development code.

1. Set `dev` as the default branch before the first release.

    - This ensures users run the latest development code by default.
    - After the first release, set the default branch to `master` so users run the latest stable release by default.

### Configure repository settings

Configure the following repository settings on GitHub:

- Add a description, the [https://nf-co.re](https://nf-co.re) URL, and keywords.
- Enable **Issues** and disable **Wiki** and **Projects**.
- Protect the `master` branch to require review and passing tests.
- Set write permissions for `nf-core/all` and admin permissions for `nf-core/admin`.

Verify these settings using the nf-core [repository health web page](https://nf-co.re/pipeline_health). This page reports the status of various checks and can fix errors automatically via the GitHub API.

:::note
When working with the nf-core repository, tests for pull requests against the `master` branch will fail because `master` should only contain code from the last release. Use the `dev` branch for new work and make all pull requests against `dev`.
:::

## Make the first release

Once your pipeline is written, tests pass, and the repository is in the nf-core organisation, you can prepare for your first release.

### Reset the default branch

When the code is stable and ready for a release, set the default branch to `master` on GitHub.

### Bump the version

To prepare for release:

1. Update the version number on the `dev` branch:

    ```bash
    nf-core pipelines bump-version --nextflow <new_version>
    ```

    - The first release should be version `1.0.0` (pre-releases are discouraged).
    - During development, use numeric versions with `dev` at the end (for example, `0.0.0dev`).
    - Release versions must be numeric and follow [Semantic Versioning](https://semver.org/).

:::tip
The version number exists in many places throughout the codebase. Use `nf-core pipelines bump-version` to update all occurrences consistently rather than editing manually.
:::

### Lint the pipeline

Run the linting tool with release checks to ensure your code meets nf-core standards:

```bash
nf-core pipelines lint --release
```

### Request review

Your pipeline requires review by the nf-core community before release.

To get your pipeline reviewed:

1. Create a pull request from the `dev` branch to `master` on the nf-core repository.

    - This is a special case where all tests should pass.
    - Ensure all automated tests complete successfully.

1. Request two reviews, with at least one from the core or maintainers team.

    - Use the [#release-review-trading](https://nfcore.slack.com/archives/C08K66XCZSL) Slack channel to coordinate reviews.

The review process depends on your `master` branch state:

- **Clean `master` branch**: If `master` contains only the initial commit, the pull request will represent all changes for review.
- **Existing commits**: The core team may merge your pull request and create a pseudo-pull request against the first commit to show all changes.

Reviewers will provide feedback on:

- README clarity and completeness
- Documentation quality, particularly output files and parameters
- Pipeline code structure and functionality

:::tip
While waiting for review, consider reviewing another community member's release pull request. This helps you understand the review process and contributes to the community. See the documentation pinned in [#release-review-trading](https://nfcore.slack.com/archives/C08K66XCZSL) for instructions.
:::

### Complete the release

Once reviewers approve your pull request:

1. Add a changelog entry describing the pipeline functionality at release.

    - Describe the general features and capabilities.
    - Tag contributors with their GitHub handles so their icons appear on the release page.

1. Follow the [pipeline release checklist](/docs/contributors/pipeline-release) to complete the release process.

The nf-core website and tools will automatically detect and display your new release.

:::note{title="Publications"}
If you publish research based on the pipeline, [give credit to nf-core](/docs/specifications/pipelines/recommendations/publication_credit).
:::

### Share your pipeline

Share your new pipeline with the community:

- Post in the [#bytesize_suggestion](https://nfcore.slack.com/archives/C081F8J2X8R) Slack channel to arrange a 15-minute introductory talk about your pipeline.

## Display your pipeline on the nf-core website

Verify that your pipeline appears on the nf-core website:

1. Check the [list of ignored repositories](https://github.com/nf-core/website/blob/main/sites/main-site/src/config/ignored_repos.yaml).

1. If your pipeline appears in this list, create a pull request to remove it.

1. Request approval from a core team member.

## Subsequent releases

After your first release, continue development following these practices:

1. Make pull requests against the `dev` branch on the nf-core repository.

    - All pull requests to `dev` should be reviewed before merging.

1. When ready for a new release:

    - Update the version number using `nf-core pipelines bump-version`.
    - Create a pull request against `master`.
    - If tests pass, the pull request can be merged and a new release made.

:::note
The `master` branch should only contain the commit from the latest release. This is important because the commit ID is used to determine whether a pipeline is up to date.
:::

## Additional considerations

### Developing pipelines privately

While nf-core encourages open development, private development is possible with some considerations:

- Private pipelines receive the same review process as public pipelines.
- Follow the same development steps but use a private repository.
- Make the repository public before release, or request reviews in the `#request-review` Slack channel if you cannot make it public until after release.

:::warning
Private pipelines receive no priority over publicly developed pipelines. If another similar pipeline is released first, your pipeline may not be accepted into nf-core. In this case, you may need to release independently or merge with the existing pipeline.
:::

If you plan to develop a pipeline privately and release through nf-core, contact the core team on Slack to discuss your strategy.
