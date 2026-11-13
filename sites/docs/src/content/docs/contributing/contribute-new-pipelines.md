---
title: Contributing a new pipeline
subtitle: Contribute a new pipeline
shortTitle: New pipelines
---

nf-core pipelines use a standardized template that provides file structure, boilerplate code, and git infrastructure for template synchronization.
By starting with the nf-core template, you benefit from automated testing, consistent structure, and easier integration with the nf-core community.

Follow these steps to create a new nf-core pipeline, from initial proposal through to development.
Start by proposing your pipeline to the community, before you write any code.
After your proposal is accepted, you create the pipeline from the template and develop it on your own GitHub account.
When the pipeline is ready, you request to move the repository to the nf-core GitHub organisation for final development, review, and release.

:::danger{title="Propose your pipeline before you build it"}
nf-core is a community of builders, not a registry of finished pipelines.

If you plan to contribute your pipeline to nf-core, discuss your idea and open a proposal **before you start development**.

See [Propose before you build](#propose-before-you-build).
:::

To contribute features, bug fixes, or improvements to existing nf-core pipelines, see [Contributing to existing pipelines](./contribute-existing-pipelines).

:::note{title="Prerequisites"}
You will need the following to get started:

- [Nextflow version 21.04.0 or later](../get_started/environment_setup/nextflow.md)
- [nf-core/tools version 2.7 or later](../get_started/environment_setup/nf-core-tools.md)
- [nf-test](https://www.nf-test.com/) for testing (included with nf-core tools)
- [prek](https://prek.j178.dev/) for code quality checks (included with nf-core tools)
- A container engine ([Docker](../get_started/environment_setup/software-dependencies#docker), [Singularity](../get_started/environment_setup/software-dependencies#singularity), or [Conda](../get_started/environment_setup/software-dependencies#condamamba))
- A [GitHub account](https://github.com/signup)

:::

## Propose before you build

Not every pipeline built with nf-core tooling needs to become an nf-core pipeline.
You can use the nf-core template, tools, and components for your own pipelines without contributing them to nf-core.
See [External use of nf-core resources](../developing/external-use).

If you plan to contribute your pipeline to nf-core, involve the community from the start.
[Discuss your idea and propose your pipeline before you start work on it](../specifications/pipelines/recommendations/build_with_community), not after the pipeline is complete.
Proposals are public so that you can find collaborators.

Propose early to:

- Avoid duplicated work.
  nf-core allows only [one pipeline for a given data and analysis type](../specifications/pipelines/requirements/workflow_specificity).
  If a pipeline with a similar purpose already exists or is proposed before yours, you will be asked to contribute to it instead - even if your pipeline is more complete or comprehensive.
- Get feedback while changes are cheap.
  Changes needed to meet the [nf-core specifications](../specifications/pipelines/overview) are easier early in development than after the pipeline is complete.
- Find collaborators who can help develop, test, and maintain the pipeline.
- Confirm that the pipeline is suitable for nf-core before you put time into development.

If you are unsure whether your idea fits, ask on the [#help](https://nfcore.slack.com/channels/help) Slack channel before you open a proposal.

### Ownership and publication

nf-core pipelines are [community owned](../specifications/pipelines/requirements/community_owned).
By adding a pipeline to nf-core, you agree to the following:

- You work with other contributors.
- The community makes decisions about the pipeline by consensus, with the pipeline maintainers as the guiding figures, but with final say from the core team and steering committee if needed.
- Other people may add features and make releases.

Community ownership means the pipeline development is not exclusive single person or group, even if it is maintained by just one person.
As people's time commitments and funding change, maintenance can shift between groups and the pipeline stays up to date.

Community ownership does not stop you from publishing your work.
The primary developers of a pipeline are explicitly permitted to publish it, and nf-core encourages publications and presentations about its pipelines.
If your pipeline is part of a paper or thesis, you can still propose it publicly before the work is published.
See [Publication credit](../specifications/pipelines/recommendations/publication_credit) for how to acknowledge nf-core and discuss co-authorship with other contributors.

### Contribute a pre-existing pipeline

If you have a Nextflow pipeline built outside nf-core, [start the discussion](https://github.com/nf-core/proposals/issues) as early as possible.

nf-core can accept pre-existing pipelines, with the following caveats:

- Your pipeline must already use the nf-core template.
- You may be asked to contribute to an existing pipeline instead if one with an identical or similar purpose exists, even if your pipeline is more mature at the time of submission.
- You may need to make major (possibly breaking) changes to your pipeline to meet [the guidelines](/docs/guidelines/pipelines/overview).
- You may need to change the name of the pipeline to follow the [pipeline specifications](../specifications/pipelines/requirements/workflow_name).

A rejected proposal is not a statement on the quality of the pipeline.
Collaboration is a core principle of nf-core, and the community promotes larger teams around a single pipeline.

:::info
If your proposal is turned down, you can still contribute to existing pipelines, and use the nf-core template, tools, and components for your own pipeline.
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

:::warning
Create your pipeline only after you have [proposed it](#request-a-new-pipeline) and discussed it with the community.
:::

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
If you have an existing Nextflow pipeline, start fresh with the nf-core template and copy your code into the relevant places.
This approach is usually easier than converting an existing pipeline.
Contact the core team on Slack if you need guidance.
:::

## Push to GitHub

To push your new pipeline to GitHub:

1. Create an empty repository on GitHub under your personal account.
   - Go to GitHub and select **+** then **New Repository**.
   - Do not initialise the repository with any files (`README`, `LICENSE`, or `.gitignore`).
     The nf-core template already includes these files.

1. Add the GitHub repository as a remote to your local git repository:

   ```bash
   git remote add origin https://github.com/<username>/<repository>.git
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
   nextflow run <pipeline-name> -profile debug
   ```

   - This enables warnings about process selectors, shows additional debug output, and disables cleanup.

## Set up test data

Test data is essential for validating your pipeline.
nf-core uses the dedicated [nf-core/test-datasets](https://github.com/nf-core/test-datasets/) repository to store test data separately from pipeline code.
This prevents repository bloat and ensures test data remains small and fast to execute.

:::note
GitHub Actions automatically run tests whenever someone updates the pipeline code.
nf-core uses nf-test to check pipeline results, and regular test runs catch syntax errors and other problems that cause Nextflow to exit with an error.
:::

### Add test data to nf-core/test-datasets

To add your test data to the nf-core/test-datasets repository:

1. Fork the [nf-core/test-datasets](https://github.com/nf-core/test-datasets/) repository to your personal account.

1. Clone the repository and create a new branch for your pipeline:

   ```bash
   git clone https://github.com/<username>/test-datasets.git
   cd test-datasets
   git checkout -b <pipeline-name>
   ```

1. Add your test data files.
   - Test data files must be very small due to GitHub file size limits and GitHub Actions timeout constraints.
   - Use minimal reference data such as PhiX, yeast, or part of a chromosome.
   - Aggressively subsample input data to be as small as possible while still testing pipeline functionality.

1. Ask for advice in the [#test-data](https://nfcore.slack.com/channels/test-data) Slack channel before adding test data.

1. Push the new files to GitHub:

   ```bash
   git add .
   git commit -m "Add test data for <pipeline-name>"
   git push --set-upstream origin <pipeline-name>
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
   nextflow run <pipeline-name> -profile test,docker --outdir <outdir>
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

You will need a full-size test dataset for the `test_full` profile.
This dataset produces realistic output when executed on each release and gets rendered on the "Results" tab of your pipeline's page on nf-co.re.
Whenever possible, resources such as datasets or databases required here should be specified as stable, secure links to the original source.
In rare cases where the data is only available via insecure http, behind a login page, or otherwise difficult to download programatically, contact the core team at [#request-core](https://nfcore.slack.com/channels/request-core.

For more information about running the tests using the `test_full` profile, watch [Bytesize 19](https://nf-co.re/events/2021/bytesize-19-aws-megatests).

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
You must be part of the nf-core organisation before transferring.
If you are not a member, add a core team member to your repository and ask them to perform the transfer.
:::

:::note
Repositories should be transferred instead of forked to nf-core.
Forking causes pull requests to default to the original user's fork rather than the nf-core repository, which can only be fixed by requesting manual detachment from GitHub support.
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

- Add a description, the [https://nf-co.re](https://nf-co.re) URL and lots of keywords
- Enable issues, disable **Wiki** and **Projects**
- Enable **Always suggest updating pull request branches**
- Enable **Automatically delete head branches**
- Add the **contributors** team with **Write** access and the **core** team with **Admin** access
- Set `master`/`main`, `dev` and `TEMPLATE` branches as protected, so that they require review and passing tests (via rulesets)

The full step-by-step checklist for these settings is in the [core team checklist](/docs/community/governance/core-team#pipeline-repository-settings).

Verify these settings using the nf-core [repository health web page](https://nf-co.re/pipeline_health).
This page reports the status of various checks and can fix errors automatically via the GitHub API.

:::note
When working with the nf-core repository, tests for pull requests against the `master` branch will fail because `master` should only contain code from the last release.
Use the `dev` branch for new work and make all pull requests against `dev`.
:::

## Make the first release

Once your pipeline is written, tests pass, and the repository is in the nf-core organisation, you can prepare for your first release.

### Reset the default branch

When the code is stable and ready for a release, set the default branch to `master` on GitHub.

### Bump the version

To prepare for release:

1. Update the version number on the `dev` branch:

   ```bash placeholder="<new-version"
   nf-core pipelines bump-version <new-version>
   ```

   - The first release should be version `1.0.0` (pre-releases are discouraged).
   - During development, use numeric versions with `dev` at the end (for example, `0.0.0dev`).
   - Release versions must be numeric and follow [Semantic Versioning](https://semver.org/).

:::tip
The version number exists in many places throughout the codebase.
Use `nf-core pipelines bump-version` to update all occurrences consistently rather than editing manually.
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
While waiting for review, consider reviewing another community member's release pull request.
This helps you understand the review process and contributes to the community.
See the documentation pinned in [#release-review-trading](https://nfcore.slack.com/archives/C08K66XCZSL) for instructions.
:::

### Complete the release

Once reviewers approve your pull request:

1. Ask a core team member on [#request-core](https://nfcore.slack.com/channels/request-core) to remove the **"Block releases"** GitHub repository ruleset.
   - All unreleased nf-core pipelines have this ruleset applied, which prevents tag creation until the pipeline is ready for its first release.

1. Ask a core team member on [#request-core](https://nfcore.slack.com/channels/request-core) to activate Zenodo archiving of the GitHub repository.

1. Add a changelog entry describing the pipeline functionality at release.
   - Describe the general features and capabilities.
   - Tag contributors with their GitHub handles so their icons appear on the release page.

1. Follow the [release procedure](./release-procedure) to complete the release process.

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

See the [release procedure](./release-procedure) for the full checklist.

:::note
The `master` branch should only contain the commit from the latest release.
This is important because the commit ID is used to determine whether a pipeline is up to date.
:::
