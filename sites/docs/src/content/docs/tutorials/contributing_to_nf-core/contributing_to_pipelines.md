---
title: Contributing to an existing pipeline
subtitle: Follow this walkthrough to add features to an existing nf-core pipeline.
---

## Before you start

So, you want to add a new pipeline to nf-core - brilliant!
Before you start typing, check that you're happy with the following points:

- You're familiar with nf-core and Nextflow (see our [introduction docs](/docs/usage/getting_started/introduction.md)).
- You're used to working with `git` and [GitHub](https://github.com)
  (see a [nice tutorial here](https://blog.scottlowe.org/2015/01/27/using-fork-branch-git-workflow/))
- The workflow you're thinking of meets the [nf-core guidelines](https://nf-co.re/docs/contributing/guidelines).

The main steps involved in adding a new nf-core pipeline covered below are:

1. [Joining the community](#join-the-community)
2. [Contribution overview](#contribution-overview)
3. [Testing](#testing)
4. [Patching bugs](#patching-bugs)
5. [Getting help](#getting-help)
6. [Pipeline contribution conventions](#pipeline-contribution-conventions)

## Join the community

At its heart, nf-core is a community. To contribute to a pipeline, you should be part of that community!
Please request to join the [nf-core GitHub organisation](/join/#github))
and join the nf-core [Slack](/join/slack).

Each nf-core pipeline has its own Slack channel that can be found by searching the channel list.

## Contribution overview

:::warning
It's good to introduce your idea early on so that it can be discussed, before you spend lots of time coding.
:::

If you'd like to write some code for an nf-core pipeline, the standard workflow is as follows:

1. Check that there isn't already an issue about your idea for that pipeline to avoid duplicating work. If there isn't one already, please create one so that others know you're working on this
2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the pipeline repository to your GitHub account
3. Make the necessary changes / additions on `dev` branch by using [git checkout](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/committing-changes-to-a-pull-request-branch-created-from-a-fork) within your forked repository following [pipeline conventions](#pipeline-contribution-conventions)
4. Use `nf-core pipelines schema build` and add any new parameters to the pipeline JSON schema (requires [nf-core tools](https://github.com/nf-core/tools) >= 1.10).
5. Submit a Pull Request against the `dev` branch and wait for the code to be reviewed and merged. See [guidelines](/docs/guidelines/pull_request_review) for this.

If you're new to working with git, you can view the [GitHub pull requests documentation](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or other [excellent `git` resources](https://try.github.io/) to get started.

## Testing

You can optionally test your changes by running the pipeline locally. Then it is recommended to use the `debug` profile to
receive warnings about process selectors and other debug info. Example: `nextflow run . -profile debug,test,docker --outdir <OUTDIR>`.

When you create a pull request with changes, [GitHub Actions](https://github.com/features/actions) will run automatic tests.
Typically, pull-requests are only fully reviewed when these tests are passing, though of course, we can help out before then.

There are typically two types of tests that run:

#### Lint tests

`nf-core` has a [set of guidelines](https://nf-co.re/developers/guidelines) which all pipelines must adhere to.
To enforce these and ensure that all pipelines stay in sync, we have developed a helper tool that runs checks on the pipeline code. This is in the [nf-core/tools repository](https://github.com/nf-core/tools) and once installed can be run locally with the `nf-core pipelines lint <pipeline-directory>` command.

If any failures or warnings are encountered, please follow the listed URL for more documentation.

#### Pipeline tests

Each `nf-core` pipeline should be set up with a minimal set of test data.
`GitHub Actions` then runs the pipeline on this data to ensure that it exits successfully.
If there are any failures then the automated tests fail.
These tests are run both with the latest available version of `Nextflow` and also the minimum required version that is stated in the pipeline code.

## Patching bugs

:::warning
Only in the unlikely and regretful event of a release happening with a bug.
:::

- On your own fork, make a new branch `patch` based on `upstream/master`.
- Fix the bug, and bump version (X.Y.Z+1).
- A PR should be made on `master` from patch to directly resolve this particular bug.

## Getting help

For further information/help, please consult the usage documentation for the particular pipeline and don't hesitate to get in touch on the nf-core Slack `#<pipeline-name>` channel. If you are not already a member you can ([join our Slack here](https://nf-co.re/join/slack)).

## Pipeline contribution conventions

To make the nf-core pipeline code and processing logic more understandable for new contributors and to ensure quality, we semi-standardize the way the code and other contributions are written.

#### Adding a new step

If you wish to contribute a new step, please use the following coding standards:

1. Define the corresponding input channel into your new process from the expected previous process channel.
2. Write the process block (see below).
3. Define the output channel if needed (see below).
4. Add any new parameters to `nextflow.config` with a default (see below).
5. Add any new parameters to `nextflow_schema.json` with help text (via the `nf-core pipelines schema build` tool).
6. Add sanity checks and validation for all relevant parameters.
7. Perform local tests to validate that the new code works as expected.
8. If applicable, add a new test command in `.github/workflow/ci.yml`.
9. Update MultiQC config `assets/multiqc_config.yml` so relevant suffixes, file name clean up and module plots are in the appropriate order. If applicable, add a [MultiQC](https://multiqc.info/) module.
10. Add a description of the output files and if relevant any appropriate images from the MultiQC report to `docs/output.md`.

#### Default values

Parameters should be initialised / defined with default values in `nextflow.config` under the `params` scope.

Once there, use `nf-core pipelines schema build` to add to `nextflow_schema.json`.

#### Default processes resource requirements

Sensible defaults for process resource requirements (CPUs / memory / time) for a process should be defined in `conf/base.config`. These should generally be specified generic with `withLabel:` selectors so they can be shared across multiple processes/steps of the pipeline. An nf-core standard set of labels that should be followed where possible can be seen in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/conf/base.config), which has the default process as a single core-process, and then different levels of multi-core configurations for increasingly large memory requirements defined with standardised labels.

The process resources can be passed on to the tool dynamically within the process with the `${task.cpu}` and `${task.memory}` variables in the `script:` block.

#### Naming schemes

Please use the following naming schemes, to make it easy to understand what is going where.

- initial process channel: `ch_output_from_<process>`
- intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

#### Nextflow version bumping

If you are using a new feature from core Nextflow, you may bump the minimum required version of Nextflow in the pipeline with: `nf-core pipelines bump-version --nextflow . [min-nf-version]`

#### Images and figures

For overview images and other documents, we follow the nf-core [style guidelines and examples](https://nf-co.re/developers/design_guidelines).
