---
title: Start contributing
subtitle: How to start contributing to nf-core and contribution etiquette
shortTitle: Overview
weight: 1
---

Now you have run your first pipeline, you may wish to start contributing!

There are many ways to contribute to nf-core.
This can be through pipelines, pipeline components, documentation, tooling, or even just feature requests or bug reports.

If this is your first time contributing, please read this page to orient yourself on how contributing works within our community.

## Code of conduct

The [code of conduct](https://nf-co.re/code_of_conduct) is a core document for ensuring everyone gets along well within our diverse, distributed group.
Before you start, take a moment to review this document.

## Communication

Our primary communication platform is our [Slack workspace](https://nfcore.slack.com/).
We also communicate using [Issues](#issues) and [Pull Requests](#pull-requests) on our GitHub repositories.

Popular common repositories include:

- [nf-core/modules](https://github.com/nf-core/modules): Shared Nextflow modules and subworkflows
- [nf-core/tools](https://github.com/nf-core/tools): nf-core toolkit, including module, subworkflow, and pipeline templates
- [nf-core/website](https://github.com/nf-core/website): nf-core documentation
- [nf-core/testdatasets](https://github.com/nf-core/testdatasets): Shared test data files for unit testing
- [nf-core/configs](https://github.com/nf-core/configs): Shared Nextflow configuration files for HPC and cloud infrastructure

:::tip
Every pipeline has a dedicated channel on [Slack](#joining-the-slack-and-github-organisations).
Search for the pipeline name in the Channel browser.
:::

## Joining the Slack and GitHub organisations

To join the Slack workspace, see our [Join nf-core](https://nf-co.re/join) page for instructions.

After joining Slack, we recommend you join our GitHub organisation.
The GitHub organisation is open to everyone and offers a range of benefits.
These include not having to wait for someone to approve CI tests, and reviewing and merging rights (on most repositories).

:::tip
You need to join our Slack workspace first before joining the GitHub organisation!
:::

## Asking for help

If you have questions and don't know where to start, ask in our dedicated Slack help channels:

- [#nostupidquestions](https://nfcore.slack.com/archives/C043FMKUNLB)
- [#help](https://nfcore.slack.com/archives/CE6SDBX2A)

## Issues

When you want to report a bug, or request a new feature, make an Issue on the relevant GitHub repository.

Ensure issues are descriptive:

- For bug reports
  - Include log files (`.nextflow.log` is a common one for pipelines)
  - Provide the _whole_ error message
  - Give enough information for a developer to reproduce the bug, such as the command used for execution
- For feature requests
  - Include links to other resources (for example, for adding a new tool to a pipeline, a link to the tools code repository)
  - Include links to existing discussions, such as on Slack

## Contributing

### Forks and branches

In most cases, we ask for contributors to develop on forks of original repositories.

When on a GitHub repository, press the **Fork** button in the top right set of buttons on the page.
This will make a copy of the repository under your own user name for you to develop on.

:::tip
We recommend naming the forked repository to `nf-core-<original name>`, to distinguish from the original repository.
For example, `nf-core/rnaseq` would become `<your username>/nf-core-rnaseq`.
:::

Once you are on your fork, please make a `git` branch to make your changes on.

### Pull requests

Once ready, open a pull request from your fork to the original repository.

Everyone is welcome to open a pull request on our GitHub repositories.
As a general rule, all pull requests must be [reviewed](#reviewing) by another community member.

Be patient when waiting for a review.
nf-core is driven by volunteers and can only review when they have time.

:::tip
Here are some tips to get faster reviews:

1. **Keep your pull requests small**: one or two modified files, or less than 100 lines of code, will more likely be reviewed than 100 files with 1000s of lines of code changes.
2. **Request reviews on Slack**: use the [#request-review](https://nfcore.slack.com/archives/request-review) channel on Slack to get visibility of your request
3. **Trade reviews**: if you request a review, review someone else's - they will likely reciprocate! Track other requests on the #request-review channel.

:::

### Contributing to shared repositories

To contribute to specific nf-core shared repositories, please see the specific guidance on the relevant pages:

- [Contribute components](../contributing/contribute-components)
- [Contribute to existing pipelines](../contributing/contribute-existing-pipelines)
- [Contribute to new pipelines](../contributing/contribute-new-pipelines)

### Contributing to new pipelines

If you wish to contribute a new pipeline to nf-core, _before you start working_, refer to [Contribute new pipelines](../contributing/contribute-new-pipelines.md).

:::note
nf-core is not a registry of any Nextflow pipelines!
Collaboration is a core philosophy of nf-core to ensure robust and long-living pipelines.
:::

## Reviewing

Reviews and comments can be technical (code), conceptual (for example, purpose/background/approach/domain specific), or both.

### Who can review

Everyone is welcome to review or leave a comment on an nf-core pull request.
You do not need to be in the nf-core community to give feedback.
You also do not need to be in the nf-core GitHub organisation (although this is recommended).

However, if you review a PR about a 'conceptual' question, a second technical review by another community member is required.
Reviews from community members who are not directly within your development team are recommended.

### How to review

Be polite and constructive in your comments.
Our community has a diverse range of backgrounds, cultures, and expertise.

Help make a smooth reviewing process by:

- Ensuring replies follow the [code of conduct](#code-of-conduct).
- Making _one_ review with multiple comments using the **Start review** button
  - Avoid the **Add single comment** button on the **GitHub Files Changed** tab
- Leaving descriptive comments with useful suggestions, not vague observations
- Making [direct code suggestions](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/reviewing-proposed-changes-in-a-pull-request#starting-a-review) on specific lines, avoid just describing your suggestions
- Preferring [**Comments** review types](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/reviewing-proposed-changes-in-a-pull-request#submitting-your-review) rather than 'Request changes'
  - **Request changes** blocks merging until the specific reviewer approves the PR
  - **Comments** allows other community members to re-review and finish the review on your behalf
  - Trying out [Conventional comments](https://conventionalcomments.org/) for effective feedback categorisation

General reviewing checklists for common nf-core processes can be seen here:

- [Reviewing nf-core components](../contributing/reviewing-pull-requests/review_checklists/component)
- [Reviewing pipeline releases](../contributing/reviewing-pull-requests/review_checklists/pipeline-releases)
- [Reviewing nf-core/tools](../contributing/reviewing-pull-requests/review_checklists/nf-core-tools.md)

## Community structure and contact

See the [Governance](https://nf-co.re/governance#governance) page for an overview of the communities' organisational structure.

As a general rule, always start by asking the entire community first!
Answers to most questions can be obtained from the entire nf-core 'hive-mind'.

If you get stuck on something that needs more expert advice or a decision, there are specific teams that may be able to help:

- [Maintainers team](https://nf-co.re/governance#maintainers) for issue with modules, pipelines, or test configs
- [Outreach team](https://nf-co.re/governance#outreach) for questions about nf-core events or the bytesize series
- [Infrastructure team](https://nf-co.re/governance#infrastructure) for problems with the nf-core template, tooling, GitHub actions, or AWS
- [Core team](https://nf-co.re/governance#core-team) for critical decisions that may affect the entire community, or for things requiring passwords
- [Safety team](https://nf-co.re/governance#safety) for inter-personal conflict with another member or member(s) of the community
