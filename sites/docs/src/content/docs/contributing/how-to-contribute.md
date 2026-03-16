---
title: How to contribute
subtitle: How to start contributing to nf-core and contribution etiquette
shortTitle: How to contribute
weight: 2
---

There are many ways to contribute to nf-core.
This can be through pipelines, pipeline components, documentation, tooling, or even just feature requests or bug reports.

This page provides a high-level overview and best practices to help you contribute to nf-core quickly.

## Communication

[Slack](https://nfcore.slack.com/) and GitHub [Issues](#issues) and [Pull Requests](#pull-requests) are the primary ways to communicate with others in the nf-core community.

:::tip
Every pipeline has a dedicated channel on Slack.
You can search pipeline names in the channel browser.
:::

## nf-core Slack and GitHub

See [Join nf-core](https://nf-co.re/join) for instructions on how to join the nf-core Slack.

After joining Slack, [join the nf-core GitHub organization](https://nf-co.re/join#github).
The GitHub organization is open to everyone and gives you reviewing and merging rights on most repositories.

:::tip
Join the Slack workspace before joining the GitHub organization.
:::

## Asking for help

If you have questions and don't know where to start, ask in the help Slack channels:

- [#nostupidquestions](https://nfcore.slack.com/archives/C043FMKUNLB)
- [#help](https://nfcore.slack.com/archives/CE6SDBX2A)

## Issues

To report a bug or request a new feature, open an issue on the relevant GitHub repository.

### Bug reports

Include enough detail for someone to reproduce the issue:

- Include log files (`.nextflow.log` is a common one for pipelines)
- Provide the _whole_ error message
- Give enough information to reproduce the bug, such as the command used for execution

### Feature requests

Include enough context for the community to evaluate your request:

- Include links to other resources (for example, for adding a new tool to a pipeline, a link to the tool's code repository)
- Include links to existing discussions, such as on Slack

## Contributing

Contributions to nf-core are made through forks, branches, and pull requests on GitHub.

### Forks

When working on nf-core projects, develop on a fork of the repository.

To fork a repository, click **Fork** in the top-right corner.
This creates a copy of the repository under your username.

:::tip
We recommend naming the forked repository `nf-core-<original name>`, to distinguish it from the original repository.
For example, `nf-core/rnaseq` would become `<your username>/nf-core-rnaseq`.
:::

### Branches

On your fork, create a `git` branch for your changes.

### Pull requests

Everyone is welcome to open a pull request from their fork into the nf-core GitHub repositories.
As a general rule, another community member must [review](#reviewing) all pull requests.

Be patient when waiting for a review.
nf-core volunteers review pull requests in their spare time.

:::tip
Here are some tips to get faster reviews:

1. **Keep your pull requests small**: one or two modified files, or less than 100 lines of code, will more likely be reviewed than 100 files with 1000s of lines of code changes.
2. **Request reviews on Slack**: use the [#request-review](https://nfcore.slack.com/archives/request-review) channel on Slack to get visibility of your request.
3. **Trade reviews**: if you request a review, review someone else's — they will likely reciprocate! Track other requests on the #request-review channel.
   :::

For detailed guidance on contributing to specific nf-core repositories, see:

- [Contribute components](../contributing/contribute-components)
- [Contribute to existing pipelines](../contributing/contribute-existing-pipelines)
- [Contribute to new pipelines](../contributing/contribute-new-pipelines)

:::caution[Starting a new pipeline?]
nf-core is not a registry of Nextflow pipelines.
Collaboration is a core requirement.
Before writing any code, see [Contribute new pipelines](../contributing/contribute-new-pipelines.md).
:::

## Reviewing

Reviews and comments can be technical (code), conceptual (for example, purpose/background/approach/domain-specific), or both.

### Who can review

Everyone is welcome to review or leave a comment on an nf-core pull request.
You do not need to be in the nf-core community to give feedback.
You also do not need to be in the nf-core GitHub organization (although this is recommended).

However, if you review a PR about a conceptual question, another community member must also provide a technical review.
Seek reviews from community members outside your development team.

### How to review

Be polite and constructive in your comments.
Our community has a diverse range of backgrounds, cultures, and expertise.

Help make a smooth reviewing process by:

- Ensuring replies follow the [code of conduct](#code-of-conduct)
- Making _one_ review with multiple comments using the **Start review** button
  - Avoid the **Add single comment** button on the **GitHub Files Changed** tab
- Leaving descriptive comments with useful suggestions, not vague observations
- Making [direct code suggestions](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/reviewing-proposed-changes-in-a-pull-request#starting-a-review) on specific lines, avoid just describing your suggestions
- Preferring [**Comments** review types](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/reviewing-proposed-changes-in-a-pull-request#submitting-your-review) rather than **Request changes**
  - **Request changes** blocks merging until the specific reviewer approves the PR
  - **Comments** allows other community members to re-review and finish the review on your behalf
- Trying out [Conventional comments](https://conventionalcomments.org/) for effective feedback categorization

### Reviewing checklists

See the following for general reviewing checklists:

- [Reviewing nf-core components](../contributing/reviewing-pull-requests/review_checklists/component)
- [Reviewing pipeline releases](../contributing/reviewing-pull-requests/review_checklists/pipeline-releases)
- [Reviewing nf-core/tools](../contributing/reviewing-pull-requests/review_checklists/nf-core-tools.md)

## When everything else fails

As a general rule, ask the wider community first.
Collectively, the nf-core community can answer most questions.

:::tip[Need specialist help?]
If the community can't help, reach out to the relevant team:

- [Maintainers team](https://nf-co.re/governance#maintainers) for issues with modules, pipelines, or test configs
- [Outreach team](https://nf-co.re/governance#outreach) for questions about nf-core events or the bytesize series
- [Infrastructure team](https://nf-co.re/governance#infrastructure) for problems with the nf-core template, tooling, GitHub Actions, or AWS
- [Core team](https://nf-co.re/governance#core-team) for critical decisions affecting the entire community, or anything requiring passwords
- [Docs team](https://nf-co.re/governance#docs) for questions about nf-core documentation
- [Safety team](https://nf-co.re/governance#safety) for interpersonal conflict with another community member
  :::

See the [Governance](https://nf-co.re/governance#governance) page for an overview of the community's organizational structure.
