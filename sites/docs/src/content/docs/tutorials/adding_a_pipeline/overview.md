---
title: Overview
subtitle: Follow this walkthrough to add a new pipeline to nf-core.
type: tutorial
parentWeight: 10
weight: 1
---

The steps covered in this tutorial are:

1. Introduction (this page!)
1. [Create a pipeline](/docs/tutorials/adding_a_pipeline/creating_a_pipeline)
1. [Running with test data](/docs/tutorials/adding_a_pipeline/test_data)
1. [Adding to the nf-core organisation](/docs/tutorials/adding_a_pipeline/move_to_nf-core_org)
1. [Making your first release](/docs/tutorials/adding_a_pipeline/first_release)
1. [Next steps](/docs/tutorials/adding_a_pipeline/next_steps)

## Before you start

So, you want to add a new pipeline to nf-core - brilliant!

Before you start typing, check that you're happy with the following points:

- You're familiar with nf-core and Nextflow
  - See our [introduction docs](/docs/usage/getting_started/introduction.md))
- You're able to work with `git` and [GitHub](https://github.com)
    - See a [nice tutorial here](https://blog.scottlowe.org/2015/01/27/using-fork-branch-git-workflow/))
- You're familiar with the [nf-core guidelines](https://nf-co.re/docs/guidelines/pipelines/overview)

## Join the community

At its heart, nf-core is a community - to add a pipeline you need to be part of that community!

To join the community, join the [nf-core Slack](https://nf-co.re/join/slack) and ask to be added to the GitHub association on the [#github-invitations](https://nfcore.slack.com/channels/github-invitations) channel. If you feel like it, go to the [#say-hello](https://nfcore.slack.com/channels/say-hello) channel and introduce yourself to the rest of the community.

The [nf-core guidelines](/docs/guidelines/pipelines/overview) state that no two pipelines should overlap too much
in their purpose and results. There may be an existing pipeline that can be extended with the
functionality that you are looking for, or there could be another group working on a pipeline similar to the one you're planning.

To avoid potential pipeline conflicts, come and discuss your plans with the wider nf-core community as early
as possible. Ideally before you start writing your pipeline!

:::info
Not all pipelines are suitable for inclusion in the nf-core community, e.g., bespoke or proprietary pipelines. However, the nf-core template and components are available for all to use. All nf-core code is under a MIT license and, where possible, the tools work with any Nextflow pipeline. See the [unofficial pipelines guidelines](/docs/guidelines/external_use) for more details.
:::

All nf-core discussion happens on the nf-core Slack, which you can [join here](https://nf-co.re/join).

These topics are specifically discussed in the [`#new-pipelines` channel](https://nfcore.slack.com/channels/new-pipelines).
