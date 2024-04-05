---
title: Adding a new pipeline
subtitle: Follow this walkthrough to add a new pipeline to nf-core.
---

# Before you start

So, you want to add a new pipeline to nf-core - brilliant!
Before you start typing, check that you're happy with the following points:

- You're familiar with nf-core and nextflow (see our [introduction docs](/docs/usage/introduction.md)).
- You're used to working with `git` and [GitHub](https://github.com)
  (see a [nice tutorial here](https://blog.scottlowe.org/2015/01/27/using-fork-branch-git-workflow/))
- The workflow you're thinking of meets the [nf-core guidelines](https://nf-co.re/docs/contributing/guidelines).

The main steps involved in adding a new nf-core pipeline covered below are:

1. [Joining the community](#join-the-community)
2. [Creating a pipeline](#create-a-pipeline-from-the-template)
3. [Running with test data](#running-with-test-data)
4. [Adding to the nf-core organisation](#adding-your-pipeline-to-the-nf-core-organisation)
5. [Making your first release](#making-the-first-release)
6. [Updates and new releases](#subsequent-releases)

# Join the community

At its heart, nf-core is a community - to add a pipeline you need to be part of that community!
Please join us on [Slack](https://nf-co.re/join/slack), and ask to be added to the GitHub association through the [#github-invitations](https://nfcore.slack.com/channels/github-invitations) channel. If you feel like it, you can go to the [#say-hello](https://nfcore.slack.com/channels/say-hello) channel and introduce yourself to the rest of the community.

:::warning
It's good to introduce your idea early on so that it can be discussed, before you spend lots of time coding.
:::

The [nf-core guidelines](/docs/contributing/guidelines) state that no two pipelines should overlap too much
in their purpose and results. There may be an existing pipeline that can be extended to give the
functionality that you are looking for, or there could be another group working on a similar to the
pipeline to the one you're planning.

To avoid problems at a later date, please come and discuss your plans with the nf-core community as early
as possible. Ideally before you make a start on your pipeline!

:::info
Not all pipelines are suitable for inclusion in the main nf-core community (eg. bespoke or proprietary workflows). However, we hope that you may still wish to use the nf-core template and/or use components of nf-core code. All nf-core code is under a MIT license and where possible we have endeavoured to make the tools work with any Nextflow pipeline. If this is the case for you, please see the [unofficial pipelines tutorial](/docs/contributing/tutorials/unofficial_pipelines.md) for more details.
:::

All nf-core discussion happens on the nf-core Slack, which you can join here:
[https://nf-co.re/join](https://nf-co.re/join)

These topics are specifically discussed in the `#new-pipelines` channel:
[https://nfcore.slack.com/channels/new-pipelines](https://nfcore.slack.com/channels/new-pipelines)
