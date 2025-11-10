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

- You're familiar with nf-core and nextflow (see our [introduction docs](/docs/usage/getting_started/introduction.md)).
- You're used to working with `git` and [GitHub](https://github.com)
  (see a [nice tutorial here](https://blog.scottlowe.org/2015/01/27/using-fork-branch-git-workflow/))
- The workflow you're thinking of meets the [nf-core guidelines](https://nf-co.re/docs/guidelines/pipelines/overview).

## Join the community

At its heart, nf-core is a community - to add a pipeline you need to be part of that community!
Please join us on [Slack](https://nf-co.re/join/slack), and ask to be added to the GitHub association through the [#github-invitations](https://nfcore.slack.com/channels/github-invitations) channel. If you feel like it, you can go to the [#say-hello](https://nfcore.slack.com/channels/say-hello) channel and introduce yourself to the rest of the community.

:::warning
It's good to introduce your idea early on so that it can be discussed, before you spend lots of time coding.
:::

The [nf-core guidelines](/docs/guidelines/pipelines/overview) state that no two pipelines should overlap too much
in their purpose and results. There may be an existing pipeline that can be extended to give the
functionality that you are looking for, or there could be another group working on a similar to the
pipeline to the one you're planning.

To avoid problems at a later date, please come and discuss your plans with the nf-core community as early
as possible. Ideally before you make a start on your pipeline!

:::info
Not all pipelines are suitable for inclusion in the main nf-core community (eg. bespoke or proprietary workflows). However, we hope that you may still wish to use the nf-core template and/or use components of nf-core code. All nf-core code is under a MIT license and where possible we have endeavoured to make the tools work with any Nextflow pipeline. If this is the case for you, please see the [unofficial pipelines guidelines](/docs/guidelines/external_use) for more details.
:::

## Pre-existing pipelines

In some cases, you might have a general-purpose Nextflow pipeline you built outside of nf-core that you believe would be a valuable addition to the community. We welcome this type of contribution. However, you should be aware of several caveats in this scenario:

- If there is already a released or in-development nf-core pipeline with an identical or highly similar purpose, we will ask you to contribute to it, even if your pipeline is more mature at the time of submission.
- You might have to make major (possibly breaking) changes to your pipeline to meet [the guidelines](/docs/guidelines/pipelines/overview).
- Due to the nf-core naming policy, the pipeline will most likely need a new, descriptive name within nf-core.

Remember to [start the discussion](https://github.com/nf-core/proposals/issues) as early as possible, and do not despair if your pipeline is not suitable for the community. You can still [contribute to existing pipelines](/docs/tutorials/contributing_to_nf-core/contributing_to_pipelines) or [use nf-core components for your project](/docs/guidelines/external_use).

## Where to start

All nf-core discussion happens on the nf-core Slack, which you can join here:
[https://nf-co.re/join](https://nf-co.re/join)

These topics are specifically discussed in the [nf-core/proposals](https://github.com/nf-core/proposals/issues) repository.
