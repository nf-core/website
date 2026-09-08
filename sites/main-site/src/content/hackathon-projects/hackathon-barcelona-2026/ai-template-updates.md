---
title: AI-assisted template updates
category: tooling
slack: https://nfcore.slack.com/channels/tools
location: Barcelona
image: "/assets/images/events/2026/hackathon-barcelona/ai-template-updates.jpg"
image_alt: "Black and white photo of a laptop showing a software update screen. Photo by Clint Patterson (@cbpsc1) on Unsplash."
leaders:
  ewels:
    name: Phil Ewels
    slack: https://nfcore.slack.com/team/UE4L64T4G
---

## Goal

Get as many nf-core pipelines as possible onto the latest template, using agentic coding tools to do the merge work.
We'll build agent skills and instructions as we go based on the experience, to make syncs easier for people in the future.

## Description

Keeping pipelines up to date with the nf-core/tools pipeline template is a central aspect of how the nf-core community works.
However, most people agree that it's a bit of a pain.
We have automated tooling to try to synchronise updates to each pipeline which creates a pull-request, but these invariably contain numerous merge-conflicts that need to be picked apart and understood.

As a result, many pipelines lag behind the latest available template (see [Pipeline Versions](https://nf-co.re/pipeline_versions)), some significantly so.
These pipelines don't get security updates, feature updates and all the other goodies that ship in the central template.
Worse still, since Nextflow is getting new syntax updates, many are becoming impossible to run with recent versions of Nextflow.

It feels like AI agents should be good at the task of template updates, but most people who have tried report that it doesn't work as well as you'd expect.
In this project we will try to write up some new agent skills for updates, then run them and see how they work.
Every time we have problems, we will fix and wrap what we learnt back into the skill and try again.

By the end of the hackathon, we hope to have updated a good number of pipelines to the latest template,
and also produced a reusable skill that others can make use of to sync pipelines in the future.
Once synced, pipelines can also be updated to work with strict syntax.

Easier sync = more up to date pipelines = better nf-core pipelines.

## Tasks

We'll start by discussing what usually goes into pipeline syncs, and co-drafting an initial skill for the process.
Once that's done, we'll move onto tackling updates in parallel:

:::info{title="Great first issues"}
Reviewing someone else's template update PR is one of the best ways to learn what the template actually contains, and it needs no agent setup.
:::

- Claim a pipeline in the project notes so nobody duplicates work
- Check for existing automated PRs, and / or run `nf-core pipelines sync` and see what comes out
- Hand the conflicts to an agent along with the pipeline's own context (`AGENTS.md`, recent commits, why the file was customised), then get CI green
- Review incoming template update PRs from other people in the group
- Write down where the agent went wrong, and what context would have prevented it

## Making PRs

Template update PRs go against each pipeline's `dev` branch.
Unless you're very confident in the code changes, make sure that you get a maintainer of that pipeline to review before merging.

Syntax and tooling fixes that apply to every pipeline should go upstream to the pipeline template in [nf-core/tools](https://github.com/nf-core/tools) instead.
