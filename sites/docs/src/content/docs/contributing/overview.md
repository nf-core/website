---
title: Overview
subtitle: Contribute to nf-core
shortTitle: Overview
weight: 1
---

nf-core is a community-driven project built on the contributions.
Whether you're fixing a bug, developing a new pipeline, writing documentation, or reviewing pull requests, every contribution helps advance reproducible analysis for the entire community.

This section covers the different ways you can contribute to nf-core.

## Pipelines

nf-core pipelines are end-to-end bioinformatics workflows built using the nf-core template and community standards.

- **[Contributing a new pipeline](./contribute-new-pipelines):** Propose, develop, and release a new nf-core pipeline, from initial proposal through to your first release
- **[Contributing to existing pipelines](./contribute-existing-pipelines):** Add features, fix bugs, or improve documentation in existing nf-core pipelines

## Components

nf-core modules and subworkflows are reusable components shared across pipelines.
Contributing components benefits the entire Nextflow community by reducing duplicated effort and enabling standardized, reproducible analyses:

- **[Contributing components](./contribute-components):** Contribute modules and subworkflows to the nf-core/modules repository

## Tools

nf-core/tool is a Python package with helper tools for developing and maintaining nf-core pipelines, modules, and subworkflows:

- **[nf-core/tools documentation](https://docs-v2--nf-core-docs.netlify.app/docs/nf-core-tools/cli/installation):** Documentation of the nf-core/tools command-line tooling
- **[Contributing tools](https://nf-co.re/docs/nf-core-tools):** Contribute to the nf-core/tool package

## Configs

nf-core/configs contains a set of centralised Nextflow configuration files that can be used with all nf-core pipelines for specific clusters or compute environments at different institutions:

- **[Contributing configs](https://github.com/nf-core/configs):** Add a new institutional configuration file with custom parameters to run nf-core pipelines

<!-- TODO: Add a section about plugins when docs become available in the future -->

## Documentation

Good documentation is essential to nf-core's mission.
Documentation contributions are just as valuable as code contributions:

- **[Contributing documentation](./documentation):** Fix typos, clarify instructions, or write new guides for nf-core documentation

## Reviewing pull requests

Reviewing contributions is an important part of nf-core's collaborative development process:

- **[Reviewing components](./reviewing-pull-requests/review_checklists/component):** Checklist for reviewing module and subworkflow pull requests
- **[Reviewing pipeline releases](./reviewing-pull-requests/review_checklists/pipeline-release):** Checklist for reviewing pipeline release pull requests
- **[Reviewing nf-core/tools](./reviewing-pull-requests/review_checklists/nf-core-tools):** Checklist for reviewing nf-core/tools pull requests
- **[nf-core-bot](./reviewing-pull-requests/nf-core-bot):** Use the nf-core bot to control automated GitHub actions on pull requests

## Project proposals

Large projects that affect a significant proportion of the community go through the nf-core RFC (Request for Comment) process:

- **[Project proposals](./project-proposals):** Propose and develop major changes through the nf-core RFC process, from initial suggestion through to implementation

## Other considerations

These pages cover additional topics that may be relevant depending on your contribution:

- **[Deprecating modules](../developing/components/deprecating-components):** Mark outdated modules or subworkflows as deprecated when a better alternative is available
- **[Contributor types](./contributor-types):** Understand how nf-core attributes and recognizes different types of contributors in pipeline manifests
- **[Contributor's list](./contributors-list):** Add your institution to the [nf-core website](https://nf-co.re/contributors#organisations) to be featured

## Ask the community

If you're unsure where to start or have questions about contributing, reach out on [nf-core Slack](https://nf-co.re/join).
The [#help](https://nfcore.slack.com/channels/help) channel is a great place to ask general questions, and pipeline-specific channels are available for more targeted discussions.
