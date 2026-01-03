---
title: nf-core specifications
subtitle: Guidelines and requirements for nf-core contributions and best practices.
weight: 1
---

The nf-core specifications define standards and best practices for developing robust, reproducible, and maintainable Nextflow components and pipelines. While these specifications are written primarily for nf-core contributors and are enforced in nf-core repositories, the principles and patterns documented here represent proven approaches to writing high-quality bioinformatics workflows that can benefit the broader Nextflow community.

These specifications are the result of collaborative development by the nf-core community. They address common challenges in workflow development, from managing software dependencies and computational resources to ensuring reproducibility and facilitating community contributions.

## Specification categories

The nf-core specifications are organized into three main categories:

### Components

Standards for developing reusable nf-core modules and subworkflows. These specifications ensure consistency across the nf-core component library and facilitate sharing and maintenance of workflow building blocks.

- **[Module specifications](/developers/specifications/components/modules):** Guidelines for developing individual process wrappers
- **[Subworkflow specifications](/developers/specifications/components/subworkflows):** Standards for combining modules into reusable workflow units

### Pipelines

Requirements and recommendations for nf-core end-to-end analysis pipelines. These specifications define what it means to be an nf-core pipeline and provide best practice guidance.

- **[Pipeline requirements](/developers/specifications/pipelines/requirements):** Mandatory standards all nf-core pipelines MUST follow
- **[Pipeline recommendations](/developers/specifications/pipelines/recommendations):** Best practices nf-core pipelines SHOULD adopt where applicable

### Test Data

Guidelines for managing test data used in continuous integration and testing of nf-core components and pipelines.

- **[Test data specifications](/developers/specifications/test-data):** Standards for test data organization, size, licensing, and documentation

## Ask the community

The instructions below are subject to interpretation and specific scenarios.
If in doubt, ask the community for feedback on the [#help slack channel](https://nfcore.slack.com/channels/help).

## If the guidelines don't fit

The guidelines are relatively rigid and may not be for everyone.
If that's the case, there is still a lot of ways that you can get involved with nf-core!

The nf-core best practices, tooling and community are helpful for anyone building Nextflow pipelines, even if they are not a good fit for being listed as official nf-core pipelines.
You are very welcome to use the helper tools and collaborate on modules / subworkflows / ideas.
Indeed, numerous pipelines outside of nf-core now extensively use and contribute back to [nf-core/modules](https://github.com/nf-core/modules).

If using nf-core tools and especially the template, you MUST NOT call your pipeline `nf-core/<pipeline>`.
Your pipeline SHOULD be described as "using" nf-core rather than "being" nf-core.
Remember that you can generate a pipeline with `nf-core pipelines create` that excludes nf-core branding.
Citation and acknowledgement of the work that goes into these tools and templates is welcome.

If a pipeline is found to be violating the nf-core guidelines _after_ it has been added to the community, issues will be addressed:

1. The core team will attempt to resolve problems with the pipeline maintainers through discussion. Hopefully the pipeline can then be updated so that it adheres to the guidelines.
1. If this is not possible, the core team will make a recommendation to the steering committee about what action to take. Such actions could include archiving the pipeline or removing it completely.

All members of the nf-core community must adhere to the [nf-core code of conduct](https://nf-co.re/code_of_conduct).
The guidelines and actions within the code of conduct take precedence over the development pipelines described in this page.
