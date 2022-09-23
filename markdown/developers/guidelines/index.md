---
title: Guidelines Overview
subtitle: Guidelines and requirements for nf-core pipelines.
menu:
  main:
    weight: 10
---

## Introduction

The aim of nf-core is to have standardised best-practice pipelines.
To ensure this standardisation, we maintain a set of guidelines which all nf-core
pipelines must adhere to.

> If you're thinking of adding a new pipeline to nf-core, please read the documentation
> about [adding a new pipeline](adding_pipelines.md).

The following lists an overview of the guidelines. Follow links to dedicated pages for more detail about a given topic.

### Ask the community

The instructions above are subject to interpretation and specific scenarios.
If in doubt, please ask the community for feedback on the [`#new-pipelines` Slack channel](https://nfcore.slack.com/channels/new-pipelines).
You can join the nf-core Slack [here](/join).

## Requirements

All nf-core pipelines _must_ follow the following guidelines:

- Nextflow: Workflows must be built using Nextflow.
- [Identity and branding](requirements/identity_branding.md): Primary development must on the nf-core organisation.
- [Workflow specificity](requirements/workflow_specificity.md): There should only be a single pipeline per data / analysis type.
- [Workflow size](requirements/workflow_size.md): Not too big, not too small.
- [Workflow name](requirements/workflow_name.md): Names should be lower case and without punctuation.
- [Use the template](requirements/use_the_template.md): All nf-core pipelines must be built using the nf-core template.
- [Software licence](requirements/mit_licence.md): Pipelines must open source, released with the MIT licence.
- [Bundled documentation](requirements/docs.md): Pipeline documentation must be hosted on the nf-core website.
- [Docker support](requirements/docker.md): Software must be bundled using Docker and versioned.
- [Continuous integration testing](requirements/ci_testing.md): Pipelines must run CI tests.
- [Semantic versioning](requirements/semantic_versioning.md): Pipelines must use stable release tags.
- [Standardised parameters](requirements/parameters.md): Strive to have standardised usage.
- [Single command](requirements/single_command.md): Pipelines should run in a single command.
- [Keywords](requirements/keywords.md): Excellent documentation and GitHub repository keywords.
- [Pass lint tests](requirements/linting.md): The pipeline must not have any failures in the `nf-core lint` tests.
- [Credits and Acknowledgements](requirements/acknowledgements.md): Pipelines must properly acknowledge prior work.
- [Minimum inputs](requirements/minimum_inputs.md): Pipelines should be able to run with as little input as possible.
- [Use nf-core git branches](requirements/git_branches.md): Use `master`, `dev` and `TEMPLATE`.

## Recommendations

All nf-core pipelines _should_ follow the following guidelines, if possible / appropriate:

- [Use Bioconda](recommendations/bioconda.md): Package software using bioconda and biocontainers.
- [File formats](recommendations/file_formats.md): Use community accepted modern file formats such as `CRAM`.
- [DOIs](recommendations/dois.md): Pipelines should have digital object identifiers (DOIs).
- [Cloud compatible](recommendations/cloud_compatible.md): Pipelines should be tested on cloud computing environments.

## If the guidelines don't fit

We appreciate that the above guidelines are relatively rigid and may not be for everyone.
If that's the case, there is still a lot of ways that you can get involved with nf-core!

### Non nf-core pipelines

If a pipeline does not conform to nf-core guidelines, it can of course exist outside of nf-core.
You are verly wecome to the helper tools and collaborate on modules / subworkflows / ideas.
If doing this, please don't call your pipeline `nf-core/<yourpipeline>` and please say that your pipeline _"uses"_ nf-core rather than rather than _"is"_ nf-core.
Remember that you can generate a pipeline with `nf-core create` that excludes nf-core branding.

Non nf-core pipelines can be added to [nextflow-io/awesome-nextflow](https://github.com/nextflow-io/awesome-nextflow) for added visibility.

### Enforcement

If a pipeline is found to be violating the nf-core guidelines after it has been added to the community, we will try to address the problems via the following steps:

First, the core team will attempt to resolve problems with the pipeline maintainers through discussion.
Hopefully the pipeline can then be updated so that it adheres to the guidelines.

If this is not possible, the core team will make a recommendation to the steering committee about what action to take.
Such actions could include archiving the pipeline or removing it completely.

All members of the nf-core community must adhere to the [nf-core code of conduct](https://nf-co.re/code_of_conduct).
The guidelines and actions within the code of conduct take precedence over the development pipelines described in this page.
