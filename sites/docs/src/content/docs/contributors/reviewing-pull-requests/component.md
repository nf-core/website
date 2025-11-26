---
title: Review checklist for nf-core components
subtitle: A checklist for reviewing PRs in the nf-core/modules repository
shortTitle: Component PR
markdownPlugin: checklist
---

<!-- TODO: Add links to other pages and guide where possible -->

When you review a PR, you examine a new component submission or proposed changes to a component. You provide constructive feedback on those changes before the maintainers merge them into the nf-core repository. Your review ensures that the code meets the coding standards of the project, maintains consistency and achieves high quality.

While the team of [maintainers](https://github.com/orgs/nf-core/teams/maintainers/members) oversees the PR review process for components, these guidelines help you review PRs consistently and effectively as a community member. The following collection of community suggestions can guide your review process.

## General

- [ ] All components adhere to the nf-core [module specifications](/docs/guidelines/components/modules) or [subworkflow specifications](/docs/guidelines/components/subworkflows).
- [ ] Ensure all checks pass, including linting, conda, singularity, and docker.

Otherwise, you can cover most of the specifications by checking for the following:

- [ ] Check that the component runs offline without assuming automatic database downloads.
- [ ] Check that Nextflow changes the `--entrypoint` to `/bin/bash` and that environment variables are sourced again when running docker containers.
- [ ] Check that it adheres to nf-core coding standards (for example, use of meta map).
- [ ] Check that the code is readable and the formatting is correct (for example, indenting, extra spaces).

## `modules/nf-core/modulename/main.nf`

- [ ] Check that all optional parameters are in the `$args` section.
- [ ] Check that the software version extraction command is optimised, if required.
- [ ] Check if the bioconda version of the tool is the latest version.
- [ ] Ensure the code removes temporary unzipped files to avoid mitigating benefits and worsening problems.
- [ ] Ensure large outputs use the correct compression tool (follow guidelines for gzip vs bzip2 vs other options).

## `../tests/modules/nf-core/modulename/main.nf` and `meta.yml`

- [ ] Check that there are tests for all outputs, including optional ones.
- [ ] Check that the `meta.yml` file has correct documentation links and patterns of files.
- [ ] Check that the `meta.yml` file has correct [`bio.tools`](https://bio.tools/) ID for the tool and correct [EDAM ontology](https://edamontology.github.io/edam-browser/#topic_0091) links of files.
- [ ] Run the tool help and check that important input (usually optional) has not been missed.
- [ ] Run nf-test (for example, on Gitpod) to check that it captures all outputs.
