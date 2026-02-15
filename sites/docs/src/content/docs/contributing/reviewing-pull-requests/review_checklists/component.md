---
title: Reviewing nf-core components
subtitle: Review component pull requests
shortTitle: Reviewing components
markdownPlugin: checklist
---

<!-- TODO: Add links to other pages and guide where possible -->

Component reviews ensure that modules and subworkflows meet nf-core standards before they reach the community.
When you review a component pull request, you examine a new component submission or proposed changes to an existing component and provide constructive feedback before maintainers merge them into the nf-core repository.

The team of [maintainers](https://github.com/orgs/nf-core/teams/maintainers/members) oversees the review process for components, but community input helps catch issues and ensures components work well across different use cases.
Your perspective as a user of the components is valuable, particularly if you've used the tool before or understand common use cases.

## General

Start by verifying that the component meets the general specifications:

- [ ] All components adhere to the nf-core [module specifications](/docs/guidelines/components/modules) or [subworkflow specifications](/docs/guidelines/components/subworkflows)
- [ ] All checks pass, including linting, conda, singularity, and docker

You can cover most of the specifications by checking for the following:

- [ ] Component runs offline without assuming automatic database downloads
- [ ] Nextflow changes the `--entrypoint` to `/bin/bash` and environment variables are sourced again when running docker containers
- [ ] Code adheres to nf-core coding standards (for example, use of meta map)
- [ ] Code is readable and the formatting is correct (for example, indenting, extra spaces)

## Review the main module file

Check `modules/nf-core/<module_name>/main.nf` for:

- [ ] All optional parameters are in the `$args` section
- [ ] Software version extraction command is optimised, if required
- [ ] Bioconda version of the tool is the latest version
- [ ] Code removes temporary unzipped files to avoid mitigating benefits and worsening problems
- [ ] Large outputs use the correct compression tool (follow guidelines for gzip vs bzip2 vs other options)

  :::tip
  Check the tool's documentation or help text to ensure important optional parameters haven't been missed.
  Run the tool help command if you're unsure about available options.
  :::

## Review tests and metadata

Check `../tests/modules/nf-core/<module_name>/main.nf` and `meta.yml` for:

- [ ] Tests exist for all outputs, including optional ones
- [ ] `meta.yml` file has correct documentation links and patterns of files
- [ ] `meta.yml` file has correct [`bio.tools`](https://bio.tools/) ID for the tool and correct [EDAM ontology](https://edamontology.github.io/edam-browser/#topic_0091) links of files
- [ ] Tool help has been checked to ensure important input (usually optional) has not been missed
- [ ] nf-test runs successfully (for example, on Gitpod) and captures all outputs
