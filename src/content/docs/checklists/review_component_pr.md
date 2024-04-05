---
title: PR Review Checklist for nf-core components
subtitle: Guidelines for reviewing PRs in the nf-core/modules repository
---

A PR review is the process of examining a new component submission or the changes proposed to a component. The reviewer provides constructive feedback on those changes before they are merged into the nf-core repository.The goal of a PR review is to ensure that the code meets the coding standards of the project, is consistent and of high-quality.

While the team of [maintainers](https://github.com/orgs/nf-core/teams/maintainers/members) is responsible for overseeing the PR review process for components, these guidelines can assist community members in reviewing PRs and ensure that the review process is consistent and effective. The following is a collection of community suggestions to have into account during the review process.

## General reviews of submissions to components:

In general, the main purpose of the review is to ensure

- All components adhere to the nf-core [module specifications](/docs/contributing/component-specifications/modules) or [subworkflow specifications](/docs/contributing/component-specifications/subworkflows).
- Ensure all checks pass, including linting, conda, singularity, and docker.

Otherwise, you can cover most of the specifications by checking for the following:

- The component is suitable for offline running, without automatic database downloads assumed.
- If running docker containers, check that Nextflow changes the `--entrypoint` to `/bin/bash` and that environment variables used by certain programs (e.g., Busco, Merqury) are sourced again to use them in container settings.
- Check that it adheres to nf-core coding standards (e.g. use of meta map).
- Check that the code is readable and the formatting is correct (e.g. indenting, extra spaces).

## In `modules/nf-core/modulename/main.nf`:

- Check that all optional parameters are in the `$args` section.
- Check that the software version extraction command is optimized, if required.
- Check if the bioconda version of the tool is the latest version.
- Ensure that temporary unzipped files are removed to avoid mitigating benefits and worsening problems.
- Ensure that large outputs are compressed with the correct tool (follow guidelines for gzip vs bzip2 vs other options).

## In `../tests/modules/nf-core/modulename/main.nf` and `../tests/modules/nf-core/modulename/meta.yml`:

- Check that there are tests for all outputs, including optional ones.
- Check that the `meta.yml` file has correct documentation links and patterns of files.
- Run the tool help and check that important input (usually optional) has not been missed.
- Check that all outputs are captured by running nf-test (e.g. on Gitpod).
