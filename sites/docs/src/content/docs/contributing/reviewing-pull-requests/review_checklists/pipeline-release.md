---
title: Reviewing pipeline releases
subtitle: Review nf-core pipeline release pull requests
shortTitle: Reviewing pipeline releases
markdownPlugin: checklist
---

Pipeline release reviews ensure that pipelines meet nf-core standards before they reach the community.
When you review a pipeline release pull request, you examine the submission and provide constructive feedback before maintainers merge the changes into the nf-core repository.

Pipeline release PRs vary widely in their contents, making it harder to provide strict universal criteria compared to reviewing [components](/docs/contributors/reviewing-pull-requests/review_checklists/component).
Your role as a reviewer is to check for adherence to the central principles of nf-core: reproducibility, excellent reporting, documentation, and template compliance.

:::tip
Use the [pipeline release review checklist](/docs/contributors/reviewing-pull-requests/review_checklists/pipeline-release) for a quick reference whilst reviewing.
:::

## Before you start

Pipeline release reviews can feel overwhelming, particularly when doing one for the first time.
Remember that you don't need to be an expert in the scientific domain to provide valuable feedback on code quality, documentation, and adherence to nf-core standards.

Focus on the general [nf-core pipeline guidelines](/docs/guidelines/pipelines/overview) rather than the scientific accuracy of the analysis.
Your perspective helps catch issues that developers might have missed.

### Review nf-core principles

Start by verifying that the pipeline meets the general nf-core guidelines.
Check that the pipeline:

- [ ] Follows the nf-core template structure
- [ ] Uses appropriate nf-core conventions and patterns
- [ ] Maintains reproducibility across different computing environments
- [ ] Provides comprehensive reporting and documentation

You can find the complete set of pipeline guidelines in the [pipeline overview](/docs/guidelines/pipelines/overview).

## Review local code and modules

Local modules and scripts require careful attention because they don't benefit from the community review process that nf-core/modules receives.

### Local scripts

Check that scripts in the `bin/` directory include:

- [ ] Author information
- [ ] Licence declarations (prefer MIT over GPL to match pipeline code)

### Local modules

Verify that all local modules in `modules/local/` include:

- [ ] Docker, Singularity, and Conda declarations
- [ ] Tool versions from bioconda/biocontainers where possible
- [ ] Specific version tags (not `latest` or `dev`)
- [ ] Version reporting for complex modules (simple operations like single `grep` commands don't need version reporting)
- [ ] Consider whether local modules belong in `nf-core/modules` instead

  :::note
  If a module could be useful across multiple pipelines, suggest moving it to the shared repository for broader community benefit.
  :::

## Review documentation

Documentation makes pipelines accessible to the community.
Focus your review on these key areas:

### Documentation location

Verify that all documentation:

- [ ] Appears only on the nf-core website (not external platforms or non-standard locations)

### Core documentation files

Check that `usage.md`, `output.md`, and `nextflow_schema.json` provide comprehensive information:

- [ ] **usage.md**: Clear instructions for running the pipeline with common use cases
- [ ] **output.md**: Detailed descriptions of all output files and directories
- [ ] **nextflow_schema.json**: Correct types for all parameters, appropriate use of `default` and `enum` fields

### Changelog and versioning

Ensure the CHANGELOG follows three-point semantic versioning (for example, 2.1.0 rather than 2.1).
Check that:

- [ ] The version number follows SemVer conventions
- [ ] The date has been updated to the release date
- [ ] All significant changes are documented

### Citations

Review `CITATIONS.md` for:

- [ ] Consistent citation formatting
- [ ] Complete citations for all tools used in the pipeline
- [ ] Proper credit to pipeline authors and contributors
- [ ] Pipeline authors list themselves in the manifest and other contributors appear in the README

### Proofreading

Proofread documentation for:

- [ ] Typos
- [ ] Grammatical errors
- [ ] Unclear explanations

  :::note
  Well-written documentation helps users get started quickly and reduces support requests.
  :::

## Review code quality

### General code review

Assess the code for:

- [ ] General readability and clarity
- [ ] Possible bugs or logic errors
- [ ] Non-template components that might cause maintenance issues

Avoid requiring major code optimisation unless it affects functionality.
Small optimisation suggestions are welcome, but large refactoring should have happened during development PRs rather than at release time.

### Parameter conventions

Check that parameters follow nf-core conventions:

- [ ] Use `snake_case` for all parameter names
- [ ] Boolean parameters default to `false` (prefer `skip_step = false` over `run_step = true`)
- [ ] Parameter names are consistent and self-explanatory

### Manifest

Verify that the manifest:

- [ ] Includes all required information, including DOIs where applicable

## What you don't need to review

To keep reviews focused and efficient, you don't need to:

- Review code from nf-core/modules (these modules have their own review process)
- Comment on scientific accuracy unless you're familiar with the topic
- Require major code optimisation at the release stage

  :::note
  You can suggest small code optimisation.
However, larger code optimisation shouldn't necessarily be required for release.
  :::

Remember that your review helps maintain quality across the nf-core ecosystem.
Thank you for contributing your time to make nf-core pipelines better for everyone.
