---
title: Review checklist for pipeline releases
subtitle: A checklist for reviewing pipeline release PRs
shortTitle: Pipelines
markdownPlugin: checklist
---

<!-- TODO: Add links to other pages and guide where possible -->

When you review a PR, you examine a pipeline release submission.
You provide constructive feedback on those changes before the maintainers merge them into the nf-core repository.
Your review ensures that the pipeline meets the coding standards of the project, maintains consistency and achieves high quality.

Use this guide to evaluate nf-core pipeline release pull requests.
Focus on checking for adherence to nf-core principles rather than applying strict universal criteria, given the diverse nature of pipeline releases.

## Do: nf-core principles

- [ ] Verify that the pipeline meets the general nf-core guidelines.

## Do: Local code and modules

- [ ] Check that local scripts in `bin/` include author and licence information.
- [ ] Check that licences are MIT (not GPL) to match pipeline code.
- [ ] Ensure all local modules have docker/singularity/conda declarations, preferably from bioconda/biocontainers.
- [ ] Verify the code specifies tool versions (avoiding `latest` or `dev` tags).
- [ ] Check that complex modules include version reporting.
- [ ] Evaluate whether local modules belong in nf-core/modules instead.

## Do: Documentation

- [ ] Verify documentation appears only on the nf-core website.
- [ ] Check that `usage.md`, `output.md`, and `nextflow_schema.json` are comprehensive.
- [ ] Check schema for correct types, appropriate use of `default` and `enum` fields.
- [ ] Proofread documentation for errors.
- [ ] Ensure CHANGELOG follows three-point semantic versioning (for example, 2.1.0) with updated dates.
- [ ] Verify consistent citation formatting in `CITATIONS.md`.
- [ ] Confirm the documentation properly credits all tools and authors.

## Do: Code

- [ ] Assess general readability and potential bugs.
- [ ] Ensure parameter consistency using `snake_case`.
- [ ] Verify boolean parameters default to `false` (prefer `skip_step = false` over `run_step = true`).
- [ ] Check manifest includes DOI when applicable.

## Don't have to do

- [ ] Review code from nf-core/modules.
- [ ] Comment on scientific accuracy unless qualified.
- [ ] Demand major code optimisation (minor suggestions acceptable).
