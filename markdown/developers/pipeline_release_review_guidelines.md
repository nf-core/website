---
title: Pipeline Release Review Guidelines
subtitle: Suggestions for reviewing pipeline release PRs
---

Pipeline release reviews are often quite an overwhelming task, particularly when doing one for the first time.

Pipeline release PRs are often quite diverse in their contents, making it harder to provide a strict set of reviewing guidelines that the PR must meet for approval, when compared to [modules](modules) or [subworkflows](subworkflows).

Overall, the role of the reviewer of pipeline release PRs is to check for adherence to the central principles of nf-core (reproducibility, execellent reporting, documented, keeping to the template etc.,). These are detailed in the pipeline guidelines [here](guidelines)

Here we provide a general set of suggestions of what a reviewer should and do not need to look for when doing pipeline releases:

## Do

### nf-core principles

- Does the pipeline meet the general [nf-core guidelines](guidelines)

### Local code and modules

- Do local scripts in `bin/` have author and license embedded?
  - Local script licenses if written for the pipeline should not GPL, but ideally MIT (matching with pipeline code, upon which is assumed if not otherwise specified)
- Do all local modules have docker/singularity/conda declarations?
  - Are they ideally in bioconda/biocontainers
- Do all local modules conda/container tool declarations have versions? (and _not_ `latest`, `dev` etc.)
- Do all local modules report versions (if applicable)?
  - Simple modules with e.g. single `grep` operations not necessary
  - It would be good to add with more complex operations such as `awk`
- Should any local modules be in nf-core/modules?

### Documentation

- Documention is only on the nf-core website (not pointed to other non-nf-core standard places, e.g. not read the docs )
- Is documentation sufficiently described (`usage.md`, `output.md`, `nextflow_schema.json`)?
  - nextflow_schema.json: check if types are correct and that `default` and `enum` are used where applicable
- Are there any typos in the documentation (`usage.md`, `output.md`, `nextflow_schema.json`)
- Is CHANGELOG sufficiently filled in?
  - Check version system is three-point SemVer e.g. 2.1.0
  - Has the date been updated?
- Check citation formatting consistency in `CIATIONS.md`
- Check that all tools are cited
- Check that (all) pipeline author(s) listed themselves in the manifest and other contributors are added in the README

### Code

- Check no overly non-template components (no read the docs, entirely custom logo etc.)
- Check for general code readability
- Check for possible code bugs
- Check for consistency in parameters
  - i.e. `snake_case`
  - All boolean parameters evaluate to `false` (e.g. bad: `params.run_step = true`, good: `params.skip_step = false` )
- Check manifest includes DOI (if present) etc.

## Don't have to do

- Comment on scientific content (unless you are familiar with the topic)
- Major code optimisation
  - You _can_ suggest small code optimisations
  - Larger ones you can recommend, but shouldn't necessarily be _require_ for release (hopefully this should have been covered during development PRs)
