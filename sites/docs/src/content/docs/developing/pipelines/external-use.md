---
title: External use of nf-core resources
subtitle: Guidelines for using nf-core code and infrastructure outside the organisation
shortTitle: External use
weight: 10
---

<!-- TODO: Add links to this page from contributing and developing pages, maybe move this to developing? -->

Whenever possible, you should contribute to existing nf-core pipelines. However, in some cases, this may not be possible or appropriate. In this case, you should still use parts of nf-core best practices for your own non-nf-core pipeline. This guide describes how to use nf-core code without making official nf-core pipelines.

If you use nf-core code or infrastructure for your own pipeline development, you must acknowledge the nf-core initiative and follow these guidelines for proper attribution and branding.

## Acknowledge nf-core

When using nf-core code or infrastructure, include this acknowledgment in your pipeline documentation:

```
This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
> In addition, references of tools and data used in this pipeline are as follows:
```

## Remove nf-core branding

Most of the following suggestions refer to when you have started your pipeline using the nf-core template.

### General

Replace the following:

- References to the pipeline name that include the `nf-core/` suffix.

  :::warning
  If you still wish to use the `nf-core tools` linting functionality, this may cause linting failures. Create a `.nf-core.yml` file that allows you to ignore or skip certain lint tests. See linting config for more information. <!-- TODO: Add link to linting config -->
  :::

### `README`

Remove the following:

- The nf-core logo
- Slack badge for help messages
- References to the [https://nf-co.re](https://nf-co.re) website
- References to the nf-core Slack workspace
- (Optional) Reference to nf-core/configs in `Quick Start`

  :::note
  In many cases it is still reasonable to use these configs for non-nf-core pipelines as they typically only modify generic Nextflow options.
  :::

Replace the following:

- The `Citations` section with the following:

```markdown
This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) initative, and reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.

In addition, references of tools and data used in this pipeline are as follows:
```

### `docs/`

Remove the following:

- `docs/usage.md`: Reference to documentation being viewable on the [nf-co.re](https://nf-co.re) website
- `docs/usage.md`: Reference to the nf-core Slack workspace
- (Optional) Reference to nf-core/configs in `-profile`

  :::note
  In many cases it is still reasonable to use these configs for non-nf-core pipelines as they typically only modify generic Nextflow options.
  :::

Replace the following:

- The nf-core logo under `docs/images`

  :::note
  You will also have to change the filename/path at the top the main repository `README`
  :::

### `.github/`

Remove the following:

- `PULL_REQUEST_TEMPLATE.md`: Reference to the nf-core/tests-datasets repository
- `CONTRIBUTING.md`: References to the nf-core Slack workspace
- `CONTRIBUTING.md`: (Optional) References to the nf-core/tools linting tests
- `workflows/`: nf-core specific GitHub actions YAMLs
  - `push_dockerhub*`
- `workflows/`: (Optional) nf-core linting tests under (Note: you can alternatively keep linting and [configure](#general) to ignore certain tests)
  - `linting.yml`
  - `linting_comment.yml`
- `workflows/`: (Optiona) AWS test GitHub action workflows
- `ISSUE_TEMPLATE/` Remove references to nf-core website in all templates and config files

Rename the following:

- `ci.yml`: The name of your pipeline (rather than `nf-core CI tests`)

<!-- TODO: Using external module repositories -->
