---
title: 'Tutorial: Using nf-core components outside nf-core'
subtitle: Guidance on how to use nf-core code and best practices in non-nf-core pipelines.
---

While we highly recommend and promote contributing to existing nf-core pipelines, in some cases this may not be possible or appropriate. In this case we hope that it will be helpful to use _parts_ of nf-core best practices for your own non-nf-core pipeline. This tutorial describes how to use nf-core code without making official nf-core pipelines.

## Acknowledging nf-core

If you use any form of nf-core code or infrastructure, we kindly ask that you acknowledge the nf-core initiative by adding the following markdown to your `README.md` and pipeline documentation.

```markdown
This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
> In addition, references of tools and data used in this pipeline are as follows:
```

## Removing nf-core only branding

> Most of the following suggestions refer to when you have started your pipeline using the nf-core template.

## General

- Replace all references to the pipeline name that include the `nf-core/` suffix.
  - :warning: If you still wish to use the `nf-core tools` linting functionality, this may result in a lot of linting failures. You can create a `.nf-core.yml` file that allows you to ignore or skip certain lint tests. You can find more information on how to do this [here](/tools/#linting-config). One example is the `pipeline_name_conventions`

### `README`

To ensure users do not get confused, we suggest removing the following components of the main pipeline `README`:

1. Remove (or replace) the nf-core logo (see [below](#docs))
2. Slack badge for help messages.
3. Replace all references of the [https://nf-co.re](https://nf-co.re) website to the location of your own documentation
4. Remove reference to the nf-core Slack workspace under `Contributions and Support
5. Replace the `Citations` section with the following markdown:

   ```markdown
   This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) initative, and reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

   > The nf-core framework for community-curated bioinformatics pipelines.
   >
   > Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
   >
   > Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.

   In addition, references of tools and data used in this pipeline are as follows:
   ```

Optional:

1. Reference to nf-core/configs in `Quick Start`
   - Note in many cases it is still reasonable to use these configs for non-nf-core pipelines as they typically only modify generic Nextflow options

### `docs/`

1. Replace the nf-core logo under `docs/images` with your own one.
   - You will also have to change the filename/path at the top the main repository README
2. `docs/usage.md`: remove reference to documentation being viewable on the [nf-co.re](https://nf-co.re) website
3. `docs/usage.md`: remove reference to the nf-core Slack workspace

Optional:

1. Reference to nf-core/configs in `-profile`
   - Note in many cases it is still reasonable to use these configs for non-nf-core pipelines as they typically only modify generic Nextflow options

### `.github/`

1. `PULL_REQUEST_TEMPLATE.md`: Remove reference to the nf-core/tests-datasets repository
2. `CONTRIBUTING.md`: remove references to the nf-core Slack workspace
3. `workflows/`: remove nf-core specific GitHub actions YAMLs
   - `push_dockerhub*`
4. `ci.yml`: rename to the name of your pipeline rather than `nf-core CI tests`
5. `ISSUE_TEMPLATE/` Remove references to nf-core website in all templates and config files

Optional:

1. `CONTRIBUTING.md`: remove references to the nf-core/tools linting tests
2. `workflows/`: remove the nf-core linting tests under (Note: you can alternatively keep linting and [configure](#general) to ignore certain tests)
   - `linting.yml`
   - `linting_comment.yml`
3. `workflows/`: remove the AWS test GitHub action workflows

<!-- TODO: Using external module repositories -->
