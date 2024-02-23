---
title: 'nf-core/tools - 2.13.0'
subtitle: 'Out with the old, in with the new'
pubDate: 2024-02-20T00:00:00+01:00
headerImage: 'https://images.unsplash.com/photo-1535261170473-22e6c4a35c43'
headerImageAlt: 'Wodden attic with empty shelves'
authors:
  - 'mashehu'
  - 'mirpedrol'
  - 'drpatelh'
  - 'maxulysse'
label:
  - 'tools'
embedHeaderImage: true
---

This release contains some template changes and bug fixes.

# Highlights

We refactored the pipelines template a bit:

- The `lib` directory is removed: :wave: Groovy code
- Instead we now use nf-core subworkflows for pipeline initialisation
  - Some nf-core pipelines adopted this already: [nf-core/fetchngs](https://github.com/nf-core/fetchngs/tree/dev) and [nf-core/rnaseq](https://github.com/nf-core/rnaseq/tree/dev)
- The [`nf-validation` plugin](https://nextflow-io.github.io/nf-validation/1.1/) is now used to create an input channel from a sample sheet file.

:::tip
If you haven't started merging in the template updates yet, it may almost be easier to ignore the template updates for now and try and remove the `lib/` directory on the pipeline `dev` branch by:

1. Individually installing the `utils_*` subworkflows with the following commands:

```bash
nf-core subworkflows install utils_nextflow_pipeline
nf-core subworkflows install utils_nfcore_pipeline
nf-core subworkflows install utils_nfvalidation_plugin
```

2. Creating a local `utils_*` subworkflow for the pipeline. You can copy the one in [rnaseq (`dev`)](https://github.com/nf-core/rnaseq/blob/dev/subworkflows/local/utils_nfcore_rnaseq_pipeline/main.nf) or the [pipeline template](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/subworkflows/local/utils_nfcore_pipeline_pipeline/main.nf) and customise this to your requirements. Make sure you move any custom functions in `lib/` directory to this file.
3. Include the `utils_*` subworkflows in the main.nf as done in [rnaseq (`dev`)](https://github.com/nf-core/rnaseq/blob/48663bffadb900e1ae4e11fb3391134cbf12ffc7/main.nf#L22-L25).
4. Including the `utils_*` subworkflows in the workflow main.nf as done in [rnaseq (`dev`)](https://github.com/nf-core/rnaseq/blob/48663bffadb900e1ae4e11fb3391134cbf12ffc7/workflows/rnaseq/main.nf#L25-L30).
5. Delete the `lib/` directory after you have confirmed everything has been ported across.
6. Once you have merged this to `dev` the template sync PR will be updated and tell you whether you have missed anything.
7. The `nf-core lint` command might complain about having to recompute checksum of subworkflow(s).
    - Be sure to check the `modules.json` file as previously installed subworkflows might have dissapeared from it.
    - The `git merge` command adds a new `subworkflows` sections with the new `utils_*` subworkflows, which doesn't create a merge conflict but clashes with any previously installed subworfklows.
    - Reinstalling the subworkflows using `nf-core subworkflows install` should fix this, otherwise one could manually edit the modules.json file (only recommended for advanced users).

:bulb: It helped to disable running the main workflow whilst wiring all of this in to speed up development.
:::

## Additional bug fixes and improvements

- Have you seen a `no space left on device` error or your CI tests lately? This is now fixed by adding a clean up step in the github actions.

- [New api-docs](https://nf-co.re/tools/docs/latest) on the nf-core website to get clarifications on linting errors etc.

- Updating a module won’t remove custom nextflow.config files anymore.

- The commands `nf-core modules test` and `nf-core subworkflows test` now have a new parameter `--profile` to specify the container to run nf-test with.

Thanks to all the contributors who made this release possible, especially to [Harshil Patel](https://github.com/drpatelh/) and his team for kicking off the `lib` restructuring and the `utils_*` subworkflows.

You can find the whole changelog [on GitHub](https://github.com/nf-core/tools/releases/tag/2.13).
