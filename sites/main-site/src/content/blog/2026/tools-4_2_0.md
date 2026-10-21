---
title: "nf-core/tools - 4.2.0"
subtitle: "Next stop: pipeline metro maps"
headerImage: "https://images.unsplash.com/photo-1599143527645-34980fcb6315"
headerImageAlt: "A close-up picture of the metro map for São Paulo"
pubDate: 2026-10-21T12:00:00+02:00
authors:
  - "mashehu"
label:
  - "tools"
maxHeadingDepth: 3
embedHeaderImage: true
---

This is a small release of nf-core/tools: mostly template maintenance and bug fixes, plus pipeline metro map support and a new helper command for container configs.

As always, if you have any problems or run into any bugs, reach out on the [#tools slack channel](https://nfcore.slack.com/channels/tools).

# Highlights

- [Metro maps with `manifest.diagram{:groovy}`](#metro-maps-with-manifestdiagram): define your pipeline's metro map in `nextflow.config`.
- [`nf-core pipelines containers create-configs` command](#nf-core-pipelines-containers-create-configs-command): regenerate the `conf/containers_*.config` files on demand.

## Metro maps with `manifest.diagram{:groovy}`

Nextflow 26.10.0 added a `manifest.diagram{:groovy}` field, which points to a workflow diagram of the pipeline. This will help different tools to easier pick up the tube maps, e.g. the RO-crates or [Seqera Platform](https://docs.seqera.io/platform-enterprise/launch/launchpad#pipeline-diagram)
The template now has a commented-out `manifest.diagram{:groovy}` in `nextflow.config`, ready for your [metro map](/docs/community/brand/workflow-schematics):

```groovy title="nextflow.config"
manifest {
    ...
    // diagram = 'docs/images/metro_map.svg'
}
```

The `nextflow_config` lint test warns if `manifest.diagram{:groovy}` is not set, but no pressure if you don't want to add one (even if [all the cool kids do it](/docs/community/brand/workflow-schematics#examples)).
If it is set, it must be a relative path to an existing image file inside the pipeline (`.svg`, `.png`, `.jpg`, `.jpeg`, `.gif` or `.webp`), or the test fails.
Older Nextflow versions ignore the field, so you can set it now.

Don't have a metro map yet? See all the ways to create one in the [workflow schematics docs](/docs/community/brand/workflow-schematics).

## `nf-core pipelines containers create-configs{:bash}` command

`nf-core pipelines lint{:bash}` regenerates the `conf/containers_*.config` files, but it refuses to run on a git tree with uncommitted changes.
That is a good safety net in general, but can be annoying.
The new command directly regenerates these files on demand:

```bash
nf-core pipelines containers create-configs
```

## Miscellaneous

- `ext.use_gpu{:groovy}` has been superseded by Nextflow's `accelerator` directive. We removed it from the `process_gpu` label in the template, and it is no longer a permitted `ext` key in module linting. Use `task.accelerator{:groovy}` instead, as described in [the docs](/docs/developing/components/gpu-modules)
- nf-schema is bumped to 2.7.3, which is needed to run pipelines with Nextflow 26.04. The nf-schema version used in the `awsfulltest` workflow is now pinned as well, [to avoid an annoying `Missing process or function paramsSummaryLog(...){:console}` error](https://github.com/nf-core/tools/issues/4487).
- The `fix_linting.yml` workflow is now `fix-linting.yml` and uses the reusable workflow from the [nf-core/actions](https://github.com/nf-core/actions) repository, so the nf-core bot token can't leak. This was already patched into the `main`/`master` branch of all nf-core pipelines, so now the `dev` branch will also have it.
- `workflow.onComplete` no longer waits for the MultiQC report. If no report was produced, the pipeline doesn't hang at the end anymore.
- nf-test files are now git-ignored in the template.
- The template's VS Code settings recommend now extensions useful for developing for nf-core.
- `strandedness{:groovy}` is now a permitted meta key in module linting (again).

# Changelog

You can find the complete changelog and technical details [on GitHub](https://github.com/nf-core/tools/releases/tag/4.2.0).
<!--
# How to resolve common merge conflicts on pipeline sync PRs

All the template changes are listed in [the Template section of the changelog](https://github.com/nf-core/tools/releases/tag/4.2.0).
Below, we have collected the most common merge conflicts that you may find and how to resolve them.

{/* TODO: add common merge conflicts once the sync PRs are out */}-->
