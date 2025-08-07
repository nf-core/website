---
title: "Refurbishing the pipeline downloads command"
subtitle: ""
headerImage: "https://images.unsplash.com/photo-1636819488537-a9b1ffb315ce"
headerImageAlt: "3D rendering of a cute little space rocket painted red and white with a fluffy smoke cloud coming out of it."
pubDate: 
authors:
    - "ErikDanielsson"
label:
    - "tools"
maxHeadingDepth: 3
---

:::tip{.fa-hourglass-clock title="Too Long; Didn't Read"}

- Pipeline downloads tests are failing since the `nf-core pipelines download` command relies on [Nextflow 25.04](https://www.nextflow.io/docs/latest/migrations/25-04.html) due to [recent changes major in the download commnd](https://github.com/nf-core/tools/pull/3634).
- The `nf-core pipelines download` command now support downloading pipeline Docker containers in addition to Singularity containers.

:::

## Introducing `nextflow inspect` in the `nf-core/tools` download command 
The `nf-core pipeline download` command, used when you want to run an `nf-core` pipeline in an offline compute environment has for a long time been reliant on complex regex logic for fetching your pipeline's container. 
Due to the variability in where pipelines define containers, this has been prone to break.
Not anymore!
As of [Nextflow 25.04](https://www.nextflow.io/docs/latest/migrations/25-04.html) the `nextflow inspect` command has been substantially extended to allow fetching of all pipeline containers by executing a dry-run of the pipeline.
This means that the `nf-core pipelines download` code can offload the tricky problem of finding pipeline containers to Nextflow, making the command more reliant at the same time.

:::warning{title=Why are my pipeline download tests failing?!}

`nf-core` pipelines use a GitHub workflow that runs the `nf-core pipelines download` to check that a pipeline is downloadable before release (see for example the [`nf-core/rnaseq` workflow](https://github.com/nf-core/rnaseq/actions/workflows/download_pipeline.yml)).
This workflow uses the **dev** branch of `nf-core/tools`, originally to allow maintainers of `nf-core/tools` to quickly push patches when the regex logic broke.

However, since the refurbished downloads command requires that the pipeline supports the 25.04 release of Nextflow, any release of a pipeline with an older Nextflow version will fail.

:::

## Added support for downloading Docker containers
In addition to the use of `nextflow inspect` for container detection, the download command now also supports downloads of docker container images. 
This means that `nf-core` pipelines can now be run on offline HPCs which only support Docker or Podman container systems, which has been [requested within the `nf-core` community](https://github.com/nf-core/tools/issues/2309)

:::info{title="Saving docker images" collapse}

Compared to Singularity containers which are just normal files, Docker images are generally handled only by the Docker daemon. 

:::

## 
