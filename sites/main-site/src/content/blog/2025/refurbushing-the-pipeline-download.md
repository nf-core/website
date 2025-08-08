---
title: "Refurbishing the pipeline downloads command"
subtitle: "Why are my pipeline download tests failing?"
headerImage: "https://images.unsplash.com/photo-1636819488537-a9b1ffb315ce"
headerImageAlt: "3D rendering of a cute little space rocket painted red and white with a fluffy smoke cloud coming out of it."
pubDate:
authors:
  - "ErikDanielsson"
label:
  - "tools"
maxHeadingDepth: 3
---

import old_nf_version from "@assets/images/blog/pipeline-download-refactor/old_nf_version.png";
import rnaseq_in_action from "@assets/images/blog/pipeline-download-refactor/rnaseq_in_action.png";

import { Image } from "astro:assets";

The `nf-core pipelines download` command, used when you want to run an `nf-core` pipeline in an offline compute environment, has recently undergone a [substantial refactor](https://github.com/nf-core/tools/pull/3634).
The command has for a long time been reliant on complex regex logic to capture container strings in modules and config files.
Due to the variability in where a pipeline can define containers, this has been prone to breaking.

Not anymore!
As of [Nextflow 25.04](https://www.nextflow.io/docs/latest/migrations/25-04.html) the `nextflow inspect` command has been substantially extended to allow fetching of all pipeline containers by executing a dry-run of the pipeline.
This has now been integrated into the `nf-core/tools` code base which means that the download command can offload the tricky problem of finding pipeline containers to Nextflow, making the command simpler and more reliant at the same time.

:::warning{title=Why are my pipeline download tests failing?!}

<Image src={old_nf_version} alt="Old Nextflow version" height="300" />

`nf-core` pipelines use a GitHub workflow that runs the `nf-core pipelines download` to check that a pipeline is downloadable before release (see for example the [`nf-core/rnaseq` workflow](https://github.com/nf-core/rnaseq/actions/workflows/download_pipeline.yml)).
This workflow uses the **dev** branch of `nf-core/tools`, originally to allow maintainers of `nf-core/tools` to quickly push patches when the regex logic broke.

However, since the refurbished downloads command requires that the pipeline supports the 25.04 release of Nextflow, any release of a pipeline with an older Nextflow version will fail.

:::

### Added support for downloading Docker containers

<Image src={rnaseq_in_action} alt="Docker container download for `rnaseq` 3.19.0" height="300" />

In addition to the use of `nextflow inspect` for container detection, the download command now also supports downloads of docker container images.
This means that `nf-core` pipelines can now be run on offline HPCs which support only Docker or Podman container systems.
This has been [requested within the `nf-core` community](https://github.com/nf-core/tools/issues/2309), and will allow `nf-core` pipelines to be able to run on more systems than before.

:::tip{title="How are Docker images saved?"}

Compared to Singularity containers which are just normal files kept on your file system, Docker images are generally handled by the Docker daemon.
However, via the [`docker image save`](https://docs.docker.com/reference/cli/docker/image/save/) command, Docker allows packaging images as `tar` archives.
The `tar` archives can subsequently be loaded into another Docker daemon with [`docker image load`](https://docs.docker.com/reference/cli/docker/image/load/), or if you are running Podman with [`podman load`](https://docs.podman.io/en/v5.0.2/markdown/podman-load.1.html).

The `nf-core pipelines downloads` command creates `tar` archives for all Docker images within the pipeline.

:::

:::warning{title="`--parallel-download n` becomes `--parallel n`"}

The `--parallel-download` flag is changed to `--parallel` to better reflect the semantics:
The flag now represents that `n` threads can be used to _either_ download or pull container images, depending on the type of container deteced.
The short flag `-d` is unchanged.

:::
