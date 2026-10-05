---
title: New container syntax and Wave adoption
category: components
location: Barcelona
leaders:
  mribeirodantas:
    name: Marcel Ribeiro-Dantas
    slack: https://nfcore.slack.com/team/U03932BSX1V
  mashehu:
    name: Matthias Hörtenhuber
    slack: https://nfcore.slack.com/team/UQZG2UCF3
---

## Goal

Help nf-core pipelines validate and complete adoption of the current container setup: use Apptainer-compatible module directives, adopt Wave-built Seqera Containers where practical, and make sure installed and local modules select the right image for each engine and platform.

## Description

There are two related but distinct updates. Apptainer support was added alongside Singularity in the module template, and the shared modules were updated to recognize both engines ([migration PR](https://github.com/nf-core/modules/pull/11260)). Separately, nf-core/tools can use [Wave](https://seqera.io/wave/) to build Docker and Singularity-format images and Conda lock files from a module's `environment.yml`. The resulting module metadata can drive generated pipeline configuration for each supported engine and architecture.

The shared modules have already had the Apptainer-compatible directive change, so this project is not a repeat of that bulk migration. We'll focus on checking pipeline-local modules, adopting Wave containers in a manageable set of eligible modules, then installing or updating them in pipelines and verifying the generated configuration. The aim is to deliver tested pull requests and document gaps in tooling or adoption guidance—not to convert every module or pipeline during the hackathon.

Useful background:

- [nf-core/tools 4.0.0: Apptainer support and generated container configuration](https://nf-co.re/blog/2026/tools-4_0_0)
- [nf-core/tools 4.1.0: creating Seqera Containers with Wave](https://nf-co.re/blog/2026/tools-4_1_0)
- [Migrating nf-core modules to Seqera Containers](https://nf-co.re/blog/2024/seqera-containers-part-2)
- [nf-core module specifications](https://nf-co.re/docs/specifications/components/modules/general)
- [Seqera Containers](https://nf-co.re/docs/developing/containers/seqera-containers)

## What the current container directive looks like

Here, “new container syntax” means nf-core's module container selection and generated metadata; it is separate from Nextflow's strict-syntax migration.

Before the Wave migration, modules commonly selected a Docker or Singularity image inline in `main.nf`, for example:

```nextflow
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
    'biocontainers/fastqc:0.12.1--hdfd78af_0' }"
```

With Wave-built Seqera Containers, the generated declaration still selects an image based on the configured container engine: the Wave-built Docker image for Docker, or a Seqera-hosted Singularity-format image for Singularity/Apptainer. For example, the current [`nf-core/wget` module](https://github.com/nf-core/modules/blob/master/modules/nf-core/wget/main.nf) uses:

```nextflow
container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
    ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b54fa9135194c72a18d00db6b399c03248103f87e43ca75e4b50d61179994b3/data'
    : 'community.wave.seqera.io/library/wget:1.21.4--8b0fcde81c17be5e'}"
```

The `nf-core modules containers create` command generates the container declaration; do not edit it by hand. The 2024 migration blog describes a simpler, single-image declaration as the intended design, but the current nf-core/tools 4.1 implementation still generates this engine-conditional directive. It also records per-engine and per-architecture image references in `meta.yml`; when a module with this metadata is installed or updated, nf-core/tools generates the pipeline's container configuration. For example, metadata can map Docker images for both `linux/amd64` and `linux/arm64`:

```yaml
containers:
  docker:
    linux/amd64:
      name: community.wave.seqera.io/library/fastqc:0.12.1--5cfd0f3cb6760c42
    linux/arm64:
      name: community.wave.seqera.io/library/fastqc:0.12.1--d3caca66b4f3d3b0
```

See the full code [here](https://github.com/nf-core/modules/blob/37d472286462f856e1bb3ce9e865850b648eac29/modules/nf-core/fastqc/meta.yml#L88).

The container creation command is a beta preview in nf-core/tools 4.1. For an eligible module with a valid `environment.yml`, the starting point is:

```bash
nf-core modules containers create <module_name>
```

## Tasks

- Agree on the scope and claim modules or pipelines to avoid duplicated work (assign yourself in the issue).
- Check selected modules and local pipeline modules for the current container directives; update them where needed.
- For suitable modules that do not yet have Wave container metadata and have a valid `environment.yml`, try `nf-core modules containers create <module_name>` to build Wave containers and update the module metadata and container definition.
- Install or update converted modules in a pipeline, then review the generated per-platform container configuration rather than editing those generated files by hand.
- Run the relevant module and pipeline checks, and open pull requests with the changes and any migration caveats.
- Record blockers or missing functionality and, where useful, open follow-up issues or pull requests against [nf-core/tools](https://github.com/nf-core/tools).

:::info
Creating Wave containers uses the Wave API. A Seqera Platform access token can help avoid API rate limits; do not include tokens in pull requests or logs.
:::

## Making PRs

Module changes should go to [nf-core/modules](https://github.com/nf-core/modules), and pipeline changes (for local modules, for example) should go to the relevant pipeline's `dev` branch. Tooling or template changes should go to [nf-core/tools](https://github.com/nf-core/tools). Ask for a review in the relevant project or repository Slack channel before merging.
