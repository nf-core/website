---
title: GPU-capable modules
subtitle: Write modules that use GPU acceleration
shortTitle: GPU modules
---

Some bioinformatics tools run substantially faster on GPUs, but GPU-enabled software is awkward to package the same way as CPU-only tools.
GPU builds (e.g., CUDA PyTorch) can be several gigabytes larger than their CPU counterparts, and some vendor containers (e.g., NVIDIA Parabricks) are proprietary with no conda equivalent.

This guide explains how to write nf-core modules that use GPU acceleration: detecting GPU requests, choosing a container strategy, pinning CUDA versions, and testing.

For the normative requirements, see the [resource requirements](/docs/specifications/components/modules/resource-requirements#gpu-acceleration), [software requirements](/docs/specifications/components/modules/software-requirements#gpu-capable-modules), and [testing](/docs/specifications/components/modules/testing#gpu-tests) specifications.

## Detecting a GPU request

Modules detect whether a GPU has been requested with `task.accelerator{:groovy}`.
A module never sets the `accelerator` directive itself; the pipeline controls GPU allocation by setting `accelerator = 1{:groovy}` in its process config (for example, via a `process_gpu` label or a `withName` block).

:::info{title="Why detection lives in the module, not a label"}
Placing GPU allocation in the pipeline config lets users control it through their pipeline config or profiles.
A label-only alternative (e.g., requiring a `process_gpu` label) would not work for modules that support both CPU and GPU modes, such as [`ribodetector`](https://github.com/nf-core/modules/tree/master/modules/nf-core/ribodetector), so the choice is left to the pipeline author.
:::

Pipelines also set GPU container flags via `containerOptions` in their process config.
Use `containerOptions` (not global `docker.runOptions`) to scope GPU flags to GPU processes only.

## Choosing a container approach

Pick a container strategy based on the tool:

- **Minimal GPU overhead, or a CPU fallback that fits in one image.** A single container is simpler and preferred.
- **Dual CPU/GPU variants where the GPU build has significant overhead** (e.g., CUDA PyTorch adds ~3 GB). Use the [dual-container pattern](#the-dual-container-pattern) so CPU-only users are not penalised. For example, [`ribodetector`](https://github.com/nf-core/modules/tree/master/modules/nf-core/ribodetector).
- **Vendor-provided GPU containers** with no conda equivalent. For example, [`parabricks/rnafq2bam`](https://github.com/nf-core/modules/tree/master/modules/nf-core/parabricks/rnafq2bam) uses NVIDIA's container. Guard these modules against conda/mamba profiles:

  ```groovy
  if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
      error("This module does not support Conda. Please use Docker / Singularity / Podman instead.")
  }
  ```

## The dual-container pattern

When the GPU container is substantially larger, switch between containers based on `task.accelerator`:

```groovy
conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    (task.accelerator ? '<singularity-gpu-url>' : '<singularity-cpu-url>') :
    (task.accelerator ? '<docker-gpu-url>' : '<docker-cpu-url>') }"
```

Provide a separate `environment.gpu.yml` for the GPU-specific dependencies, and leave the CPU `environment.yml` unchanged so that non-GPU users are unaffected.
Build the GPU container from `environment.gpu.yml` using [Wave](https://wave.seqera.io), and provide both Docker and Singularity URLs.

:::note{title="Wave build template"}
Pass `--build-template conda/micromamba:v2` to Wave when building GPU environments. This is required for now, until it becomes the default.
:::

## Pinning the CUDA version

Pin `cuda-version` exactly (nf-core does not allow version ranges).
The pin sets the minimum NVIDIA driver version required on hosts, so pick the lowest value the GPU package actually supports on conda-forge — that gives the widest host compatibility.

```yaml
dependencies:
  - "bioconda::ribodetector=0.3.3"
  - "conda-forge::pytorch-gpu=2.1.0"
  - "conda-forge::cuda-version=11.2"
```

Use `micromamba search -c conda-forge '<package>=<version>'` to see which CUDA minor versions the GPU package was built against on conda-forge — this is usually the real floor.
For instance, `pytorch-gpu=2.1.0` has `cuda112` builds, while `pytorch-gpu=2.10.0` only has `cuda128`/`cuda129`/`cuda130` builds, raising the driver floor accordingly.

NVIDIA drivers are backward compatible with older CUDA versions, so pinning lower widens host reach.
The reverse is not true: a CUDA 12.x container cannot run on a host with a CUDA 11.x-only driver.
Modules that specifically want broader reach on legacy hosts MAY provide an alternative `environment.gpu.yml` pinned to CUDA 11.

## Capturing the CUDA runtime version

Emit the CUDA runtime version on the `versions` topic channel so it appears in provenance reports alongside the tool version.
One simple pattern uses the pytorch dependency that most CUDA-aware conda environments already pull in:

```groovy
tuple val("${task.process}"), val('cuda'), eval('python -c "import torch; print(torch.version.cuda or \'no CUDA available\')"'), emit: versions_cuda, topic: versions
```

This reports the actual CUDA minor the container was built with on the GPU path, and `no CUDA available` on the non-GPU path of dual-container modules.
Prefer a descriptive string over something like `cpu`, which reviewers reasonably flag as not a version.

## GPU-enabled Python packages

A Python tool may delegate GPU work to a separately-packaged backend whose CUDA build is selected on conda-forge with a build-string match.
Pin both the wrapper and the backend version, and select the CUDA build of the backend with a build glob (for example, [`llama-cpp-python`](https://github.com/abetlen/llama-cpp-python) over the `llama.cpp` backend):

```yaml
dependencies:
  - "conda-forge::llama-cpp-python=0.3.16"
  - "conda-forge::llama.cpp=6191=*cuda*"
```

The `=6191=*cuda*` spec pins the backend version and matches a CUDA build of it (the `name=version=build` form is required to constrain the build string).
It pulls in the corresponding `cuda-version`, `libcublas` and `cuda-cudart` packages and carries RPATHs in its binaries, so it needs no extra runtime configuration.
The build glob targets the CUDA variant rather than a single build hash; pin `cuda-version` as well (see [Pinning the CUDA version](#pinning-the-cuda-version)) to fix the CUDA minor and the resulting driver floor.
`environment.gpu.yml` requests the `=*cuda*` build while a plain `environment.yml` requests the `=*cpu*` build of the same pinned version.

The build glob is not optional: without it the solver selects the CPU backend of the same version and silently produces a non-accelerated container from `environment.gpu.yml`.
The explicit `=*cuda*` match is what pins the variant to the CUDA backend.
Building a CUDA backend with Wave requires the v2 build template (`--build-template conda/micromamba:v2`), because the stock build image cannot resolve the `__cuda`-gated CUDA packages.

When a tool has no conda packaging and is distributed only as a pre-built wheel from a custom pip index, pin the full wheel URL in a `pip:` block and pull the CUDA runtime from conda-forge alongside:

```yaml
dependencies:
  - python=3.11
  - pip
  - "conda-forge::cuda-version=12.4"
  - "conda-forge::cuda-runtime"
  - pip:
      - "https://github.com/<owner>/<project>/releases/download/v<version>-cu124/<wheel>.whl"
```

Build the container with Wave's `--config-env` so the wheel's binary can resolve the conda-provided CUDA libraries at `dlopen` time:

```bash
wave --conda-file environment.gpu.yml --freeze --await --singularity \
     --config-env 'LD_LIBRARY_PATH=/opt/conda/lib'
```

`activate.d` hooks do not fire under `docker run`, so the library path is set at image-build time.
Backends installed from conda-forge ship RPATHs in their binaries and do not need this.

## Selecting binaries and GPU count in scripts

Tools that provide separate GPU and CPU binaries select between them based on `task.accelerator`.
For example, [`ribodetector`](https://github.com/nf-core/modules/blob/20423f58f6ff54c4bc851dfd143d75f9b9f86f41/modules/nf-core/ribodetector/main.nf):

```groovy
def binary = task.accelerator ? "ribodetector" : "ribodetector_cpu"
```

Tools that accept a GPU count specify it in the command using `task.accelerator.request`, allowing users to override via their pipeline config (e.g., `accelerator = 2`).
Do not hardcode this value.
For example, [`parabricks/rnafq2bam`](https://github.com/nf-core/modules/blob/0c44f69eefe8a8c373c7cdd9528b3e1d60cb895f/modules/nf-core/parabricks/rnafq2bam/main.nf):

```groovy
def num_gpus = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''
```

## Testing GPU modules

Modules that support both CPU and GPU modes include a separate GPU test file (`main.gpu.nf.test`).
GPU-only modules may use a single test file — see [`parabricks`](https://github.com/nf-core/modules/tree/master/modules/nf-core/parabricks) for an example.

Tag GPU tests with `"gpu"` or `"gpu_highmem"` so the GPU CI workflow discovers and runs them on GPU-enabled runners.

:::note
The `"gpu"` tag runs on smaller AWS GPU instances (e.g., [`g4dn.xlarge`](https://aws.amazon.com/ec2/instance-types/g4/)), while `"gpu_highmem"` runs on larger instances (e.g., [`g4dn.2xlarge`](https://aws.amazon.com/ec2/instance-types/g4/)) for tools with higher memory requirements such as [Parabricks](https://github.com/nf-core/modules/tree/master/modules/nf-core/parabricks).
:::

Include a `nextflow.gpu.config` that sets `accelerator = 1` on the process, use the same assertions as the CPU tests so that GPU and CPU modes are verified to produce equivalent results, and include both a real test and a stub test.

```
modules/nf-core/<tool>/tests/
  main.nf.test           # CPU tests
  main.gpu.nf.test       # GPU tests (tag "gpu")
  nextflow.gpu.config    # Sets accelerator = 1
```

For an example, see the [`ribodetector` GPU tests](https://github.com/nf-core/modules/tree/master/modules/nf-core/ribodetector/tests).

:::caution{title="GPU concurrency under Singularity"}
Multiple concurrent GPU processes sharing a single GPU can deadlock under Singularity.
Docker's NVIDIA runtime handles GPU memory arbitration, but Singularity does not.
This can happen in CI (where all tasks share one GPU), on HPC nodes with a local executor, or any setup where multiple GPU tasks land on the same machine.
Set `maxForks = 1` on GPU processes to serialise access when this is a risk.
:::
