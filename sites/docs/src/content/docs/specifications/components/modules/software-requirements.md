---
title: Software requirements
subtitle: Specify software dependencies
markdownPlugin: addNumbersToHeadings
shortTitle: Software
weight: 7
---

The keywords "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

[BioContainers](https://biocontainers.pro/#/) is a registry of Docker and Singularity containers automatically created from all of the software packages on [Bioconda](https://bioconda.github.io/).
Where possible nf-core will use BioContainers to fetch pre-built software containers and Bioconda to install software using Conda.

## Use of container directives

Software requirements SHOULD be declared within the module file using the Nextflow `container` directive.
For single-tool BioContainers, the `nf-core modules create` command will automatically fetch and fill in the appropriate Conda / Docker / Singularity definitions by parsing the information provided in the first part of the module name:

```groovy
conda "bioconda::fastqc=0.11.9"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
    'biocontainers/fastqc:0.11.9--0' }"
```

## Use of conda directive

If the software is available on Conda, also define it in an `environment.yml` file alongside the `main.nf` of the module, and pass it to the Nextflow `conda` directive within `main.nf`.

Using `bioconda::bwa=0.7.17` as an example, pin software to the channel (i.e., `bioconda`) and version (i.e., `0.7.17`).

Do not pin Conda packages to a build because they can vary on different platforms.

## Re-use of multi-tool containers

Multi-tool containers are also available on BioContainers.
For example, [`bwa` and `samtools`](https://biocontainers.pro/#/tools/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40).
Install and use the [`galaxy-tool-util`](https://anaconda.org/bioconda/galaxy-tool-util) package to search for both single- and multi-tool containers available in Conda, Docker and Singularity format.
For example, to search for Docker (hosted on Quay.io) and Singularity multi-tool containers with both `bowtie` and `samtools` installed, use the following command:

```console
mulled-search --destination quay singularity --channel bioconda --search bowtie samtools | grep "mulled"
```

:::note
Build information for all tools within a multi-tool container can be obtained in the `/usr/local/conda-meta/history` file within the container.
:::

## Creation of new multi-tool containers

It is also possible for a new multi-tool container to be built and added to BioContainers by submitting a pull request on their [`multi-package-containers`](https://github.com/BioContainers/multi-package-containers) repository.

- Fork the [multi-package-containers repository](https://github.com/BioContainers/multi-package-containers)
- Make a change to the `hash.tsv` file in the `combinations` directory see [here](https://github.com/aunderwo/multi-package-containers/blob/master/combinations/hash.tsv#L124) for an example where `pysam=0.16.0.1,biopython=1.78` was added.
- Commit the code and then make a pull request to the original repo, for [example](https://github.com/BioContainers/multi-package-containers/pull/1661)
- Once the PR has been accepted a container will get built and you can find it using a search tool in the `galaxy-tool-util conda` package

  ```console
  mulled-search --destination quay singularity conda  --search pysam biopython  | grep "mulled"
  quay         mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f  185a25ca79923df85b58f42deb48f5ac4481e91f-0  docker pull quay.io/biocontainers/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0
  singularity  mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f  185a25ca79923df85b58f42deb48f5ac4481e91f-0  wget https://depot.galaxyproject.org/singularity/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0
  ```

- Copy and paste the `mulled-*` path into the relevant Docker and Singularity lines in the Nextflow `process` definition of the module
- To confirm that this is correct, spin up a temporary Docker container:

  ```console
  docker run --rm -it quay.io/biocontainers/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0  /bin/sh
  ```

  And in the command prompt type:

  ```console
  $ grep specs /usr/local/conda-meta/history
  # update specs: ['biopython=1.78', 'pysam=0.16.0.1']
  ```

  The packages should reflect those added to the multi-package-containers repo `hash.tsv` file.

- If the multi-tool container already exists, use [this](https://midnighter.github.io/mulled) helper tool to obtain the `mulled-*` path.

## GPU-capable modules

GPU-enabled software has two properties that make it awkward to package the same way as CPU-only tools: GPU builds (e.g., CUDA PyTorch) can be several GB larger than their CPU counterparts, and some vendor containers (e.g., NVIDIA Parabricks) are proprietary with no conda equivalent.
The specification therefore allows three container approaches, chosen according to the tool:

- **Dual CPU/GPU variants of a tool**, where the GPU build has significant overhead (e.g., CUDA PyTorch adds ~3 GB): use the [dual-container pattern](#dual-container-pattern) below so CPU-only users are not penalised.
  For example, [`ribodetector`](https://github.com/nf-core/modules/tree/master/modules/nf-core/ribodetector).
- **Minimal GPU overhead or CPU fallback within one container**: a single container is simpler and preferred.
- **Vendor-provided GPU containers** with no conda equivalent.
  For example, [`parabricks/rnafq2bam`](https://github.com/nf-core/modules/tree/master/modules/nf-core/parabricks/rnafq2bam) uses NVIDIA's container.
  Vendor-provided container modules SHOULD guard against conda/mamba profiles:

  ```groovy
  if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
      error("This module does not support Conda. Please use Docker / Singularity / Podman instead.")
  }
  ```

### Dual-container pattern

When the GPU container is substantially larger, modules SHOULD switch between containers based on `task.accelerator`:

```groovy
conda "${ task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml" }"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    (task.accelerator ? '<singularity-gpu-url>' : '<singularity-cpu-url>') :
    (task.accelerator ? '<docker-gpu-url>' : '<docker-cpu-url>') }"
```

A separate `environment.gpu.yml` SHOULD be provided for GPU-specific dependencies.
The CPU `environment.yml` MUST remain unchanged so that non-GPU users are unaffected.

GPU containers SHOULD be built using [Wave](https://wave.seqera.io) from the `environment.gpu.yml` file.
Both Docker and Singularity URLs MUST be provided.

:::note{title="Wave build template"}
Until it becomes the default, GPU environments that pull packages requiring the `__cuda` virtual package (most post-2022 CUDA-aware conda builds) MUST be built with `--build-template conda/micromamba:v2` so Wave's [automatic `CONDA_OVERRIDE_CUDA` retry](https://github.com/seqeralabs/wave/pull/1027) kicks in. Solves exercising the retry can exceed the default 15-minute `--await`; pass `--await 60m` to cover them.
:::

### CUDA version pinning

Pin `cuda-version` exactly (nf-core does not allow version ranges). The pin sets the minimum NVIDIA driver version required on hosts, so pick the lowest value the GPU package actually supports on conda-forge - that gives the widest host compatibility.

```yaml
dependencies:
  - "bioconda::ribodetector=0.3.3"
  - "conda-forge::pytorch-gpu=2.1.0"
  - "conda-forge::cuda-version=11.2"
```

Use `micromamba search -c conda-forge '<package>=<version>'` to see which CUDA minor versions the GPU package was built against on conda-forge - this is usually the real floor. For instance `pytorch-gpu=2.1.0` has `cuda112` builds, while `pytorch-gpu=2.10.0` only has `cuda128`/`cuda129`/`cuda130` builds, raising the driver floor accordingly.

NVIDIA drivers are backward compatible with older CUDA versions, so pinning lower widens host reach. The reverse is not true: a CUDA 12.x container cannot run on a host with a CUDA 11.x-only driver. Modules that specifically want broader reach on legacy hosts MAY provide an alternative `environment.gpu.yml` pinned to CUDA 11.

### Capturing the CUDA runtime version

GPU modules SHOULD emit the CUDA runtime version on the `versions` topic channel so it appears in provenance reports alongside the tool version. One simple pattern, using the pytorch dependency that most CUDA-aware conda envs already pull in:

```groovy
tuple val("${task.process}"), val('cuda'), eval('python -c "import torch; print(torch.version.cuda or \'cpu\')"'), emit: versions_cuda, topic: versions
```

This reports the actual CUDA minor the container was built with on the GPU path, and `cpu` on the non-GPU path of dual-container modules.

### Pip-based GPU packages

Some GPU-compiled Python tools ship only as pre-built wheels from custom pip indexes rather than conda (e.g. [`llama-cpp-python`](https://abetlen.github.io/llama-cpp-python/whl/)). Pin the full wheel URL in the `pip:` block and pull the CUDA runtime from conda-forge alongside:

```yaml
dependencies:
  - python=3.11
  - pip
  - "conda-forge::cuda-version=12.4"
  - "conda-forge::cuda-runtime"
  - pip:
      - "https://github.com/<owner>/<project>/releases/download/v<version>-cu124/<wheel>.whl"
```

Build the container with Wave's `--config-env` so the wheel's binary can resolve the conda-provided CUDA libs at `dlopen` time:

```bash
wave --conda-file environment.gpu.yml --freeze --await --singularity \
     --config-env 'LD_LIBRARY_PATH=/opt/conda/lib'
```

Conda's `activate.d` hooks don't fire under `docker run`, so the library path has to be set at image-build time.
Conda-forge-native GPU packages (e.g. `pytorch-gpu`) ship RPATHs in their binaries and don't need this.

### Binary and GPU count selection in scripts

Tools that provide separate GPU and CPU binaries SHOULD select between them based on `task.accelerator`.
For example, [`ribodetector`](https://github.com/nf-core/modules/blob/20423f58f6ff54c4bc851dfd143d75f9b9f86f41/modules/nf-core/ribodetector/main.nf):

```groovy
def binary = task.accelerator ? "ribodetector" : "ribodetector_cpu"
```

Tools that accept a GPU count SHOULD specify this in the command using `task.accelerator.request`, allowing users to override via their pipeline config (e.g., `accelerator = 2`).
This parameter SHOULD NOT be hardcoded.

For example, [`parabricks/rnafq2bam`](https://github.com/nf-core/modules/blob/0c44f69eefe8a8c373c7cdd9528b3e1d60cb895f/modules/nf-core/parabricks/rnafq2bam/main.nf):

```groovy
def num_gpus = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''
```

## Software not on Bioconda

If the software is not available on Bioconda a `Dockerfile` MUST be provided within the module directory. nf-core will use GitHub Actions to auto-build the containers on the [GitHub Packages registry](https://github.com/features/packages).
