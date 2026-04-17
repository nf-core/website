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
- **Vendor-provided GPU containers** with no conda equivalent: use the vendor container directly.
  For example, [`parabricks/rnafq2bam`](https://github.com/nf-core/modules/tree/master/modules/nf-core/parabricks/rnafq2bam) uses NVIDIA's container.
  These modules SHOULD guard against conda/mamba profiles:

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

### CUDA version pinning

The `environment.gpu.yml` SHOULD pin the CUDA major version to avoid the conda solver selecting builds for unreleased CUDA versions:

```yaml
dependencies:
  - "bioconda::ribodetector=0.3.3"
  - "conda-forge::pytorch-gpu=2.10.0"
  - "conda-forge::cuda-version>=12,<13"
```

NVIDIA drivers are backward compatible, so a container built against CUDA 12.x will run on any CUDA 12.0+ host driver.
A CUDA 12.x container will not run on an older CUDA 11.x driver, so modules that need to support older hosts MAY provide an alternative `environment.gpu.yml` pinned to CUDA 11.

### Script patterns

Tools that provide separate GPU and CPU binaries SHOULD select between them based on `task.accelerator`. For example, [`ribodetector`](https://github.com/nf-core/modules/blob/master/modules/nf-core/ribodetector/main.nf):

```groovy
def binary = task.accelerator ? "ribodetector" : "ribodetector_cpu"
```

Tools that accept a GPU count SHOULD read it from `task.accelerator.request`, allowing users to override via their pipeline config (e.g., `accelerator = 2`). For example, [`parabricks/rnafq2bam`](https://github.com/nf-core/modules/blob/master/modules/nf-core/parabricks/rnafq2bam/main.nf):

```groovy
def num_gpus = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''
```

## Software not on Bioconda

If the software is not available on Bioconda a `Dockerfile` MUST be provided within the module directory. nf-core will use GitHub Actions to auto-build the containers on the [GitHub Packages registry](https://github.com/features/packages).
