---
title: Installation
subtitle: Install the software requirements needed to run nf-core pipelines.
menu:
  main:
    weight: 20
---

## Nextflow

All nf-core pipelines use Nextflow, so this must be installed on the system where you launch your analysis.
See [nextflow.io](https://www.nextflow.io/docs/latest/getstarted.html#installation) for the latest installation instructions. Once installed, see [_Nextflow configuration_](configuration.md) to set up Nextflow to run on your system.

### Typical installation

Nextflow runs on most POSIX systems (Linux, macOS, etc) and can typically be installed by running these commands:

```console
# Make sure that Java v8+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your user's PATH:
mv nextflow ~/bin/
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
```

### Bioconda installation

You can also install Nextflow using [Bioconda](https://bioconda.github.io/).

First, set up Bioconda according to the [Bioconda documentation](https://bioconda.github.io/#usage), notably setting up channels:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

> **Warning**
> This is important - if not done, dependencies (such as Java) will be installed from
> the wrong channels and things may break in strange ways.

A best practice with conda is to create an environment and install the tools in it.
Therefore you will prevent version conflicts and keep everything clean.
To do so use the following command:

```bash
conda create --name env_nf nextflow
conda activate env_nf
```

To deactivate the conda environment, run the following command:

```bash
conda deactivate
```

> If you're already in the conda environment you want to use, you can just install Nextflow directly:
>
> ```bash
> conda install nextflow
> ```

If you want to develop your own pipelines or collaborate with others using the nf-core guidelines you will need the nf-core packages.
Please follow the instructions [here](tutorials/nf_core_usage_tutorial.md)

### Particularity for Windows

For Windows the installation procedure is more complex and is fully described on the [Nextflow website](https://nextflow.io/blog/2021/setup-nextflow-on-windows.html).

The main steps will be the following:

- Install Windows PowerShell
- Configure the Windows Subsystem for Linux (WSL2)
- Install a Linux distribution (on WSL2)

The step to install Nextflow in itself will afterwards be the same as previously mentioned.

### Edge releases

Stable releases will become less frequent as Nextflow shifts to a more dynamic development model with the use of plugins. This will allow functionality to be added as an extension to the core codebase with a release cycle that is potentially independent of Nextflow itself. As a result of the reduction in stable releases, some pipelines may require the use of Nextflow `edge` releases in order to exploit cutting "edge" features. For example, version 3.0 of the nf-core/rnaseq pipeline requires Nextflow `>=20.11.0-edge`, which allows Singularity containers to be downloaded over `http` (see [nf-core/rnaseq#496](https://github.com/nf-core/rnaseq/issues/496)).

There are several ways to install Nextflow `edge` releases. The main difference to stable releases is that you have to `export` the version you want to install before running the installation/execution commands highlighted below.

- If you have Nextflow installed already, you can flag the version to use in the pipeline command and it will be fetched (if required) before pipeline execution.

```bash
NXF_VER="20.11.0-edge" nextflow run nf-core/rnaseq -profile test,docker -r 3.0 --outdir <OUTDIR>
```

- If you have Nextflow installed already, another option is to `export` it as an environment variable before you run the pipeline command:

```bash
export NXF_VER="20.11.0-edge"
nextflow run nf-core/rnaseq -profile test,docker -r 3.0 --outdir <OUTDIR>
```

- If you would like to download and install a Nextflow `edge` release from scratch with minimal fuss:

```bash
export NXF_VER="20.11.0-edge"
wget -qO- get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow run nf-core/rnaseq -profile test,docker -r 3.0 --outdir <OUTDIR>
```

> Note: if you don't have the `sudo` privileges required for the command above, you can move the `nextflow` binary to another directory and export that directory to `$PATH` instead. One way of doing that on Linux would be to add `export PATH=$PATH:/path/to/nextflow/binary/` to your `~/.bashrc` file so that it is available every time you log into your system.

- Manually download and install Nextflow from the available [assets](https://github.com/nextflow-io/nextflow/releases) on Github. See the [Nextflow installation docs](https://www.nextflow.io/docs/latest/getstarted.html#installation).

### Updating Nextflow

We recommend using a personal installation of Nextflow where possible, instead of using a system-wide installation. This makes it easier to update.

Updating nextflow is as simple as running `nextflow self-update`
or `conda update nextflow`, depending on how it was installed.

## Pipeline software

An analysis pipeline chains the execution of multiple tools together.
Historically, all tools would have to be manually installed — often a source of great frustration and a key step where reproducibility between analyses is lost.
nf-core pipelines utilise the built-in support for software packaging that Nextflow offers: all can work with Docker and Singularity, and most pipelines also support Conda.

- [Docker](https://docs.docker.com/install/)
  - Typically used locally, on single-user servers, and the cloud
  - Analysis runs in a _container_, which behaves like an isolated operating system
  - Previously required system root access, though a "rootless mode" has been available since 2019
- [Singularity](https://www.sylabs.io/)
  - Often used as an alternative to Docker on multi-user systems, such as HPC systems
  - Also runs _containers_ and can create these from Docker images
  - Does not need root access or any daemon processes - images are built from files
- [Podman](https://podman.io/)
  - Often used as an alternative to Docker on multi-user systems, such as HPC systems
  - Is a daemonless container engine that can serve as a drop-in replacement for Docker
- [Charliecloud](https://hpc.github.io/charliecloud/)
  - Often used as an alternative to Docker on multi-user systems, such as HPC systems
  - Uses Linux user namespaces to run containers with no privileged operations or daemons
- [Shifter](https://www.nersc.gov/research-and-development/user-defined-images/)
  - An experimental implementation of a container system that can convert from a wide range of other images
- [Conda](https://conda.io/)
  - Packaging system that manages environments instead of running analysis in containers
  - Poorer reproducibility than Docker / Singularity
    - There can be changes in low-level package dependencies over time
    - The software still runs in your native operating system environment and so core system functions can differ

## Pipeline code

### Automatic

The pipeline needs no installation - Nextflow will automatically fetch it from GitHub if `nf-core/<pipeline-name>` is specified as the pipeline name.

This method requires an internet connection. If you're running on a system that has no internet connection, please see [Running Offline](offline.md).

### Development

If you would like to make changes to the pipeline, fork the GitHub repository and then clone the files. Once cloned, you can run the pipeline with `nextflow run <path-to-repo>`.

Note that you should _only_ do this if you intend to make significant changes to the pipeline. All configuration options can be changed without editing the pipeline code. Forking the pipeline repositories means that you cannot use stable releases and you will fall behind new updates.

## Reference genomes

Some pipelines come with built-in support for iGenomes references.
We recommend downloading these references locally to avoid fetching the same reference many times.
For more information, see [_Reference genomes_](reference_genomes.md).
