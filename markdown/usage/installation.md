---
title: Installation
subtitle: Installing the software requirements needed for running nf-core pipelines.
---

## Nextflow

All nf-core pipelines use Nextflow, so this must be present on the system where you launch your analysis.
See [nextflow.io](https://www.nextflow.io/docs/latest/getstarted.html#installation) for the latest installation instructions. Once installed you will probably need to configure Nextflow to run on your system. For instructions, see [_Nextflow configuration_](configuration.md).

### Typical installation

Generally speaking, Nextflow runs on most POSIX systems (Linux, Mac OSX etc) and can typically be installed by running the following commands:

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

First, set up Bioconda according to the [Bioconda documentation](https://bioconda.github.io/user/install.html), notably setting up channels:

```console
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Then install Nextflow:

```console
conda install nextflow
```

### Edge releases

Stable releases will be becoming more infrequent as Nextflow shifts its development model to becoming more dynamic via the usage of plugins. This will allow functionality to be added as an extension to the core codebase with a release cycle that could potentially be independent to that of Nextflow itself. As a result of the reduction in stable releases, some pipelines may be required to use Nextflow `edge` releases in order to be able to exploit cutting "edge" features e.g. version 3.0 of the nf-core/rnaseq pipeline requires Nextflow `>=20.11.0-edge` in order to be able to directly download Singularity containers over `http` (see [nf-core/rnaseq#496](https://github.com/nf-core/rnaseq/issues/496)).

There are a number of ways you can install Nextflow `edge` releases, the main difference with stable releases being that you have to `export` the version you would like to install before issuing the appropriate installation/execution commands as highlighted below.

* If you have Nextflow installed already, you can issue the version you would like to use on the same line as the pipeline command and it will be fetched if required before the pipeline execution.

```bash
NXF_VER="20.11.0-edge" nextflow run nf-core/rnaseq -profile test,docker -r 3.0
```

* If you have Nextflow installed already, another alternative to the option above is to `export` it as an environment variable before you run the pipeline command:

```bash
export NXF_VER="20.11.0-edge"
nextflow run nf-core/rnaseq -profile test,docker -r 3.0
```

* If you would like to download and install a Nextflow `edge` release from scratch with minimal fuss:

```bash
export NXF_VER="20.11.0-edge"
wget -qO- get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow run nf-core/rnaseq -profile test,docker -r 3.0
```

> Note if you don't have `sudo` privileges required for the last command above then you can move the `nextflow` binary to somewhere else and export that directory to `$PATH` instead. One way of doing that on Linux would be to add `export PATH=$PATH:/path/to/nextflow/binary/` to your `~/.bashrc` file so that it is available every time you login to your system.

* Manually download and install Nextflow from the available [assets](https://github.com/nextflow-io/nextflow/releases) on Github. See [Nextflow installation docs](https://www.nextflow.io/docs/latest/getstarted.html#installation).

### Updating

We recommend using a personal installation of Nextflow where possible, instead of using a system-wide installation. This makes it easier to update.

Updating nextflow is as simple as running `nextflow self-update`
or `conda update nextflow`, depending on how it was installed.

## Pipeline software

An analysis pipeline chains the execution of multiple tools together.
Historically, all tools would have to be manually installed, often a source of great frustration and a key step where reproducibility between analyses is lost.
nf-core pipelines utilise the built-in support for software packaging that Nextflow offers: all can work with Docker and Singularity, and most pipelines also have support for Conda.

* [Docker](https://docs.docker.com/install/)
  * Typically used locally / on single-user servers and the cloud.
  * Analysis runs in a _container_, which behaves like an isolated operating system
  * Previously required system root access, though a "rootless mode" is available since late 2019
* [Singularity](https://www.sylabs.io/)
  * Often used as an alternative to Docker on multi-user systems such as HPC systems.
  * Also runs _containers_ and can create these from Docker images
  * Does not need root access or any daemon processes - images built from files
* [Podman](https://podman.io/)
  * Often used as an alternative to Docker on multi-user systems such as HPC systems.
  * Is a daemonless container engine that can serve as a drop-in replacement for Docker.
* [Charliecloud](https://hpc.github.io/charliecloud/)
  * Often used as an alternative to Docker on multi-user systems such as HPC systems.
  * Uses Linux user namespaces to run containers with no privileged operations or daemons
* [Shifter](https://www.nersc.gov/research-and-development/user-defined-images/)
  * An experimental implementation of container system that can convert from a wide range of other images
* [Conda](https://conda.io/)
  * Packaging system that manages environments instead of running analysis in containers.
  * Poorer reproducibility than Docker / Singularity
    * There can be changes in low-level package dependencies over time
    * The software still runs in your native operating system environment and so core system functions can differ

## Pipeline code

### Automatic

This pipeline itself needs no installation - nextflow will automatically fetch it from GitHub if `nf-core/<pipeline-name>` is specified as the pipeline name.

This method requires an internet connection. If you're running on a system that has no internet connection, please see the [Running Offline](offline.md) documentation.

### Development

If you would like to make changes to the pipeline, fork the GitHub repository and then clone the files. Once cloned you can run the pipeline with `nextflow run <path-to-repo>`.

Note that you should _only_ do this if you intend to make significant changes to the pipeline. All configuration options can be changed without editing the pipeline code. Forking the pipeline repositories means that you cannot use stable releases and you will fall behind new updates.

## Reference genomes

Some pipelines come with built-in support for iGenomes references.
It may be preferable for you to download a local copy of these to your system to avoid fetching the same reference many times.
For more information, see [_Reference genomes_](reference_genomes.md).
