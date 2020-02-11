---
title: Installation
subtitle: Installing the software requirements needed for running nf-core pipelines.
---

This documentation assumes that you have already read the [introduction](introduction.md) and are familiar with the tools described below.

## Table of contents
<!-- Install Atom plugin markdown-toc-auto for this ToC -->
<!-- TOC START min:2 max:3 link:true asterisk:true -->
* [Install nextflow](#install-nextflow)
* [Install the pipeline](#install-the-pipeline)
  * [Automatic](#automatic)
  * [Offline](#offline)
  * [Development](#development)
* [Pipeline configuration](#pipeline-configuration)
  * [Docker](#docker)
  * [Singularity](#singularity)
  * [Conda](#conda)
  * [Configuration profiles](#configuration-profiles)
* [Reference genomes](#reference-genomes)
<!-- TOC END -->

## Install nextflow
Nextflow runs on most POSIX systems (Linux, Mac OSX etc). It can be installed by running the following commands:
=======
## Nextflow

All nf-core pipelines use Nextflow, so this must be present on the system where you launch your analysis.
See [nextflow.io](https://www.nextflow.io/) for the latest installation instructions.

Generally speaking, Nextflow runs on most POSIX systems (Linux, Mac OSX etc) and can typically be installed by running the following commands:
>>>>>>> upstream/master

```bash
# Make sure that Java v8+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your user's PATH:
mv nextflow ~/bin/
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
```

You can also install Nextflow using [Bioconda](https://bioconda.github.io/):

```bash
conda install -c bioconda nextflow
```

We recommend using a personal installation of Nextflow where possible, instead of using a system-wide installation. This makes it easier to update.

Updating nextflow is as simple as running `nextflow self-update`
or `conda update nextflow`, depending on how it was installed.

Once installed you will probably need to configure Nextflow to run on your system. For instructions, see [_Nextflow configuration_](configuration.md).


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
* [Conda](https://conda.io/)
  * Packaging system that manages environments instead of running analysis in containers.
  * Poorer reproducibility than Docker / Singularity
    * There can be changes in low-level package dependencies over time
    * The software still runs in your native operating system environment and so core system functions can differ

## Pipeline code

### Automatic

This pipeline itself needs no installation - nextflow will automatically fetch it from GitHub if `nf-core/<pipeline-name>` is specified as the pipeline name.

This method requires an internet connection. If you're running on a system that has no internet connection, please see the [Running Offline](/usage/offline) documentation.

### Development

If you would like to make changes to the pipeline, fork the GitHub repository and then clone the files. Once cloned you can run the pipeline with `nextflow run <path-to-repo>`.

<<<<<<< HEAD
See [Adding your own configuration profile](adding_own_config.md).

## Reference genomes

See [Reference genomes](reference_genomes.md)
=======
Note that you should _only_ do this if you intend to make significant changes to the pipeline. All configuration options can be changed without editing the pipeline code. Forking the pipeline repositories means that you cannot use stable releases and you will fall behind new updates.

## Reference genomes

Some pipelines come with built-in support for iGenomes references.
It may be preferable for you to download a local copy of these to your system to avoid fetching the same reference many times.
For more information, see [_Reference genomes_](reference_genomes).
>>>>>>> upstream/master
