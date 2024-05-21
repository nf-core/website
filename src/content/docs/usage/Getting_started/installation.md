---
title: Dependency installation
subtitle: Install the software dependencies to run nf-core pipelines.
shortTitle: Dependency installation
weight: 2
---

## Nextflow

All nf-core pipelines require Nextflow to be installed on the system where you launch your analysis.

Use a personal installation of Nextflow where possible, instead of a system-wide installation.
This makes it easier to update and control versions.

To supplement the official [Nextflow installation instructions](https://www.nextflow.io/docs/latest/install.html), this page describes '[quick start](#quick-start-installation)', [Conda](#bioconda-installation), and [Windows](#installation-on-windows) installation instructions.

Any instructions on this page are provided for convinence, and may not be up-to-date.

:::note
You only need Nextflow to run nf-core pipelines.
However, the `nf-core` command line tools offer a number of helpful functions essential for pipeline development with nf-core components.
See [Tools](/tools) for more information.
:::

### Official Nextflow installation.

See [Nextflow installation](https://www.nextflow.io/docs/latest/install.html) for the latest official instructions.

### Quick start installation

Nextflow runs on most POSIX systems (Linux, macOS) and can typically be installed with these commands:

```bash
java -version                           # Check that Java v11+ is installed
curl -s https://get.nextflow.io | bash  # Download Nextflow
chmod +x nextflow                       # Make executable
mv nextflow ~/bin/                      # Add to user's $PATH
```

### Bioconda installation

To install Nextflow with [Bioconda](https://bioconda.github.io/), first set up Bioconda according to the [Bioconda documentation](https://bioconda.github.io/#usage), notably setting up channels:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

:::warning
If channels are not set up this way, dependencies (such as Java) may be installed from
the wrong channels and cause unexpected errors.
:::

It is best practice to use a dedicated Conda environment.
This helps prevent version conflicts and keep your installation clean:

```bash
conda create --name env_nf nextflow
conda activate env_nf
```

To deactivate the Conda environment, run:

```bash
conda deactivate
```

### Windows installation

Windows installation is performed with Windows Subsystem for Linux (WSL), as described on the [Nextflow website](https://nextflow.io/blog/2021/setup-nextflow-on-windows.html).

The main steps include:

- Windows PowerShell installation
- Windows Subsystem for Linux (WSL2) configuration
- Linux distribution installation (on WSL2)

The remaining steps are similar to Linux Nextflow installations.

### Update Nextflow

To update Nextflow, run the self-update command:

```bash
nextflow self-update
```

or run `conda update nextflow`, depending on your original installation method.

### Install specific Nextflow versions

You can install a specific version of Nextflow with the `NXF_VER` environment variable.
Nextflow will automatically install the specified version when the next command is run
(`nextflow self-update` or `nextflow run`).

It is recommended to `export` this bash environment variable. This is helpful for full reproducibility in analysis scripts:

```bash
export NXF_VER=23.10.1
nextflow run hello-world
```

Or prepend the variable to a Nextflow command:

```bash
NXF_VER=23.10.1 nextflow run hello-world
```

### Edge releases

Nextflow has two `stable` releases per year and monthly `edge` releases.
Some pipelines may require the use of Nextflow `edge` releases to use cutting edge features.

Nextflow installs the latest stable version by default.
Specify an edge release either with the exact version in `NXF_VER`,
or with the `NXF_EDGE` environment variable:

```bash
NXF_VER=24.02.0-edge nextflow run hello-world
NXF_EDGE=1 nextflow self-update
```

## Pipeline software

An analysis pipeline chains the execution of multiple tools together.
Historically, all tools would have to be manually installed — often a source of great frustration and a key step where reproducibility between analyses is lost.
nf-core pipelines utilize Nextflow's built-in support for software packaging tools — all can work with Docker and Singularity, and most pipelines also support Conda.

Run your nf-core pipeline with `-profile <type>` (such as `-profile docker` or `-profile conda`) to specify your software packaging tool:

- [Docker](https://docs.docker.com/install/)
  - Typically used locally, on single-user servers, and the cloud
  - Analysis runs in a container, which behaves like an isolated operating system
  - Typically requires system root access, though a "rootless mode" is available
- [Singularity](https://www.sylabs.io/)
  - Often used as an alternative to Docker on HPC systems
  - Analysis runs in a container, optionally created from Docker images
  - Does not require root access or any daemon processes
- [Apptainer](https://apptainer.org/)
  - Open source version of Singularity (split from Singularity in 2021)
  - nf-core pipelines run with `-profile apptainer` build using
    Docker containers instead of pre-built Singularity containers.
  - To use the Singularity containers, use `-profile singularity` instead.
    This works because `apptainer` simply defines `singularity` as an alias
    to the `apptainer` command.

- [Podman](https://podman.io/), [Charliecloud](https://hpc.github.io/charliecloud/), and [Shifter](https://www.nersc.gov/research-and-development/user-defined-images/)
  - All alternatives to Docker, often used on HPC systems
- [Conda](https://conda.io/)
  - Packaging system that manages environments instead of running analyses in containers
  - Poorer reproducibility than Docker or Singularity
    - Low-level changes in package dependencies may occur over time
    - Software still runs in your native operating system environment — core system functions may differ
- [Mamba](https://mamba.readthedocs.io/)
  - A faster implementation of Conda

## Pipeline code

### Automatic

Nextflow automatically fetches pipeline code from GitHub if you specify `nf-core/<pipeline-name>` as the pipeline name.

This method requires an internet connection. See [Running offline](offline.md) to run on systems with no internet access.

### Development

To make changes to a pipeline prior to a run, fork the pipeline GitHub repository and then clone the files. You can then run the pipeline with `nextflow run <path-to-repo>`.

Only do this if you intend to make significant changes to the pipeline. All configuration options can be changed without editing the pipeline code. Forking the pipeline repository means that you cannot use stable releases and will fall behind pipeline updates.
