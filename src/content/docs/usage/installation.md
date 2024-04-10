---
title: Installation
subtitle: Install the software requirements needed to run nf-core pipelines.
weight: 2
---

## Nextflow

All nf-core pipelines use Nextflow, so this must be installed on the system where you launch your analysis.
Once installed, see [_Nextflow configuration_](configuration.md) to set up Nextflow to run on your system.

We recommend using a personal installation of Nextflow where possible, instead of using a system-wide installation.
This makes it easier to update and control versions.

:::tip{collapse title="Official Nextflow docs"}
If in doubt, please ee the [Nextflow installation docs](https://www.nextflow.io/docs/latest/getstarted.html#installation)
for the latest instructions.
Installation details are included here for convenience and may not be up to date.
:::

:::note
You don't need to install the `nf-core` command line tools to run nf-core pipelines, you only need Nextflow.
However, they offer a number of helpful functions and are essential for pipeline development using nf-core components.
See the [tools page](/tools) for more information.
:::

### Typical installation

Nextflow runs on most POSIX systems (Linux, macOS, etc) and can typically be installed by running these commands:

```bash
java -version                           # Check that Java v11+ is installed
curl -s https://get.nextflow.io | bash  # Download Nextflow
chmod +x nextflow                       # Make executable
mv nextflow ~/bin/                      # Add to user's $PATH
```

### Bioconda installation

You can also install Nextflow using [Bioconda](https://bioconda.github.io/).
First, set up Bioconda according to the [Bioconda documentation](https://bioconda.github.io/#usage), notably setting up channels:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

:::warning
This is important - if not done, dependencies (such as Java) will be installed from
the wrong channels and things may break in strange ways.
:::

A best practice with conda is to use a dedicated conda environment.
This can help to prevent version conflicts and keep everything clean:

```bash
conda create --name env_nf nextflow
conda activate env_nf
```

To deactivate the conda environment, run:
```bash
conda deactivate
```

### Installation on Windows

For Windows the installation procedure is more complex and is fully described on the [Nextflow website](https://nextflow.io/blog/2021/setup-nextflow-on-windows.html).

The main steps will be the following:

- Install Windows PowerShell
- Configure the Windows Subsystem for Linux (WSL2)
- Install a Linux distribution (on WSL2)

The step to install Nextflow in itself will afterwards be the same as previously mentioned.

### Updating Nextflow

Updating nextflow is as simple as running 
```bash
nextflow self-update
```
or `conda update nextflow`, depending on how it was installed.

### Specific Nextflow versions

You can install a specific version of Nextflow by using the `NXF_VER` environment version.
Nextflow will automatically install that version at the start of the next command
(be it `nextflow self-update` or `nextflow run`).

You can `export` this bash env variable, which is good to do in analysis scripts for full reproducibility:

```bash
export NXF_VER=23.10.1
nextflow run hello-world
```

Or even just prepend it to a Nextflow command:

```bash
NXF_VER=23.10.1 nextflow run hello-world
```

### Edge releases

Nextflow has two `stable` releases per year and monthly `edge` releases.
Some pipelines may require the use of Nextflow `edge` releases in order to exploit _cutting edge_ features.

Nextflow installs the latest stable version by default.
You can get an edge release either by defining the exact version with `NXF_VER`
or by using the `NXF_EDGE` environment variable:

```bash
NXF_VER=24.02.0-edge nextflow run hello-world
NXF_EDGE=1 nextflow self-update
```

## Pipeline software

An analysis pipeline chains the execution of multiple tools together.
Historically, all tools would have to be manually installed — often a source of great frustration and a key step where reproducibility between analyses is lost.
nf-core pipelines utilise the built-in support for software packaging that Nextflow offers: all can work with Docker and Singularity, and most pipelines also support Conda.

To use any of the below, simply run your nf-core pipeline with `-profile <type>`.
For example, `-profile docker` or `-profile conda`.

- [Docker](https://docs.docker.com/install/)
  - Typically used locally, on single-user servers, and the cloud
  - Analysis runs in a _container_, which behaves like an isolated operating system
  - Typically requires system root access, though a _"rootless mode"_ is available
- [Singularity](https://www.sylabs.io/)
  - Often used as an alternative to Docker on HPC systems
  - Also runs _containers_, and can optionally create these from Docker images
  - Does not need root access or any daemon processes
- [Apptainer](https://apptainer.org/)

  - Open source version of Singularity (split from Singularity in 2021)
  - :::warning
    Currently, nf-core pipelines run with `-profile apptainer` will build using
    docker containers instead of using pre-built singularity containers.

    To use the singularity containers, use `-profile singularity` instead.
    This works because `apptainer` simply defines `singularity` as an alias
    to the `apptainer` command.
    :::

- [Podman](https://podman.io/), [Charliecloud](https://hpc.github.io/charliecloud/) and [Shifter](https://www.nersc.gov/research-and-development/user-defined-images/)
  - All alternatives to Docker, often used on HPC systems
- [Conda](https://conda.io/)
  - Packaging system that manages environments instead of running analysis in containers
  - Poorer reproducibility than Docker / Singularity
    - There can be changes in low-level package dependencies over time
    - The software still runs in your native operating system environment and so core system functions can differ
- [Mamba](https://mamba.readthedocs.io/)
  - A faster reimplementation of Conda

## Pipeline code

### Automatic

The pipeline needs no installation - Nextflow will automatically fetch it from GitHub if `nf-core/<pipeline-name>` is specified as the pipeline name.

This method requires an internet connection. If you're running on a system that has no internet connection, please see [Running Offline](offline.md).

### Development

If you would like to make changes to the pipeline, fork the GitHub repository and then clone the files. Once cloned, you can run the pipeline with `nextflow run <path-to-repo>`.

Note that you should _only_ do this if you intend to make significant changes to the pipeline. All configuration options can be changed without editing the pipeline code. Forking the pipeline repositories means that you cannot use stable releases and you will fall behind new updates.

## Reference genomes

Some pipelines come with built-in support for genome reference files.
We recommend downloading these references locally to avoid fetching the same reference many times.
For more information, see [_Reference genomes_](reference_genomes.md).
