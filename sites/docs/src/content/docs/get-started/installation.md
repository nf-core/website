---
title: Installation
subtitle: Install the software requirements needed to run nf-core pipelines.
shortTitle: Installation
weight: 2
---

## Installation

Nextflow must be installed on the system where you launch an nf-core pipeline.

Nextflow is executable on any POSIX-compatible system (Linux, macOS, etc.), and on Windows through [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux). It requires Bash 3.2 (or later) and [Java 11 (or later, up to 22)](https://www.oracle.com/java/technologies/downloads/?er=221886).

:::note
Use a personal installation of Nextflow where possible, instead of a system-wide installation.
:::

### Install Java

If you don’t have a compatible version of Java installed, it is recommended that you install the latest LTS version of Temurin through [SDKMAN!](https://sdkman.io/) See [Which Version of JDK Should I Use?](https://whichjdk.com/) for more information.

To install Java with SDKMAN:

1. Install SDKMAN:

```bash
curl -s https://get.sdkman.io | bash
```

2. Open a new terminal.

3. Install Java:

   ```bash
   sdk install java 17.0.10-tem
   ```

4. Confirm that Java is installed correctly:

   ```bash
   java -version
   ```

### Install Nextflow

Nextflow is distributed as a self-installing package that can be installed in a few simple steps.

1. Install Nextflow:

   ```bash
   curl -s https://get.nextflow.io | bash
   ```

2. Make Nextflow executable:

   ```
   chmod +x nextflow
   ```

3. Move `nextflow` to a location on your systems `$PATH`. For example:

   ```bash
   mkdir -p $HOME/.local/bin/
   mv nextflow $HOME/.local/bin/
   ```

   Ensure the directory `$HOME/.local/bin/` is included in your `PATH` variable. You can set export `PATH="$PATH:$HOME/.local/bin"` to temporarily add a directory to your `PATH`. You can also add this command to your shell configuration file (i.e. `~/.bashrc`, `~/.bash_profile`, or `~/.zshrc`) to add it to your `PATH` at every login.

   :::warning
   Nextflow will update its executable during the self-update process and it should not be placed in a directory with restricted permissions.
   :::

4. Confirm that Nextflow is installed correctly:

   ```bash
   nextflow info
   ```

Alternatively, you can also install Nextflow with [Bioconda](https://bioconda.github.io/).

To install Nextflow with Bioconda:

1. Add conda channels:

   ```bash
   conda config --add channels bioconda
   conda config --add channels conda-forge
   ```

   :::warning
   Conda channels are the locations where packages are distributed from. If the channels are not configured accurately your dependencies may be installed from
   the wrong channels and be incompatible.
   :::

2. Create and activate a dedicated Nextflow conda environment.

   ```bash
   conda create --name env_nf nextflow
   conda activate env_nf
   ```

   :::note
   Use `conda deactivate env_nf` deactivate the `env_nf` environment.
   :::

### Update Nextflow

A stable version of Nextflow is released in the 4th and 10th month of each year.

Run the self-update command to update to the latest version of Nextflow:

```bash
nextflow self-update
```

An edge version of Nextflow is released monthly and can be used to access the latest updates and experimental features.

Set NXF_EDGE=1 when updating Nextflow to use the latest edge release:

```bash
NXF_EDGE=1 nextflow self-update
```

The `NXF_VER` environmental variable can be set to use a specific version of Nextflow.

Add `NXF_VER` to your execution command to temporarily use a different version of Nextflow:

```bash
NXF_VER=24.04.4 nextflow info
```

:::note
Set the `NXF_VER` variable in your `.bashrc` or `.zshrc` if you want to fix your version of Nextflow for all sessions.
:::

Alternatively, set NXF_VER to pin the version of Nextflow:

```bash
export NXF_VER=24.04.4
```

### Install nf-core tools

You don't need to install nf-core tools to run nf-core pipelines. However, nf-core tools offer a number of helpful commands for users and are essential for pipeline developers.

For users:

- List nf-core pipelines
- Download nf-core pipelines
- Configure pipeline runs
- Launch pipeline runs

For developers:

- Create nf-core components from nf-core templates
- Lint nf-core components
- Manage nf-core pipeline schemas
- And more...

See [nf-core tools](/docs/nf-core-tools) for installation instructions.

## Software dependencies

Analysis pipelines often chain together the execution of multiple tools.
Historically, all tools would need to be installed manually and was often a source of great frustration and irreproducibility.

nf-core pipelines utilise Nextflow's built-in support for software packaging that Nextflow offers.
Using profiles, software dependencies can be managed through various packaging (e.g., Docker, Singularity, and Conda).
These ensure that you have the correct tool and version of the tool that each pipeline needs to correctly execute.

The respective tooling for a profile must be installed prior to pipeline execution. Follow the links below to install the required profile tooling:

- Typically requires system root access, though a _"rootless mode"_ is available
- Needs system administrator for install
- [Docker](https://docs.docker.com/install/)
  - Typically used locally, on single-user servers, and the cloud
  - Analysis runs in a container, which behaves like an isolated operating system
  - Typically requires system root access, though a _"rootless mode"_ is available
  - Requires system administrator for install
- [Singularity](https://www.sylabs.io/)
  - Often used as an alternative to Docker on HPC systems
  - Also runs _containers_, and can optionally create these from Docker images
  - Does not require root access or any daemon processes
  - Requires system administrator for install
- [Apptainer](https://apptainer.org/)

  - Open source version of Singularity (split from Singularity in 2021)

:::warning
Currently, nf-core pipelines run with `-profile apptainer` will build using
docker containers instead of the pre-built singularity containers. To use the singularity containers, use `-profile singularity` instead.
This works because `apptainer` simply defines `singularity` as an alias
to the `apptainer` command.
:::

- [Podman](https://podman.io/), [Charliecloud](https://hpc.github.io/charliecloud/), and [Shifter](https://www.nersc.gov/research-and-development/user-defined-images/)
  - All alternatives to Docker, often used on HPC systems
- [Conda](https://conda.io/)
  - Packaging system that manages environments instead of running analysis in containers
  - Can be installed by any user
  - Can be installed by any user
  - Less reproducible than Docker / Singularity
    - Low-level changes in package dependencies may occur over time
    - Software still runs in your native operating system environment — core system functions may differ
- [Mamba](https://mamba.readthedocs.io/)
  - A faster implementation of Conda

## Windows installation

The installation procedure for Windows computers is more complex.

The main steps include:

- Installing Windows PowerShell
- Configuring the Windows Subsystem for Linux (WSL2)
- Installing a Linux distribution (on WSL2)

See the [guide for setting up a Nextflow environment on Windows 10](https://nextflow.io/blog/2021/setup-nextflow-on-windows.html) for more information.

:::warning
Some information in the [guide for setting up a Nextflow environment on Windows 10](https://nextflow.io/blog/2021/setup-nextflow-on-windows.html) may be out of date.
:::
