---
title: Installation of nf-core dependencies
subtitle: Install the software requirements needed to run nf-core pipelines.
shortTitle: Dependency installation
weight: 2
---

## Installation requirements

Nextflow must be installed on the system where you launch an nf-core pipeline.

Nextflow can be used on any POSIX-compatible system (Linux, macOS, etc), and on Windows through [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux). It requires Bash 3.2 (or later) and [Java 11 (or later, up to 22)](https://www.oracle.com/java/technologies/downloads/?er=221886).

A personal installation of Nextflow is recommended to simplify updates and version control.

:::note
You don't need to install the `nf-core` command line tools to run nf-core pipelines. However, `nf-core` tools offer a number of helpful commands for users and are essential for pipeline developers. See the [tools page](/tools) for more information.
:::

### Java

You can see which version of Java you have installed using the following command:

```bash
java -version
```

If you don’t have a compatible version of Java installed in your computer, it is recommended that you install it through [SDKMAN!](https://sdkman.io/), and that you use the latest LTS version of Temurin. See [this website](https://whichjdk.com/) for more information.

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

### Nextflow

Nextflow is distributed as a self-installing package, in order to make the installation process as simple as possible:

1. Install Nextflow:

    ```bash
    curl -s https://get.nextflow.io | bash
    ```

This will create the `nextflow` executable in the current directory.

2. Make Nextflow executable:

    ```
    chmod +x nextflow
    ```

3. Move Nextflow into an executable path. For example:

    ```bash
    mkdir -p $HOME/.local/bin/
    mv nextflow $HOME/.local/bin/
    ```

    Ensure the directory `$HOME/.local/bin/` is included in your `PATH` variable. If it is not, add it by setting `export PATH="$PATH:$HOME/.local/bin"`. Setting your `PATH` variable in your `.bashrc` or `.zshrc` file will fix your version of Nextflow across sessions. Alternatively, you could move Nextflow to a directory already in your `PATH`.

    :::warning
    Nextflow will update its executable during the self update process and it should not be placed in a directory with restricted permissions.
    :::

4. Confirm that Nextflow is installed correctly:

    ```bash
    nextflow info
    ```

Alternatively, Nextflow can be installed using [Bioconda](https://bioconda.github.io/):

1. Add conda channels:

    ```bash
    conda config --add channels bioconda
    conda config --add channels conda-forge
    ```

    :::warning
    This is important - if not done, dependencies (such as Java) will be installed from
    the wrong channels and things may break in strange ways.
    :::

2. Create and activate a dedicated Nextflow conda environment.

    ```bash
    conda create --name env_nf nextflow
    conda activate env_nf
    ```

    :::note
    To deactivate the `env_nf` conda environment, run `conda deactivate`.
    :::

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

## Updating Nextflow

With Nextflow installed in your environment, you can update to the latest version using the following command:

```bash
nextflow self-update
```

Specific versions of Nextflow can be installed using the `NXF_VER` environment version.
Nextflow will automatically install that version at the start of the next command.

```bash
export NXF_VER=23.10.1
nextflow run hello-world
```

Setting the `NXF_VER` variable in your `.bashrc` or `.zshrc` file will fix your version of Nextflow across sessions.

You can also temporarily switch to a specific version of Nextflow with the `NXF_VER` environment variable. For example:

```bash
NXF_VER=23.10.0 nextflow run hello
```

:::note
The `conda update nextflow` command can be used to update Nextflow Bioconda installations.
:::

## Edge releases

A stable version of Nextflow is released every six months, in the 4th and 10th month of each year. Additionally, an edge version is released on a monthly basis. The edge releases can be used to access the latest updates and experimental features.

To use the latest edge release, set `NXF_EDGE=1` when updating:

```bash
NXF_EDGE=1 nextflow self-update
```

You can also use `NXF_VER` to temporarily switch to any edge release. For example:

```bash
NXF_VER=24.06.0-edge nextflow info
```

## Software dependencies

Analysis pipelines often chain together the execution of multiple tools.
Historically, all tools would need to be installed manually - often a source of great frustration and a major source of irreproducibility.

nf-core pipelines utilise the built-in support for software packaging that Nextflow offers.
Using profiles, software dependencies can be managed through various packaging (e.g., container runtimes).
To use any of the below, simply execute your nf-core pipeline with the `-profile` option.
For example, `-profile docker` or `-profile singularity`.
The respective tooling for each profile (e.g., [Docker](https://docs.docker.com/install/)) must be installed prior to execution.

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
    :::warning
    Currently, nf-core pipelines run with `-profile apptainer` will build using
    docker containers instead of using pre-built singularity containers.

    To use the singularity containers, use `-profile singularity` instead.
    This works because `apptainer` simply defines `singularity` as an alias
    to the `apptainer` command.
    :::

- [Podman](https://podman.io/), [Charliecloud](https://hpc.github.io/charliecloud/), and [Shifter](https://www.nersc.gov/research-and-development/user-defined-images/)
  - All alternatives to Docker, often used on HPC systems
- [Conda](https://conda.io/)
  - Packaging system that manages environments instead of running analysis in containers
  - Poorer reproducibility than Docker / Singularity
    - There can be changes in low-level package dependencies over time
    - The software still runs in your native operating system environment and so core system functions can differ
- [Mamba](https://mamba.readthedocs.io/)
  - A faster implementation of Conda
