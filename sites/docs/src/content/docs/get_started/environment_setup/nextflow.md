---
title: Installing Nextflow
subtitle: Install Nextflow
shortTitle: Nextflow
weight: 2
---

Nextflow is the workflow management system that runs all nf-core pipelines.
This page will walk you through installing and configuring Nextflow on your system.

:::note{title="Prerequisites"}
You will need the following to get started:

- Bash version 3.2 or later
- Java version 17 or later (up to Java 25)

For the latest instructions, see the [Nextflow installation guide](https://www.nextflow.io/docs/latest/install.html).
:::

## Install Java

The easiest way to install Java is using SDKMAN:

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

## Install Nextflow

Nextflow is distributed as an easy to use self-installing package.
It is also distributed via Conda and as a standalone distribution.

The self-installing package is recommended as it ensures you always get the latest version and handles dependencies automatically.

### Self-installing package

To install Nextflow with the self-installing package:

1. Download and install Nextflow:

   ```bash
   curl -s https://get.nextflow.io | bash
   ```

2. Make the binary executable:

   ```bash
   chmod +x nextflow
   ```

3. Move Nextflow to a directory in your `$PATH`:

   ```bash
   mkdir -p $HOME/.local/bin/
   mv nextflow $HOME/.local/bin/
   ```

   :::note
   Common directories on Linux or OSX operating systems for executables include `$HOME/.local/bin/`, `/usr/local/bin/`, or `$HOME/bin/`.
   Make sure the directory you choose is in your `$PATH` environment variable.
   :::

4. Verify the installation:

   ```bash
   nextflow info
   ```

### Conda installation

To install Nextflow with Conda:

1. Create an environment with Nextflow:

   ```bash
   conda create --name nf-env bioconda::nextflow
   ```

2. Activate the environment:

   ```bash
   source activate nf-env
   ```

3. Confirm Nextflow is installed correctly:

   ```bash
   nextflow info
   ```

:::warning
The Conda installation method is not recommended as it may install outdated versions and cause dependency conflicts.
Use the self-installing package method when possible.
:::

### Standalone distribution

For offline or restricted environments, you can download a pre-built executable from the [Nextflow GitHub releases](https://github.com/nextflow-io/nextflow/releases).

To use the standalone distribution:

1. Download the standalone distribution from Assets section of the GitHub releases page.

2. Grant execution permissions to the downloaded file. For example:

   ```bash
   chmod +x nextflow-25.10.0-dist
   ```

3. Move Nextflow to a directory in your `$PATH`:

   ```bash
   mkdir -p $HOME/.local/bin/
   mv nextflow-25.10.0-dist $HOME/.local/bin/
   ```

4. Use it as a drop-in replacement for nextflow command. For example:

   ```bash
   nextflow-25.10.0-dist info
   ```

Replace `25.10.0` with your desired version number.

## Update Nextflow

To use the latest Nextflow features or run specific versions of pipelines, you may need to update or pin your version of Nextflow.

### Update to the latest stable version

To update Nextflow to the latest stable release, run:

```bash
nextflow self-update
```

### Pin a specific version

To pin a specific Nextflow version:

- For a single command, prefix the Nextflow command with the `NXF_VER` variable:

  ```bash
  NXF_VER=23.10.0 nextflow run nf-core/rnaseq
  ```

- For your current terminal session, export the `NXF_VER` variable:

  ```bash
  export NXF_VER=23.10.0
  ```

- For persistent usage, add the export `NXF_VER` variable command to `~/.bashrc` or `~/.zshrc`:

  ```bash
  export NXF_VER=23.10.0
  ```

:::tip
Pinning a specific Nextflow version can help ensure reproducibility across different runs and environments.
:::
