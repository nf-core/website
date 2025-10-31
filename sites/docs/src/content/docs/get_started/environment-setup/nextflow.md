---
title: Installing Nextflow
subtitle: Learn how to install Nextflow
shortTitle: Nextflow
weight: 1
---

Nextflow is the workflow management system that runs all nf-core pipelines. This page will walk you through installing and configuring Nextflow on your system.

:::note{title="Prerequisites"}
You will need the following to get started:

- Bash version 3.2 or later
- Java version 17 or later (up to Java 25)
:::

## Install Nextflow

Nextflow is distributed as an easy to use self-installing package. It is also distributed via Conda and as a standalone distribution.

The self-installing package is recommended as it ensures you always get the latest version and handles dependencies automatically.

### Self-installing package

To install Nextflow with the self-installing package:

1. Download and install Nextflow:

    ```bash
    curl -s https://get.nextflow.io | bash
    ```

1. Make the binary executable:

    ```bash
    chmod +x nextflow
    ```

1. Move Nextflow to a directory in your `$PATH`:

    ```bash
    mkdir -p $HOME/.local/bin/
    mv nextflow $HOME/.local/bin/
    ```

    :::note
    Common directories for executables include `$HOME/.local/bin/`, `/usr/local/bin/`, or `$HOME/bin/`. Make sure the directory you choose is in your `$PATH` environment variable.
    :::

1. Verify the installation:

    ```bash
    nextflow info
    ```

### Conda installation

To install Nextflow with Conda:

1. Create an environment with Nextflow:

    ```bash
    conda create --name nf-env bioconda::nextflow
    ```

1. Activate the environment:

    ```bash
    source activate nf-env
    ```

1. Confirm Nextflow is installed correctly:

    ```bash
    nextflow infoactivate nf-env
    ```

:::warning
The Conda installation method is not recommended as it may install outdated versions and cause dependency conflicts. Use the self-installing package method when possible.
:::

### Standalone distribution

For offline or restricted environments, you can download a pre-built executable from the [Nextflow GitHub releases](https://github.com/nextflow-io/nextflow/releases).

To use the standalone distribution:

1. Download the standalone distribution from Assets section of the GitHub releases page.

1. Grant execution permissions to the downloaded file. For example:

    ```bash
    chmod +x nextflow-25.10.0-dist
    ```

1. Use it as a drop-in replacement for nextflow command. For example:

    ```bash
    ./nextflow-25.10.0-dist run info
    ```


Replace `25.10.0` with your desired version number.

## Install Java

The easiest way to install it is using SDKMAN:

To install Java with SDKMAN:

1. Install SDKMAN:

    ```bash
    curl -s https://get.sdkman.io | bash
    ```

1. Open a new terminal.

1. Install Java:

    ```bash
    sdk install java 17.0.10-tem
    ```

1. Confirm that Java is installed correctly:

    ```bash
    java -version
    ```

## Update Nextflow

To use the latest Nextflow features or run specific versions of pipelines, you may need to update or pin your version of Nextflow.

### Update to the latest stable version

To update Nextflow to the latest stable release, run:

```bash
nextflow self-update
```

### Update to the latest edge version

To update to the latest edge release, run:

```bash
NXF_EDGE=1 nextflow self-update
```

To update to the latest stable release, run:

```bash
NXF_EDGE=0 nextflow self-update
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

## Troubleshooting

### Java version issues

**Problem:** Error about incompatible Java version

**Solution:** Ensure you have Java 17 or later installed. Check your version with `java -version` and install a compatible version using SDKMAN (see [Install Java](#installing-java)).

### Permission denied errors

**Problem:** Cannot execute Nextflow after installation

**Solution:** Make sure the Nextflow executable has execution permissions:

```bash
chmod +x /path/to/nextflow
```

### Update failures

**Problem:** `nextflow self-update` fails

**Solution:** The update can fail if Nextflow is installed in a directory with restricted permissions. Either:

- Request elevated permissions to update in the current location
- Reinstall Nextflow in a directory where you have write permissions

### Command not found

**Problem:** `nextflow: command not found`

**Solution:** Ensure Nextflow is in a directory included in your `$PATH`. Add the directory to your PATH by adding this line to `~/.bashrc` or `~/.zshrc`:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

Then reload your shell configuration:

```bash
source ~/.bashrc  # or source ~/.zshrc
```

### Conda installation issues

**Problem:** Outdated Nextflow version or dependency conflicts with Conda

**Solution:** Uninstall the Conda version and use the self-installing package method instead:

```bash
conda deactivate
conda remove --name nf-env --all
curl -s https://get.nextflow.io | bash
```

### Checking for migration requirements

**Problem:** Pipeline fails after updating Nextflow

**Solution:** When updating across major stable releases, consult the [Nextflow migration guides](https://nextflow.io/docs/latest/migration.html) and [changelog](https://github.com/nextflow-io/nextflow/releases) for breaking changes that may affect your pipelines.

## Additional resources

For more information about installing Nextflow, see:

- [Nextflow installation](https://www.nextflow.io/docs/latest/install.html)
