---
title: Troubleshooting
subtitle: Troubleshoot common issues
shortTitle: Troubleshooting
weight: 11
---

This page covers common issues you may encounter when setting up your nf-core environment.

## Nextflow

### Error: incompatible Java version

An incompatible Java version error can occur when running Nextflow.

This issue occurs when your installed Java version is older than the minimum required version (Java 17). Check your version with `java -version` and install a compatible version using SDKMAN if needed (see [Install Java](./environment_setup/nextflow.md#install-java)).

### Error: `permission denied` when running Nextflow

A `permission denied` error can occur when trying to execute Nextflow after installation.

This issue occurs when the Nextflow executable lacks execution permissions. To resolve this issue, run:

```bash
chmod +x /path/to/nextflow
```

### `nextflow self-update` fails

`nextflow self-update` can fail without applying the update.

This issue occurs when Nextflow is installed in a directory with restricted write permissions. To resolve this issue, either request elevated permissions to update in the current location, or reinstall Nextflow in a directory where you have write permissions.

### Error: `nextflow: command not found`

A `nextflow: command not found` error can occur when running Nextflow from the terminal.

This issue occurs when the Nextflow binary is not in a directory included in your `$PATH`. To resolve this issue, add the directory to your PATH by appending this line to `~/.bashrc` or `~/.zshrc`:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

Then reload your shell configuration with `source ~/.bashrc` (or `source ~/.zshrc`).

### Outdated Nextflow version when installed via Conda

The Nextflow version installed via Conda may be outdated or cause dependency conflicts.

This issue occurs because the Conda-packaged version of Nextflow may lag behind the official release. To resolve this issue, uninstall the Conda version and use the self-installing package instead:

```bash
conda deactivate
conda remove --name nf-env --all
curl -s https://get.nextflow.io | bash
```

### Pipeline fails after updating Nextflow

A pipeline that previously worked can fail after updating Nextflow.

This issue occurs when a major Nextflow release introduces breaking changes. Consult the [Nextflow migration guides](https://nextflow.io/docs/latest/migrations/index.html) and [changelog](https://github.com/nextflow-io/nextflow/releases) for breaking changes that may affect your pipelines.

## nf-core tools

### Error: `nf-core: command not found`

An `nf-core: command not found` error can occur when running nf-core tools from the terminal.

This issue occurs when nf-core tools is not installed in the active environment or is not in your `$PATH`. To resolve this issue:

| Installation method | Solution                                                     |
| ------------------- | ------------------------------------------------------------ |
| Conda               | Run `conda activate nf-core-env` to activate the environment |
| pip                 | Ensure the Python scripts directory is in your `$PATH`       |
| Docker              | Verify your alias is correctly configured in your shell      |

### Installation fails due to Python version

Installation of nf-core tools can fail with an error about your Python version.

This issue occurs because nf-core tools requires Python 3.8 or later. Check your version with `python --version`. To resolve this issue, install a newer Python version using Conda:

```bash
conda create --name nf-core-env python=3.11 nf-core
```

### Error: `permission denied` when installing with pip

A `permission denied` error can occur when installing nf-core tools with pip.

This issue occurs when you do not have write permissions to the system Python directory. To resolve this issue, install to your user directory:

```bash
pip install --user nf-core
```

Alternatively, use a virtual environment:

```bash
python -m venv nf-core-venv
source nf-core-venv/bin/activate
pip install nf-core
```

### Files created by Docker have wrong ownership

Files written by nf-core tools inside a Docker container may be owned by `root` instead of your user.

This issue occurs when the `-u $(id -u):$(id -g)` flag is not passed to Docker, causing the container to run as root. To resolve this issue, always include the flag when running Docker commands, or use the alias provided in the [Docker installation section](./environment_setup/nf-core-tools.md#install-with-docker).

## Dev Containers

### SKU name error when creating a Codespace

A SKU name error can occur when creating a GitHub Codespace for a Dev Container.

This issue occurs when the wrong hardware option is selected. To resolve this issue, select the 4-CPU hardware option.

### Docker daemon connection fails in Dev Containers

Nextflow may fail to connect to the Docker daemon inside a Dev Container.

This issue occurs when using `-profile docker`, which conflicts with the Dev Containers environment. The workaround is to use `-profile singularity` instead.

### Singularity permission errors in Dev Containers

Singularity permission errors can occur when running pipelines inside a Dev Container.

This issue occurs when the Dev Container is missing the required privileges to run Singularity. To resolve this issue, add `"privileged": true` to the `devcontainer.json` file.
