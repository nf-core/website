---
title: Troubleshooting
subtitle: Common issues and solutions for getting started with nf-core
shortTitle: Troubleshooting
weight: 10
---

This page covers common issues you may encounter when setting up your nf-core environment.

## Nextflow

### Java version issues

If you encounter errors about an incompatible Java version, ensure you have Java 17 or later installed. Check your version with `java -version` and install a compatible version using SDKMAN if needed (see [Install Java](./environment_setup/nextflow.md#install-java)).

### Permission denied errors

If you cannot execute Nextflow after installation, the executable likely lacks execution permissions. Make sure the Nextflow executable has execution permissions by running `chmod +x /path/to/nextflow`.

### Update failures

If `nextflow self-update` fails, the update may be blocked because Nextflow is installed in a directory with restricted permissions. You can either request elevated permissions to update in the current location, or reinstall Nextflow in a directory where you have write permissions.

### Command not found

If you see `nextflow: command not found`, Nextflow is not in a directory included in your `$PATH`. Add the directory to your PATH by adding this line to `~/.bashrc` or `~/.zshrc`:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

Then reload your shell configuration with `source ~/.bashrc` (or `source ~/.zshrc`).

### Conda installation issues

If you experience an outdated Nextflow version or dependency conflicts with Conda, uninstall the Conda version and use the self-installing package method instead:

```bash
conda deactivate
conda remove --name nf-env --all
curl -s https://get.nextflow.io | bash
```

### Pipeline fails after updating

If your pipeline fails after updating Nextflow, you may need to check for breaking changes. When updating across major stable releases, consult the [Nextflow migration guides](https://nextflow.io/docs/latest/migrations/index.html) and [changelog](https://github.com/nextflow-io/nextflow/releases) for breaking changes that may affect your pipelines.

## nf-core tools

### Command not found

If you see `nf-core: command not found`, ensure nf-core tools is installed in an active environment. For Conda installations, run `conda activate nf-core-env` if you installed in a separate environment. For pip installations, ensure the Python scripts directory is in your `$PATH`. For Docker, check that your alias is properly configured in your shell.

### Python version issues

If installation fails due to your Python version, note that nf-core tools requires Python 3.8 or later. Check your version with `python --version`. If necessary, install a newer Python version using Conda with `conda create --name nf-core-env python=3.11 nf-core`.

### Permission errors with pip

If you encounter permission denied errors when installing with pip, install in your user directory instead using `pip install --user nf-core`. Alternatively, use a virtual environment:

```bash
python -m venv nf-core-venv
source nf-core-venv/bin/activate
pip install nf-core
```

### Docker file permissions

If files created by Docker have the wrong ownership, always include the `-u $(id -u):$(id -g)` flag when running Docker commands, or use the alias provided in the [Docker installation section](./environment_setup/nf-core-tools.md#install-with-docker).

## Dev Containers

### SKU name errors when creating Codespaces

If you encounter SKU name errors when creating a Codespace, you likely selected the incorrect hardware option. Select the 4-CPU hardware option instead.

### Docker daemon connection failures

If you experience connection failures to the Docker daemon in the Dev Containers environment, this is caused by using `-profile docker`. Use `-profile singularity` instead.

### Singularity permission errors

If you see Singularity permission errors, the Dev Containers configuration is missing required privileges. Add `"privileged": true` to the `devcontainer.json` file.
