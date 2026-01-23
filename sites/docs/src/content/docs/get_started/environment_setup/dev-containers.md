---
title: Dev Containers
subtitle: Learn how to use and create Dev Containers
shortTitle: Dev Containers
weight: 6
---

nf-core repositories offer [Development Container](https://containers.dev/) (dev container) configurations.
Dev containers offer a portable development environment that can be executed in a web browser, such as cloud-based [GitHub Codespaces](https://github.com/codespaces) or local Visual Studio Code.
These environments package all required software in isolated, containerised spaces for development.

:::note
Devcontainers are not required to run or develop nf-core pipelines.
However, some may find dev containers useful for these purposes.
:::

## Set up in GitHub Codespaces

[GitHub Codespaces](https://github.com/codespaces) provides a browser-based development platform that resembles local VS Code.
The free tier allows up to 120 hours monthly for basic compute instances.
Pipeline repositories use 4-core machines with 16GB RAM and 32GB storage.

To launch a Codespace:

1. Navigate to any nf-core pipeline repository on GitHub
2. Click the green **Code** button
3. Select the **Codespaces** tab
4. Click **Create codespace**

Select the 4-CPU hardware option for adequate performance.

No additional setup is required beyond a GitHub account.
Git comes preinstalled in the environment and Codespaces automatically configures authentication via GitHub.

## Set up in Visual Studio Code

To run Dev Containers locally in VS Code:

1. Install Docker Desktop
2. Install VS Code
3. Install the VS Code Remote Development Extension
4. Clone the repository in a container volume

Local VS Code Dev Containers require additional SSH key forwarding setup for remote repository access.

## Testing your environment

Once your environment is running, verify the setup by opening a terminal in VS Code or Codespaces and running:

```bash
nf-core --version
```

And try [running a test pipeline](../get_started/run-your-first-pipeline.md).

## Further Configuration

The `.devcontainer/devcontainer.json` file serves as the main configuration manifest, that describes the container environment to install and how it should be configured.

### Pre-built image contents

The Dev Containers uses `nfcore/devcontainer:latest` and includes:

- Python
- nf-core tools
- Nextflow
- nf-test
- Apptainer
- Docker (via docker-outside-of-docker feature)
- Pre-installed VS Code extensions for Python, linting, and nf-core development

See [nf-core extension pack](./vs-code.md#nf-core-extension-pack) for more information about VS Code packages.

## Troubleshooting

### SKU name errors when creating Codespaces

**Problem:**

SKU name errors when creating a Codespace.

**Cause:**

Incorrect hardware option selected.

**Solution:**

Select the 4-CPU hardware option.

### Docker daemon connection failures

**Problem:**

Connection failures to the Docker daemon.

**Cause:**

Using `-profile docker` in the Dev Containers environment.

**Solution:**

Use `-profile singularity` instead.

### Singularity permission errors

**Problem:**

Singularity permission errors.

**Cause:**

Missing privileges in the Dev Containers configuration.

**Solution:**

Add `"privileged": true` to the `devcontainer.json` file.
