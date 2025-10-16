---
title: Configuration
subtitle: Configuring your Devcontainer environment
weight: 2
type: "tutorial"
---

## Configuration of Codespaces Dev Container

nf-core repositories use [GitHub Codespaces](https://github.com/features/codespaces) with a pre-configured Dev Container for development.

The nf-core Dev Container uses a custom Docker image (`nfcore/devcontainer:latest`) that includes:

- Python 3.13
- nf-core tools from the current workspace
- Nextflow
- nf-test
- Apptainer
- Docker (via the [docker-outside-of-docker feature](https://github.com/devcontainers/features/tree/main/src/docker-outside-of-docker) using the host machine's docker daemon)

VS Code extensions for Python, linting, and nf-core development are pre-installed.

### The Devcontainer Configuration

The main configuration file is located in `.devcontainer/devcontainer.json`, which defines the environment and tools available.

Codespaces and your local VS-Code installation use the configuration file `.devcontainer/devcontainer.json` to setup a docker container.
This configuration uses a base image pre-built from the devcontainer configured in `.devcontainer/build-devcontainer`.
In the pre-built image, `nf-core/tools`, Nextflow, and nf-test is installed as well as various means to install and run further dependencies.
This allows starting the devcontainer quickly.

See this example configuration file, similar to what is used in the tools repository:

```json title=".devcontainer/devcontainer.json"
{
  "name": "nfcore",

  // Base image that installs:
  // nextflow, nf-test, tools, conda, apptainer
  "image": "nfcore/devcontainer:latest",

  // Run options for the container
  "remoteUser": "root",
  "privileged": true,

  // Runs a script once the container is ready
  "onCreateCommand": "./.devcontainer/setup.sh",

  // Plug-and-play devcontainer features
  "features": {
    "ghcr.io/devcontainers/features/docker-outside-of-docker:1.6.3": {}
  }
}
```

### Customizing your Dev Container

You can further customize your environment by editing `.devcontainer/devcontainer.json` to add vs-code extensions, install additional packages or features, as well as configure the environment. Common customizations include:

- Adding VS Code extensions
- Installing additional [devcontainer features](https://containers.dev/implementors/features/)
- Installing additional packages in `setup.sh`
- Adjusting resource requirements

For more details, see the [VS Code Dev Containers documentation](https://code.visualstudio.com/docs/devcontainers/containers).

> **Note:** Codespaces provides a browser-based VS Code experience with all tools pre-installed. You can also use the Dev Container locally with VS Code and Docker.
