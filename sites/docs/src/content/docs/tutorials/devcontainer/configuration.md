---
title: Configuration
subtitle: Configuring your Devcontainer environment
weight: 2
type: "tutorial"
---

## Configuration of Devcontainers

nf-core repositories use pre-configured Devcontainers that can be run in [GitHub Codespaces](https://github.com/features/codespaces) or locally with VS Code.
These environments ensure all necessary software is available for developing and testing Nextflow pipelines and modules, while maintaining fast startup times.

The Devcontainer utilizes a pre-built custom Docker image (`nfcore/devcontainer:latest`).
This image includes the core tools for nf-core development and various means to access further dependencies:

- Python
- nf-core tools
- Nextflow
- nf-test
- Apptainer
- Docker (via the [docker-outside-of-docker feature](https://github.com/devcontainers/features/tree/main/src/docker-outside-of-docker) using the host machine's docker daemon)

VS Code extensions for Python, linting, and nf-core development are pre-installed.

### The Devcontainer Configuration

The file `.devcontainer/devcontainer.json` is the main configuration manifest, which defines the environment, the base image, features, and more.

This example configuration file contains some key settings, similar to those used in nf-core repositories:

```json title=".devcontainer/devcontainer.json"
{
  "name": "nfcore",

  // Base image that installs:
  // nextflow, nf-test, tools, conda, apptainer
  "image": "nfcore/devcontainer:latest",

  // Run options for the container
  "remoteUser": "root",

  // Essential for running Apptainer (used with '-profile singularity')
  // Allows to perform the necessary mount operations within the container
  "privileged": true,

  // Executes the ./.devcontainer/setup.sh script after
  // the environment is created to perform final configuration
  "onCreateCommand": "./.devcontainer/setup.sh",

  // Install modular "plug-and-play" devcontainer features
  // Specifically the docker-outside-of-docker feature
  // (https://containers.dev/implementors/features/)
  "features": {
    "ghcr.io/devcontainers/features/docker-outside-of-docker:1.6.3": {}
  }
}
```

For a comprehensive list of settings see the [relevant docs](https://containers.dev/implementors/json_reference/#general-properties) of the devcontainer specification

### Customizing your Dev Container

You can further customize your environment by editing `.devcontainer/devcontainer.json` to add vs-code extensions, install additional tools, or adjust the environment.
Common customizations include:

- Adding VS Code extensions for specific languages or linters
- Installing additional packages by modifying the `setup.sh`
- Integrate additional [devcontainer features](https://containers.dev/implementors/features/)
- Adjusting resource requirements for local Devcontainers and Codespaces
- Mounting paths from the host machien when running locally

For more details, see the [VS Code Dev Containers documentation](https://code.visualstudio.com/docs/devcontainers/containers).
