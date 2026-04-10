---
title: nf-core tools
subtitle: Install nf-core/tools
shortTitle: nf-core tools
weight: 4
---

nf-core tools is a Python package that provides command-line utilities for working with nf-core pipelines.
While optional, it offers helpful features for downloading, launching, and developing pipelines.

The nf-core tools package provides commands for:

- Listing available nf-core pipelines and components
- Downloading pipelines for offline use
- Launching pipelines with customized parameters
- Creating and developing new pipelines and components
- Linting and validating pipeline code

## Installing nf-core tools

nf-core tools can be installed using Conda, pip, or Docker.
Choose the method that best fits your environment.

### Install with Conda

Conda is the recommended installation method as it handles all dependencies automatically.

To install nf-core tools in a dedicated environment with nf-core tools and Nextflow:

1. Create a new environment:

   ```bash
   conda create --name nf-core-env nf-core nextflow
   ```

2. Activate the environment:

   ```bash
   conda activate nf-core-env
   ```

3. Verify the installation:

   ```bash
   nf-core --version
   nextflow -version
   ```

### Install with pip

To install nf-core tools with pip:

1. Install the package:

   ```bash
   pip install nf-core
   ```

1. Verify the installation:

   ```bash
   nf-core --version
   ```

:::note
When using pip, ensure you have Python 3.8 or later installed.
You may need to use `pip3` instead of `pip` depending on your system configuration.
:::

### Install with Docker

To use nf-core tools with Docker:

1. Pull the `nfcore/tools` Docker image:

   ```bash
   docker pull nfcore/tools
   ```

1. Run nf-core tools commands using Docker:

   ```bash
   docker run -itv `pwd`:`pwd` -w `pwd` -u $(id -u):$(id -g) nfcore/tools --help
   ```

1. (Optional) Create an alias to simplify commands:
   1. Add an alias in your `~/.bashrc` or `~/.zshrc`

      ```bash
      alias nf-core="docker run -itv \`pwd\`:\`pwd\` -w \`pwd\` -u $(id -u):$(id -g) nfcore/tools"
      ```

   1. Run `nf-core` directly:

      ```bash
      nf-core --help
      ```

## Updating nf-core tools

To keep nf-core tools up to date with the latest features and bug fixes, update regularly using your installation method.

### Update Conda install

To update your nf-core tools Conda install:

1. Activate the environment:

   ```bash
   conda activate nf-core-env
   ```

2. Run the update command:

```bash
conda update nf-core
```

### Update pip install

To update your nf-core tools pip install, run:

```bash
pip install --upgrade nf-core
```

### Update Docker install

To update nf-core tools docker image, run:

```bash
docker pull nfcore/tools
```
