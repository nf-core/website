---
title: "Installation"
subtitle: Installation
weight: 1
---

nf-core tools is a Python package that provides command-line utilities for working with nf-core pipelines. While optional, it offers helpful features for downloading, launching, and developing pipelines.

The nf-core tools package provides commands for:

- Listing available nf-core pipelines and components
- Downloading pipelines for offline use
- Launching pipelines with customized parameters
- Creating and developing new pipelines and components
- Linting and validating pipeline code

## Installing nf-core tools

nf-core tools can be installed using Conda, pip, or Docker. Choose the method that best fits your environment.

### Install with Conda

Conda is the recommended installation method as it handles all dependencies automatically.

To install nf-core tools with Conda:

1. Install nf-core tools in your current environment:

   ```bash
   conda install nf-core
   ```

1. Verify the installation:

   ```bash
   nf-core --version
   ```

Alternatively, create a dedicated environment with both nf-core tools and Nextflow:

1. Create a new environment:

   ```bash
   conda create --name nf-core-env nf-core nextflow
   ```

1. Activate the environment:

   ```bash
   conda activate nf-core-env
   ```

1. Verify the installation:

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
When using pip, ensure you have Python 3.8 or later installed. You may need to use `pip3` instead of `pip` depending on your system configuration.
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

   - `-i` and `-t` are needed for the interactive CLI prompts to work (this tells Docker to use a pseudo-tty with stdin attached)
   - The `-v` argument tells Docker to bind your current working directory (`pwd`) to the same path inside the container, so that files created there will be saved to your local file system outside of the container
   - `-w` sets the working directory in the container to this path, so that it's the same as your working directory outside of the container
   - `-u` sets your local user account as the user inside the container, so that any files created have the correct ownership permissions

1. (Optional) Create an alias to simplify commands:

   1. Add an alias in your `~/.bashrc` or `~/.zshrc`:

      ```bash
      alias nf-core="docker run -itv \`pwd\`:\`pwd\` -w \`pwd\` -u $(id -u):$(id -g) nfcore/tools"
      ```

   1. Run `nf-core` directly:

      ```bash
      nf-core --help
      ```

#### Docker version tags

You can use Docker image tags to specify the version you would like to use. For example, `nfcore/tools:dev` for the latest development version of the code, or `nfcore/tools:1.14` for version `1.14` of tools. If you omit this, it will default to `:latest`, which should be the latest stable release.

If you need a specific version of Nextflow inside the container, you can build an image yourself. Clone the repo locally and check out whatever version of nf-core/tools that you need. Then build using the `--build-arg NXF_VER` flag:

```bash
docker build -t nfcore/tools:dev . --build-arg NXF_VER=20.04.0
```

## Updating nf-core tools

To keep nf-core tools up to date with the latest features and bug fixes, update regularly using your installation method.

### Update Conda install

To update your nf-core tools Conda install, run:

```bash
conda update nf-core
```

### Update pip install

To update your nf-core tools pip install, run:

```bash
pip install --upgrade nf-core
```

### Update Docker install

To update nf-core tools Docker image, run:

```bash
docker pull nfcore/tools
```

## Development version

If you would like to use the latest development version of tools, install directly from the GitHub repository:

```bash
pip install --upgrade --force-reinstall git+https://github.com/nf-core/tools.git@dev
```

If you intend to make edits to the code, first make a fork of the repository and then clone it locally. Go to the cloned directory and install with pip (this also installs development requirements):

```bash
pip install --upgrade -r requirements-dev.txt -e .
```

## Advanced usage

### Using a specific Python interpreter

You can run nf-core tools with a specific Python interpreter. The command line usage and flags are exactly the same as if you ran with the `nf-core` command. Note that the module is `nf_core` with an underscore, not a hyphen like the console command.

For example:

```bash
python -m nf_core --help
python3 -m nf_core list
~/my_env/bin/python -m nf_core create --name mypipeline --description "This is a new skeleton pipeline"
```

### Using with your own Python scripts

The nf-core tools functionality can be imported into your own Python scripts. For example, to get a list of all available nf-core pipelines:

```python
import nf_core.list
wfs = nf_core.list.Workflows()
wfs.get_remote_workflows()
for wf in wfs.remote_workflows:
    print(wf.full_name)
```

For detailed function documentation, see the [nf-core tools API reference](https://nf-co.re/tools/docs/).

### Automatic version check

nf-core tools automatically checks for new versions when run. If you would prefer to skip this check, set the `NFCORE_NO_VERSION_CHECK` environment variable:

```bash
export NFCORE_NO_VERSION_CHECK=1
```

### Shell completions

Auto-completion for the `nf-core` command is available for bash, zsh, and fish. To activate it, add the following lines to the respective shell config files:

| Shell | Config file                               | Command                                             |
| ----- | ----------------------------------------- | --------------------------------------------------- |
| bash  | `~/.bashrc`                               | `eval "$(_NF_CORE_COMPLETE=bash_source nf-core)"`  |
| zsh   | `~/.zshrc`                                | `eval "$(_NF_CORE_COMPLETE=zsh_source nf-core)"`   |
| fish  | `~/.config/fish/completions/nf-core.fish` | `eval (env _NF_CORE_COMPLETE=fish_source nf-core)` |

After restarting your shell session, you should have auto-completion for the `nf-core` command and all its sub-commands and options.

:::note
The added line will run the `nf-core` command on shell startup (which may slow down startup time). You should either have nf-core tools installed globally, or wrap the command in a conditional:

- For bash and zsh: `if type nf-core > /dev/null; then eval "$(_NF_CORE_COMPLETE=bash_source nf-core)"; fi`
- For fish: `if command -v nf-core &> /dev/null eval (env _NF_CORE_COMPLETE=fish_source nf-core) end`
  :::

:::tip
If you see the error `command not found compdef`, ensure your config file contains `autoload -Uz compinit && compinit` before the eval line.
:::
