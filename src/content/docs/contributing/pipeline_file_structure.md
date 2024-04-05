---
title: Pipeline file structure
subtitle: Description of the elements of an nf-core template
---

When you create a new pipeline from the nf-core template, you are presented with
a directory full of many files and folders. Some are boilerplate code that you'll
probably never need to touch. Some are important to understand whilst developing a pipeline.
Here we list the outputs of a pipeline from the template with all of the options ticked,
describing what each file is and what it does.

## Files

You will find the following files in each nf-core pipeline. They are automatically generated, when running `nf-core create`.

- `main.nf`: This is the main nextflow file which will get executed if the pipeline is run. Typically, parameters are initialized and validated in this script before a workflow from the `workflows/` directory is called for execution.

* `nextflow.config`: The main nextflow configuration file. It contains the default pipeline parameters, nextflow configuration options and information like pipeline and minimum nextflow version, among others.
  The `nextflow.config` also defines different configuration profiles that can be used to run the pipeline. See the [Configuration docs](/docs/usage/configuration) for more information.

- `README.md`: Basic information about the pipeline and usage

- `nextflow_schema.json`: The JSON schema file is used for pipeline parameter specification. This is automatically created using the `nf-core schema build` command. It is used for printing command-line help, validating input parameters, building the website docs and for building pipeline launch interfaces (web and cli).

- `CHANGELOG.md`: Information about the changes made to the pipeline for each release.

- `LICENSE`: The license - should be MIT

- `CODE_OF_CONDUCT.md`: The nf-core code of conduct.

- `CITATIONS.md`: All citations needed when using the pipeline

- `.gitattributes`: Git settings, primarily getting the `.config` files to render with Nextflow syntax highlighting on <github.com>

- `.gitignore`: Files that should be ignored by git.

- `.editorconfig`: Editorconfig file that helps assuring consistent coding style

- `.prettierrc.yml`: Prettier lint configuration file to assure consistent markdown files

- `.prettierignore`: Files that should be ignored by prettier

- `modules.json`: This file holds information (e.g. version) about all the modules in the pipeline that have been installed from `nf-core/modules`

- `.nf-core.yml`: Indicates the type of repository (pipeline or modules repo)

- `.gitpod.yml`: Config file to define online working environment with <https://www.gitpod.io>

- `pyproject.toml`: Config file for Python. Mostly used to configure linting of `bin/check_samplesheet.py` with Black

## Directories

- `.devcontainer`: Configuration to work with the [GitHub Codespaces](https://github.com/features/codespaces) online editing environments.

- `.github/`: Other GitHub specific files, e.g. for specifying templates and GitHub actions

- `assets/`: Any additional files needed for the pipeline

- `bin/`: Directory for scripts that must be directly accessible within a pipeline process. Anything in this directory can be directly called from within Nextflow processes.

- `conf/`: Configuration files, including a `base.config` file which is always loaded into `nextflow.config` and describes basic pipeline configurations, like CPU and memory usage for processes with low, medium and high requirements. Additionally, most pipelines also have a `igenomes.config` file which describes the locations of popular genomes that can be automatically downloaded for a pipeline run. Finally, the `test.config` and `test_full.config` files are test configurations that are loaded during test runs. Since DSL2, it also contains a `modules.config` file, which defines module-specific configurations and is explained further down in the "DSL2 and modules" section.

- `docs/`: Markdown files for documenting the pipeline

- `modules/`: Contains pipeline-specific and common nf-core modules

- `workflows/`: Contains the main pipeline workflows to be executed in the `main.nf` file

- `subworkflows/`: Contains smaller subworkflows that typically consist out of a few modules chained together
