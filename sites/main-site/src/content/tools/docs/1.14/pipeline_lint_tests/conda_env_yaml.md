# conda_env_yaml

#### `PipelineLint.conda_env_yaml(){:python}`

Checks that the conda environment file is valid.

:::note
This test is ignored if there is not an `environment.yml`
file present in the pipeline root directory.
:::

DSL1 nf-core pipelines use a single Conda environment to manage all software
dependencies for a workflow. This can be used directly with `-profile conda`
and is also used in the `Dockerfile` to build a docker image.

This test checks the conda `environment.yml` file to ensure that it follows nf-core guidelines.
Each dependency is checked using the [Anaconda API service](https://api.anaconda.org/docs).
Dependency sublists are ignored with the exception of `- pip`: these packages are also checked
for pinned version numbers and checked using the [PyPI JSON API](https://wiki.python.org/moin/PyPIJSON).

Specifically, this lint test makes sure that:

- The environment `name` must match the pipeline name and version
  > - The pipeline name is defined in the config variable `manifest.name`
  > - Replace the slash with a hyphen as environment names shouldn’t contain that character
  > - Example: For `nf-core/test` version 1.4, the conda environment name should be `nf-core-test-1.4`
- All package dependencies have a specific version number pinned
  > :::warning
  > Remember that Conda package versions should be pinned with one equals sign (`toolname=1.1`),
  > but pip uses two (`toolname==1.2`)
  > :::
- That package versions can be found and are the latest available
  > - Test will go through all conda channels listed in the file, or check PyPI if `pip`
  > - Conda dependencies with pinned channels (eg. `conda-forge::openjdk`) are ok too
  > - In addition to the package name, the pinned version is checked
  > - If a newer version is available, a warning will be reported
