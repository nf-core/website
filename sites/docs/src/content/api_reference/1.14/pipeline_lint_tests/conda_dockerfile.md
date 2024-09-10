# conda_dockerfile

#### `PipelineLint.conda_dockerfile(){:python}`

Checks the Dockerfile for use with Conda environments

:::note
This test only runs if there is both an `environment.yml`
and `Dockerfile` present in the pipeline root directory.
:::

If a workflow has a conda `environment.yml` file, the `Dockerfile` should use this
to create the docker image. These files are typically very short, just creating the conda
environment inside the container.

This linting test checks for the following:

- All of the following lines are present in the file (where `PIPELINE` is your pipeline name):
  > ```Dockerfile
  > FROM nfcore/base:VERSION
  > COPY environment.yml /
  > RUN conda env create --quiet -f /environment.yml && conda clean -a
  > RUN conda env export --name PIPELINE > PIPELINE.yml
  > ENV PATH /opt/conda/envs/PIPELINE/bin:$PATH
  > ```
- That the `FROM nfcore/base:VERSION` is tagged to the most recent release of nf-core/tools
  > - The linting tool compares the tag against the currently installed version of tools.
  > - This line is not checked if running a development version of nf-core/tools.

:::note
Additional lines and different metadata can be added to the `Dockerfile`
without causing this lint test to fail.
:::
