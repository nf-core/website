# actions_ci

#### PipelineLint.actions_ci()

Checks that the GitHub Actions pipeline CI (Continuous Integration) workflow is valid.

The `.github/workflows/ci.yml` GitHub Actions workflow runs the pipeline on a minimal test
dataset using `-profile test` to check that no breaking changes have been introduced.
Final result files are not checked, just that the pipeline exists successfully.

This lint test checks this GitHub Actions workflow file for the following:

- Workflow must be triggered on the following events:
  ```yaml
  on:
    push:
    branches:
      - dev
    pull_request:
    release:
    types: [published]
  ```
- The minimum Nextflow version specified in the pipeline’s `nextflow.config` matches that defined by `nxf_ver` in the test matrix:

  ```yaml
  strategy:
    matrix:
      # Nextflow versions: check pipeline minimum and current latest
      nxf_ver: ['19.10.0', '']
  ```

  #### NOTE

  These `matrix` variables run the test workflow twice, varying the `nxf_ver` variable each time.
  This is used in the `nextflow run` commands to test the pipeline with both the latest available version
  of the pipeline (`''`) and the stated minimum required version.

- The Docker container for the pipeline must use the correct pipeline version number:

  > - Development pipelines:
  >   ```bash
  >   docker tag nfcore/<pipeline_name>:dev nfcore/<pipeline_name>:dev
  >   ```
  > - Released pipelines:
  >   ```bash
  >   docker tag nfcore/<pipeline_name>:dev nfcore/<pipeline_name>:<pipeline-version>
  >   ```
  > - Complete example for a released pipeline called _nf-core/example_ with version number `1.0.0`:
  >   ```yaml
  >   - name: Build new docker image
  >     if: env.GIT_DIFF
  >     run: docker build --no-cache . -t nfcore/example:1.0.0
  >   ```

  > - name: Pull docker image
  >   if: ${{ !env.GIT_DIFF }}
  >   run: |
  >   docker pull nfcore/example:dev
  >   docker tag nfcore/example:dev nfcore/example:1.0.0
  >
  > ```
  >
  > ```
