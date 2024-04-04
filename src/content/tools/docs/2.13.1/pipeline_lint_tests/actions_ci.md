# actions\_ci

#### `PipelineLint.actions_ci(){:python}`

Checks that the GitHub Actions pipeline CI (Continuous Integration) workflow is valid.

The `.github/workflows/ci.yml` GitHub Actions workflow runs the pipeline on a minimal test
dataset using `-profile test` to check that no breaking changes have been introduced.
Final result files are not checked, just that the pipeline exists successfully.

This lint test checks this GitHub Actions workflow file for the following:

* Workflow must be triggered on the following events:
  ```yaml
  on:
      push:
      branches:
          - dev
      pull_request:
      release:
      types: [published]
  ```
* The minimum Nextflow version specified in the pipeline’s `nextflow.config` matches that defined by `NXF_VER` in the test matrix:

  ```yaml
  strategy:
    matrix:
      # Nextflow versions: check pipeline minimum and current latest
      NXF_VER: ['19.10.0', '']
  ```

  :::note
  These `matrix` variables run the test workflow twice, varying the `NXF_VER` variable each time.
  This is used in the `nextflow run` commands to test the pipeline with both the latest available version
  of the pipeline (`''`) and the stated minimum required version.
  :::
