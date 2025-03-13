# actions_awsfulltest

#### `PipelineLint.actions_awsfulltest(){:python}`

Checks the GitHub Actions awsfulltest is valid.

In addition to small test datasets run on GitHub Actions, we provide the possibility of testing the pipeline on full size datasets on AWS.
This should ensure that the pipeline runs as expected on AWS and provide a resource estimation.

The GitHub Actions workflow is called `awsfulltest.yml`, and it can be found in the `.github/workflows/` directory.

:::warning
This workflow incurs AWS costs, therefore it should only be triggered for pipeline releases:
`workflow_run` (after the docker hub release workflow) and `workflow_dispatch`.
:::

:::note
You can manually trigger the AWS tests by going to the Actions tab on the pipeline GitHub repository and selecting the
nf-core AWS full size tests workflow on the left.
:::

:::note
For tests on full data prior to release, [Nextflow Tower](https://tower.nf) launch feature can be employed.
:::

The `.github/workflows/awsfulltest.yml` file is tested for the following:

- Must be turned on `workflow_dispatch`.
- Must be turned on for `workflow_run` with `workflows: ["nf-core Docker push (release)"]` and `types: [completed]`
- Should run the profile `test_full` that should be edited to provide the links to full-size datasets. If it runs the profile `test`, a warning is given.
