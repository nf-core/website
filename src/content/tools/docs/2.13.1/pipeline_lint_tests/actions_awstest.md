# actions_awstest

#### `PipelineLint.actions_awstest(){:python}`

Checks the GitHub Actions awstest is valid.

In addition to small test datasets run on GitHub Actions, we provide the possibility of testing the pipeline on AWS.
This should ensure that the pipeline runs as expected on AWS (which often has its own unique edge cases).

:::warning
Running tests on AWS incurs costs, so these tests are not triggered automatically.
Instead, they use the `workflow_dispatch` trigger, which allows for manual triggering
of the workflow when testing on AWS is desired.
:::

:::note
You can trigger the tests by going to the Actions tab on the pipeline GitHub repository
and selecting the nf-core AWS test workflow on the left.
:::

The `.github/workflows/awstest.yml` file is tested for the following:

- Must _not_ be turned on for `push` or `pull_request`.
- Must be turned on for `workflow_dispatch`.
