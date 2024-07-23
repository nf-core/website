---
title: Test data
subtitle: Follow this walkthrough to add a new pipeline to nf-core.
short_title: Test data setup
weight: 20
type: 'tutorial'
---

## Running with test data

Whilst the code formatting and nf-core linting tests are helpful, they're not sufficient by themselves.
It's also important actually run the pipeline on a minimal dataset.

It's recommended to use the `debug` profile when running the pipeline, so that you receive warnings about process selectors and other debug information.

A typical command to run an nf-core pipeline test is as follows:

```bash
nextflow run . -profile debug,test,docker --outdir <OUTDIR>
```

Tests are also automatically run using GitHub Actions anytime pipeline code is updated.
The [nf-test](https://www.nf-test.com/) plugin is used to check pipeline outputs. However, regular testing using the test profile will also catch
syntax errors and other serious problems that cause Nextflow to exit with an error before they are pushed to GitHub.

### Putting the test data on GitHub

To avoid bloating pipeline repositories, test data is kept and shared in the [nf-core/test-datasets](https://github.com/nf-core/test-datasets/) repository.

To create a new test dataset, make a fork of that repo to your personal account.
Clone the repository and check out a new branch for your pipeline:

```bash
git clone https://github.com/YOUR_USERNAME/test-datasets.git
cd test-datasets
git checkout -b MY_PIPELINE
```

:::note
Test data files must be **very small**. GitHub has a low file size limit, and the GitHub Actions will time out with anything
that's too big. Common test date includes PhiX, Yeast, or part of a chromosome, and are aggressively subsampled.
As a rule, test data should be as small as possible and as large as necessary.
We highly recommend that you ask in the [#test-data channel](https://nfcore.slack.com/channels/test-data) on the [nf-core slack](https://nf-co.re/join) for advice before adding any test data!
:::

Push your new files to GitHub once they are added:

```bash
git add .
git commit -m "Added test data for MY_PIPELINE"
git push --set-upstream origin MY_PIPELINE
```

Make a pull-request against the main [nf-core/test-datasets](https://github.com/nf-core/test-datasets/) repository.
You want this repo to use a branch with the name of your pipeline, so first go
to the [repository GitHub web page](https://github.com/nf-core/test-datasets/) and use the UI to create
a new branch with the name of your pipeline.
Once created, open a pull request and select this as the target branch.

:::note
Use ([Slack](https://nf-co.re/join/slack) to ssk for [#help](https://nfcore.slack.com/channels/help) channel) if you have any doubts!
:::

### Setting up a test profile

A stub `test` profile already exist in `conf/test.config` and can be modified to point to your new test data.
Replace the example URLs for the new files you have just added by viewing the files on GitHub. Click _Raw_ to get the URL for each file.

Add any other required parameters to the test profile so that the pipeline runs with as few flags as possible. Note that the `test` profile can be combined with other profiles such as `docker`
or `conda`, so your config should not specify any hardware requirements.

Run the pipeline to see if it works:

```bash
nextflow run MY_PIPELINE -profile test,docker --outdir <OUTDIR>
```

Note that if you do need to adjust the `nextflow run` command, you'll need to update it
in the `.github/workflows/` YAML files for GitHub Actions to work.
