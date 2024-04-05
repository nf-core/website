---
title: Test data
subtitle: Follow this walkthrough to add a new pipeline to nf-core.
weight: 20
---

## Running with test data

Whilst the code formatting and nf-core linting tests are helpful, they're not sufficient by themselves.
It's also good actually run your pipeline on a minimal dataset.

It's recommended to use the `debug` profile when running the pipeline, so that you receive warnings about process selectors and other debug info.
A typical command to run an nf-core pipeline test is as follows:

```bash
nextflow run . -profile debug,test,docker --outdir <OUTDIR>
```

We also automatically run tests with GitHub Actions anytime someone updates the pipeline code.
We use nf-test to check the pipeline results, but regular test runs also often catch
syntax errors and other serious problems that cause nextflow to exit with an error.

### Putting the test data on GitHub

To avoid bloating the workflow, we don't keep test data in the same repository as
nf-core workflows.
Instead, we use the dedicated [nf-core/test-datasets](https://github.com/nf-core/test-datasets/) repository.

To set this up, make a fork of that repo to your personal account.
Clone the repository and check out a new branch for your workflow:

```bash
git clone https://github.com/YOUR_USERNAME/test-datasets.git
cd test-datasets
git checkout -b MY_WORKFLOW
```

Now add your test data files - note that they must be **very small**.
GitHub has quite a low file size limit, and the GitHub Actions will time out with anything
that's not tiny. We typically use PhiX / Yeast / part of a chromosome as a reference
and aggressively subsampled input data. I.e., as small as possible, as large as necessary.
We highly recommend that you ask in the [#test-data channel](https://nfcore.slack.com/channels/test-data) on the [nf-core slack](https://nf-co.re/join) for advice
before adding any test data!

Once added, push these new files to GitHub:

```bash
git add .
git commit -m "Added test data for MY_WORKFLOW"
git push --set-upstream origin MY_WORKFLOW
```

Finally, make a pull-request against the main [nf-core/test-datasets](https://github.com/nf-core/test-datasets/) repository with your files.
You want this repo to also use a branch with the name of your workflow, so first go
to the [repository GitHub web page](https://github.com/nf-core/test-datasets/) and create
this new branch using the UI there.
Once created, you can open a pull request and select this as the target branch.

Don't forget to ask for help if you have any doubts!
([Slack](https://nf-co.re/join/slack),alternatively if you already joined Slack you can go directly to the [#help](https://nfcore.slack.com/channels/help) channel)

### Setting up a test workflow

Now that your test data is hosted on the web, you can set up a `test` config profile in your
workflow that points to it.
A stub `test` profile should already exist in `conf/test.config`, so you just need to edit that file.
Switch out the example URLs for the ones you added (view the files on GitHub and click 'Raw' to get the URL).

Add any other required parameters so that running the pipeline runs with as few extra
flags as possible. Note that the `test` profile can be combined with other profiles such as `docker`
or `conda`, so your config should not specify a hardware environment.

Have a go at running the pipeline and see if it works:

```bash
nextflow run MY_WORKFLOW -profile test,docker --outdir <OUTDIR>
```

Note that if you do need to adjust this `nextflow run` command, you'll need to update it
in the `.github/workflows/` YAML files too.

### Test standards summary

- The `test` profile must exist and work
- The test should be as comprehensive as possible
- The test should run as much of the pipeline as possible
- The test dataset should be as small as possible, but big enough to test the pipeline
