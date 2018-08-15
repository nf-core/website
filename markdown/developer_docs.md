
# Introduction
Welcome to the nf-core developer documentation!

Some quick links to help you find your way around:

* Getting Help
    * The quickest place to get help is on the nf-core Gitter channel: [https://gitter.im/nf-core/](https://gitter.im/nf-core/Lobby) - a free chat interface that integrates nicely with GitHub.
    * If you've found a bug in a pipeline, please create a GitHub issue in that pipeline repository.
* [Lint test error codes](errors)
* [Pipeline syncing](sync)

# Guidelines for nf-core pipelines
For now, the guidelines are pretty simple. The below gives an outline of what is required.
For more detail, please see the [**list of lint test error codes**](errors).

**We highly recommend using the [`nf-core create`](tools#create) command** to begin your project.
This tool builds a skeleton pipeline from scratch with all of the features required to pass the below tests.
It also makes it trivial to keep your pipeline in sync with core updates to the nf-core pipeline template (this will be automated, see [pipeline syncing](sync)).

## Features
All pipelines must adhere to the following:

* Be built using Nextflow
* Have an [MIT licence](https://choosealicense.com/licenses/mit/)
* Have software bundled using [Docker](https://www.docker.com/)
* Continuous integration testing
* Stable release tags
* Common pipeline structure and usage
* Excellent documentation
* A responsible contact person / GitHub username
    * This will typically be the main person behind the pipeline development
    * This person should be responsible for basic maintenance and questions
* The pipeline must not have any failures in the `nf-core lint` tests
    * These tests are run by the [nf-core/tools](https://github.com/nf-core/tools) package and validate the requirements listed on this page.
    * You can see the list of tests and how to pass them on the [error codes page](errors).

_(any point can be skipped, given a good enough reason...)_

If possible, it's great if pipelines can also have:

* All software bundled using [bioconda](https://bioconda.github.io/)
    * See below for how to use a bioconda env script automatically build docker and singularity containers, meaning you have a single file to maintain.
* Optimised output file formats
    * Pipelines should generate `CRAM` alignment files by default, but have a `--bam` option to generate `BAM` outputs if required by the user.
* Explicit support for running in cloud environments
    * For example, use of [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) and aws-batch
* Benchmarks from running on cloud environments such as [AWS](https://aws.amazon.com/)

## Pipeline organisation
It's highly recommended that pipelines are built using the [`nf-core create`](tools#create) command, as this will ensure future community updates to be propagated to your pipeline (read more [here](sync)).

## Coding style
The nf-core style requirements are growing and maturing over time.
Typically, as we agree on a new standard we try to build a test for it into the `nf-core lint` command.
As such, to get a feel for what's expected, please read the [lint test error codes](errors).

However, in general, pipelines must:

* Use config profiles to organise anything that is not common for all users
* Run with as little input as possible
    * Metadata (eg. be able to run with just FastQ files, where possible)
    * Reference files (eg. auto-generate missing reference files, where possible)
* Keep only code for the latest stable on the main `master` branch.
    * The main development code should be kept in a branch called `dev`
* Use GitHub releases and keep a detailed changelog file

# How to add a new pipeline

1. Join the [nf-core GitHub organisation](https://github.com/nf-core/nf-co.re/issues/3)
2. Create a pipeline repository in the organisation
    * _If starting from scratch_
        * Ask an admin to create a new pipeline repository for you and add you as a collaborator.
        * This is the best option, as this repository will then be recognised as the "head fork" by GitHub.
    * _If you already have a pipeline_
        * Talk to us about it on the [nf-core gitter chat](https://gitter.im/nf-core/Lobby)
        * Fork your pipeline to the nf-core GitHub organisation
        * Go through [this guide](sync#setup) in order to set up the template syncing
3. Make sure that your pipeline `README.md` file has a big warning on the front saying that it's under development
4. Work on the `dev` branch until the automatic lint tests are passing
5. When you're happy, ping [@nf-core/core](https://github.com/orgs/nf-core/teams/core) for a code review
5. Once the pipeline has been approved, you can remove the development warning on the main readme and the pipeline. Tag a release and the website will be updated.

# Plans for the future
There's always more cool stuff we want to do!

Currently, our focus is on setting up [automatic syncing](sync) for all pipelines so that we can keep core pipeline structure and boilerplate code up to date across all pipelines. This is a huge task, but we made excellent progress at the first nf-core hackathon and hope to get it running soon.

We of course continue to pursue the core aim of nf-core: to assemble a large collection of high quality analysis pipelines covering a wide range of different applications.

Other big and as-yet unstarted pipe dreams include:
* Automated [@nf_core twitter](https://twitter.com/nf_core) tweets when pipelines are released
* A graphical user interface and API to launch pipelines (see [tools#61](https://github.com/nf-core/tools/issues/61))
* A web-based monitoring tool to listen to nextflow pipeline event broadcasts and visualise running workflows

Just about every nf-core repository has a healthy list of GitHub issues listing planned improvements, so if you're hungry for more, go and have an explore!
