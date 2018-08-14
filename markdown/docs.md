
# Introduction
These pages are still very much under construction and are likely to change a lot in the near future. If you have thoughts please join the discussion!

## Getting Help {#getting-help}
The quickest place to get help is on the nf-core Gitter channel: [https://gitter.im/nf-core/](https://gitter.im/nf-core/Lobby) - a free chat interface that integrates nicely with GitHub.

# Guidelines for New Pipelines
For now, the guidelines are pretty simple. The below gives an outline of what is required.
For more detail, please see the [**list of lint test error codes**](errors).

However, we recommend using **nf-core tool's** [`create`](tools.md#create) subcommand, that builds a pipeline from scratch. It generates a skeleton pipeline with all of the features
required to pass the below tests.

## Features
All pipelines must adhere to the following:

* Be built using Nextflow
* An [MIT licence](https://choosealicense.com/licenses/mit/)
* Software bundled using [Docker](https://www.docker.com/)
    * This must be at least one `Dockerfile` in the repository
    * Automatic builds will be set up with tagged versions for GitHub releases
    * Pipelines should ship with Nextflow profiles for singularity that pull from the docker repository.
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

* Software bundled using [bioconda](https://bioconda.github.io/)
    * See below for how to use a bioconda env script automatically build docker and singularity containers, meaning you have a single file to maintain.
* Optimised output file formats
    * Pipelines should generate `CRAM` alignment files by default, but have a `--bam` option to generate `BAM` outputs if required by the user.
* Explicit support for running in cloud environments
    * For example, use of [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/)
* Benchmarks from running on cloud environments such as [AWS](https://aws.amazon.com/)

## Pipeline organisation
It's highly recommended that pipelines are built using the [cookiecutter](https://github.com/nf-core/cookiecutter) starter template, as future developments are likely to be based on this assumption (see [_future plans_](#plans-for-the-future) below).

## Coding style

Pipelines must:

* Use config profiles to organise anything that is not common for all users
* Run with as little input as possible
    * Metadata (eg. be able to run with just FastQ files, where possible)
    * Reference files (eg. auto-generate missing reference files, where possible)
* Keep only code for the latest stable on the main `master` branch.
    * The main development code should be kept in a branch called `development`
* Use GitHub releases and keep a detailed changelog file
* _..and many more :)_

# How to add a new pipeline

1. Join the [nf-core GitHub organisation](https://github.com/nf-core/nf-co.re/issues/3)
2. Create a pipeline repository in the organisation
    * _If starting from scratch_
        * Ask an admin to create a new pipeline repository for you and add you as a collaborator.
        * This is the best option, as this repository will then be recognised as the "head fork" by GitHub.
    * _If you already have a pipeline_
        * Just fork your pipeline to the nf-core GitHub organisation
3. Make sure that your pipeline `README.md` file has a big warning on the front saying that it's under development
4. When you're happy with it, ping [@nf-core/admin](https://github.com/orgs/nf-core/teams/admin) for a code review
5. Once the pipeline has been approved, you can remove the development warning on the homepage and the pipeline will be added to the website.

# Plans for the future
## Base change automation
In the future it would be great to have automated bots listen to changes in a base pipeline cookiecutter template. Changes (eg. in pull-requests) could be compared against all other pipelines; any that also apply elsewhere could then be flagged or even modified with automated pull-requests. The reverse could also apply - changes in a downstream pipeline that are shared in the core template could be acted upon.

## Helper scripts
It's also possible that we could write a python package with helper scripts to make it easier to create custom config files, list available pipelines, check for updates, pull singularity images and so on. This could then be packaged using bioconda and the python package index as its own stand-alone tool.

Helper scripts could also be very useful when developing and testing pipelines. For example, linting required features and code style.
