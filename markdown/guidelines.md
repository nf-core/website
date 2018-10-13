# Requirements for nf-core pipelines
If you're thinking of adding a new pipeline to nf-core, please read the documentation
about [adding a new pipeline](/adding_pipelines).

## Workflow size
We aim to have a _"not too big, not too small"_ rule with nf-core pipelines.
This is a little fuzzy, but as a rule of thumb workflows should contain at
least three different processes and be simple enough to run that a new user
can realistically run the pipeline after spending ten minutes reading the docs.

## Minimum requirements
All nf-core pipelines _must_ adhere to the following:

* Be built using Nextflow
* Have an [MIT licence](https://choosealicense.com/licenses/mit/)
* Have software bundled using [Docker](https://www.docker.com/)
* Continuous integration testing
* Stable release tags
* Common pipeline structure and usage
* Run in a single command (not multiple sub-workflows)
    * It is ok to have workflows that use the output of _another_ nf-core pipeline as input
* Excellent documentation and GitHub repository keywords
* A responsible contact person / GitHub username
    * This will typically be the main person behind the pipeline development
    * This person should be responsible for basic maintenance and questions
* The pipeline must not have any failures in the `nf-core lint` tests
    * These tests are run by the [nf-core/tools](https://github.com/nf-core/tools) package and validate the requirements listed on this page.
    * You can see the list of tests and how to pass them on the [error codes page](errors).

## Recommended features
If possible, it's great if pipelines can also have:

* All software bundled using [bioconda](https://bioconda.github.io/)
    * Nearly all nf-core pipelines use a conda `env` script to list their software requirements.
    The pipeline Docker images are then built using this, meaning that with a single file your pipeline can support nextflow users running with conda, docker or singularity.
    * The [nf-core template](/tools#creating-a-new-workflow) comes with all required code to support this setup.
* Optimised output file formats
    * Pipelines should generate `CRAM` alignment files by default, but have a `--bam` option to generate `BAM` outputs if required by the user.
* Explicit support for running in cloud environments
    * For example, use of [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) and aws-batch
    * The [nf-core template](/tools#creating-a-new-workflow) comes with all required code to support this setup.
* Benchmarks from running on cloud environments such as [AWS](https://aws.amazon.com/)

## Coding style
The nf-core style requirements are growing and maturing over time.
Typically, as we agree on a new standard we try to build a test for it into the `nf-core lint` command.
As such, to get a feel for what's expected, please read the [lint test error codes](errors).

However, in general, pipelines must:

* Use [config profiles](https://www.nextflow.io/docs/latest/config.html) to organise hardware-specific options
* Run with as little input as possible
    * Metadata (eg. be able to run with just FastQ files, where possible)
    * Reference files (eg. auto-generate missing reference files, where possible)
* Keep only code for the latest stable on the main `master` branch.
    * The main development code should be kept in a branch called `dev`
* Use GitHub releases and keep a detailed changelog file
