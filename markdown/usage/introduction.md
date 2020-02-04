---
title: Getting started
subtitle: How to run your first nf-core pipeline.
---

## What is nf-core?
nf-core is a community effort to collect a curated set of analysis pipelines built using [Nextflow](https://www.nextflow.io/docs/latest/index.html).

nf-core has three target audiences: facilities, single users and developers.
For facilities it provides highly automated and optimized pipelines that guaranty reproducibility of results for their users.
Single users profit from portable, documented and easy to use workflows.
But you can also become a developer and write your own pipeline in Nextflow using already available templates and helper tools.

## What is Nextflow?
Nextflow is a *workflow manager*.
It has been developed specifically to ease the creation and execution of bioinformatics pipelines.
The benefits of having your pipeline in Nextflow include:

* Built-in GitHub support.
* Compatibility with virtually all computational infrastructures, including all major cluster job schedulers.
* Integrated software dependency management (Docker, Singularity, Conda).
* Portability so you can run your pipeline anywhere: laptop, cluster or cloud.
* Reproducibility of analyses independent of time and computing platform.

Whether your pipeline is a simple BLAST execution or a complex genome annotation pipeline, you can build it with Nextflow.

## How to run a pipeline

Nextflow works best when you have an active internet connection, as it is able to fetch all pipeline requirements. If you need to run offline, please see [_Running offline_](offline.md).

1. First, make sure that you have all required software installed (Nextflow + Docker / Singularity / Conda). See the [installation docs](installation.md) for more information.

2. Choose a pipeline to run. See the available pipelines at [https://nf-co.re/pipelines](https://nf-co.re/pipelines). If you have [nf-core/tools](tools.md) installed, run `nf-core list`

3. Configure Nextflow to run on your system.

    * The simplest way to run is with `-profile docker` (or `singularity` /  `conda`) which will tell Nextflow to execute jobs locally using Docker to fulfil the software requirements.

    * For more complex configuration of Nextflow for your system, please see the [_Nextflow configuration_](https://nf-co.re/usage/configuration) documentation.

4. To test that everything is working properly, try running the tests for your pipeline of interest in the terminal:

    ```bash
    nextflow run nf-core/<pipeline_name> -profile test,docker
    ```

    * Replace `<pipeline_name>` with the name of an nf-core pipeline.

    * There is no need to download anything first - nextflow will pull the code for you from the GitHub repository automatically and fetch the software requirements too.

    * If you don't have Docker installed, replace `docker` in the command with either `singularity` or `conda`.

    * If the pipeline fails, check the [troubleshooting docs](/usage/troubleshooting.md) and ask for help on the nf-core Slack channel for that particular pipeline (see [https://nf-co.re/join](https://nf-co.re/join)).

5. Read the pipeline documentation to see which command-line parameters are required. These will be specific to your data type and usage.

6. Launch the pipeline with some real data by omitting the `test` config profile and providing the required pipeline-specific parameters. For example, if you want to run the `methylseq` pipeline, you might use the following command:

    ```bash
    nextflow run nf-core/methylseq -profile docker -reads 'input_data/*.fastq.gz' --outdir myproj/results --genome GRCh38
    ```

7. Once complete, check the pipeline execution and quality control reports. Each pipeline comes with documentation describing the different outputs.

## Tips and tricks

* Hyphens matter! Core Nextflow command-line options use one (`-`) whereas pipeline specific parameters use two (`--`)
* Specify `--email your@email.com` to receive emails when your pipeline completes
* Always specify `-r <version-number>` when running to explicitly use a specific release. Then an identical command can be used in the future to give identical results.
* Use `-resume` to restart pipelines that did not complete. This ensures that successful tasks from the previous run wont be re-executed.
* Use `nextflow log` to find names of all previous runs in your directory. These can be used with `-resume` to restart specific runs.
* Be clever with multiple Nextflow configuration locations.. For example, use `-profile` for your cluster configuration, `~/.nextflow/config` for your personal config such as `params.email` and a working directory `nextflow.config` file for reproducible run-specific configuration.

## Helper tools

To help you manage your nf-core pipelines and discover updates, we have written some command-line helper tools.
These allow you to list all available pipelines and versions, with information about what versions you're running locally.
There are also commands to help downloading pipelines for use offline.

To find out more about these tools, read the [Tools](/tools) page.
