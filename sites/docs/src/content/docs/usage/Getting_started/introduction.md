---
title: Getting started
subtitle: What is nf-core and how to run a pipeline
weight: 1
parentWeight: 10
---

## What is nf-core?

nf-core is a community effort to collect a curated set of analysis pipelines built with [Nextflow](https://www.nextflow.io/docs/latest/index.html).

nf-core has three target audiences: facilities, single users, and developers.
For facilities, it provides highly automated and optimised pipelines that guarantee reproducibility of results for users.
Single users benefit from portable, documented, and easy-to-use workflows.
Developers can write Nextflow pipelines with nf-core's ready-made templates and tools.

## What is Nextflow?

Nextflow is a _workflow manager_.
It has been developed specifically to ease the creation and execution of bioinformatics pipelines.
The benefits of creating your pipeline with Nextflow include:

- Built-in GitHub support.
- Compatibility with virtually all computational infrastructures, including all major cluster job schedulers.
- Integrated software dependency management (Docker, Singularity, Conda).
- Portability to run your pipelines anywhere: laptop, cluster, or cloud.
- Reproducibility of analyses, independent of time and computing platform.

Whether your pipeline is a simple BLAST execution or a complex genome annotation pipeline, you can build it with Nextflow.

## How to run a pipeline

Nextflow works best with an active internet connection, as it is able to fetch all pipeline requirements. See [Running offline](offline.md) if you need to run Nextflow pipelines without an internet connection.

To run a pipeline:

1. All software dependencies must be installed (Nextflow + Docker / Singularity / Conda). See [Dependency installation](installation.md) for more information.
   - Run Nextflow's "hello world" example to confirm all tools are installed and operational:

     ```bash
     nextflow run hello
     ```

2. Choose a pipeline to run. See the available [nf-core pipelines](https://nf-co.re/pipelines). If you have [nf-core tools](https://nf-co.re/tools) installed, run `nf-core pipelines list`.

3. Configure Nextflow to run on your system:
   - The simplest way to run is with `-profile docker` (or `singularity`). This instructs Nextflow to execute jobs locally, with Docker (or Singularity) to fulfill software dependencies.

   - Conda is also supported with `-profile conda`. However, this option is not recommended, as reproducibility of results can't be guaranteed without containerization.

   - If you are a member of one of the [listed institutions](https://github.com/nf-core/configs#documentation), use the [institutional config file](https://github.com/nf-core/configs/tree/master/conf) created for your institution.

   - See [Nextflow configuration](https://nf-co.re/docs/usage/getting_started/configuration) for advanced Nextflow configuration options.

4. Run the tests for your pipeline in the terminal to confirm everything is working:

   ```bash
   nextflow run nf-core/<pipeline_name> -profile test,docker --outdir <OUTDIR>
   ```

   - Replace `<pipeline_name>` with the name of an nf-core pipeline.

   - If you don't have Docker installed, replace `docker` with either `singularity` or `conda`.

   - Nextflow will pull the code from the GitHub repository and fetch the software requirements automatically, so there is no need to download anything first.

   - If the pipeline fails, see [Troubleshooting](/docs/usage/troubleshooting/overview) or ask for help on the [nf-core Slack channel](https://nf-co.re/join) for your pipeline.

5. Read the pipeline documentation to see which command-line parameters are required. This will be specific to your data type and usage.

6. To launch the pipeline with real data, omit the `test` config profile and provide the required pipeline-specific parameters. For example, to run the `methylseq` pipeline, your command will be similar to this:

   ```bash
   nextflow run nf-core/methylseq -profile docker --input 'input_data/*.fastq.gz' --outdir myproj/results --genome GRCh38
   ```

7. Once complete, check the pipeline execution and quality control reports (such as `multiqc_report.html` files for [MultiQC](https://multiqc.info/docs/usage/pipelines/#nextflow) reports). Each pipeline's documentation describes the outputs to expect.

## Tips and tricks

- Hyphens matter! Core Nextflow command-line options use one (`-`), whereas pipeline-specific parameters use two (`--`).
- Specify `--email your@email.com` to receive emails when your pipeline completes (requires Nextflow [mail and notification configuration](https://www.nextflow.io/docs/latest/mail.html#mail-configuration)).
- To receive notifications from your pipeline (e.g. on completion or failure), see the [nf-slack plugin](https://github.com/seqeralabs/nf-slack) for Slack or the [nf-teams plugin](https://github.com/nvnieuwk/nf-teams) for Microsoft Teams.
- Include `-r <version-number>` to specify a release version explicitly. This guarantees the same run command will give the same results in future runs.
- Use `-resume` to restart pipelines that did not complete. This uses cached results for successful tasks from the previous run, instead of executing all tasks from scratch.
- Use `nextflow log` to find the names of all previous runs in your directory, then use `nextflow run <pipeline> -resume <run-name>` to restart a specific run.
- Use multiple Nextflow configuration locations to your benefit. For example, use `-profile` for your cluster configuration, `~/.nextflow/config` for your personal configuration (with `params.email`, for example), and a working directory `nextflow.config` file for reproducible run-specific configuration.
- If you use Singularity, we recommend that you specify a cache directory with the [Nextflow environment variable](https://www.nextflow.io/docs/latest/config.html#environment-variables) `NXF_SINGULARITY_CACHEDIR` in your `~./bash_profile` or `~/.bashrc` during the installation. This will store all your container images in one place, instead of downloading an image each time you run a pipeline. Specify only the base directory — Nextflow handles the folders and file names for you.

## Helper tools

nf-core includes command-line tools to help you manage your nf-core pipelines and discover updates.
These allow you to list all available pipelines and versions, with information about the versions you're running locally.
You can also download pipelines for offline use.

See [nf-core/tools](/docs/nf-core-tools) for more information.
