---
title: Getting started
subtitle: How to run your first nf-core pipeline.
menu:
  main:
    weight: 10
---

## What is nf-core?

nf-core is a community effort to collect a curated set of analysis pipelines built using [Nextflow](https://www.nextflow.io/docs/latest/index.html).

nf-core has three target audiences: facilities, single users, and developers.
For facilities, it provides highly automated and optimised pipelines that guarantee reproducibility of results for their users.
Single users profit from portable, documented, and easy-to-use workflows.
However, you can also become a developer and write your own pipelines in Nextflow using ready-made templates and helper tools.

## What is Nextflow?

Nextflow is a _workflow manager_.
It has been developed specifically to ease the creation and execution of bioinformatics pipelines.
The benefits of having your pipeline in Nextflow include:

- Built-in GitHub support.
- Compatibility with virtually all computational infrastructures, including all major cluster job schedulers.
- Integrated software dependency management (Docker, Singularity, Conda).
- Portability to run your pipelines anywhere: laptop, cluster, or cloud.
- Reproducibility of analyses independent of time and computing platform.

Whether your pipeline is a simple BLAST execution or a complex genome annotation pipeline, you can build it with Nextflow.

## How to run a pipeline

Nextflow works best when you have an active internet connection, as it is able to fetch all pipeline requirements. However, if you need to run offline, see [_Running offline_](offline.md).

1. First, make sure that you have all software dependencies installed (Nextflow + Docker / Singularity / Conda). See the [installation docs](installation.md) for more information.

   - Try running the Nextflow "hello world" example to make sure that the tools are working properly:

     ```bash
     nextflow run hello
     ```

2. Choose a pipeline to run. See the available pipelines at [https://nf-co.re/pipelines](https://nf-co.re/pipelines). If you have [https://nf-core/tools](https://nf-co.re/tools) installed, run `nf-core list`.

3. Configure Nextflow to run on your system.

   - The simplest way to run is with `-profile docker` (or `singularity`) which instructs Nextflow to execute jobs locally using Docker to fulfill the software dependencies.

   - Conda is also supported with `-profile conda`. However, this option is not recommended, as reproducibility of results can't be guaranteed without containerization.

   - If you are a member of one of the institutions listed in the [documentation](https://github.com/nf-core/configs#documentation), please use the custom config file [nf-core/configs/conf](https://github.com/nf-core/configs/tree/master/conf) that has been created for your institution.

   - For advanced Nextflow configuration options, see [_Nextflow configuration_](https://nf-co.re/docs/usage/configuration).

4. To test that everything is working properly, try running the tests for your pipeline of interest in the terminal:

   ```bash
   nextflow run nf-core/<pipeline_name> -profile test,docker --outdir <OUTDIR>
   ```

   - Replace `<pipeline_name>` with the name of an nf-core pipeline.

   - If you don't have Docker installed, replace `docker` in the command with either `singularity` or `conda`.

   - There is no need to download anything first — nextflow will pull the code from the GitHub repository and fetch the software requirements automatically.

   - If the pipeline fails, check the [troubleshooting docs](troubleshooting.md) and ask for help on the nf-core Slack channel for that particular pipeline (see [https://nf-co.re/join](https://nf-co.re/join)).

5. Read the pipeline documentation to see which command-line parameters are required. These will be specific to your data type and usage.

6. Launch the pipeline with some real data by omitting the `test` config profile and providing the required pipeline-specific parameters. For example, if you want to run the `methylseq` pipeline, you might use the following command:

   ```bash
   nextflow run nf-core/methylseq -profile docker --input 'input_data/*.fastq.gz' --outdir myproj/results --genome GRCh38
   ```

7. Once complete, check the pipeline execution and quality control reports. Each pipeline's documentation describes the different outputs to expect.

## Tips and tricks

- Hyphens matter! Core Nextflow command-line options use one (`-`), whereas pipeline-specific parameters use two (`--`)
- Specify `--email your@email.com` to receive emails when your pipeline completes
- Specify `--hook_url YOUR-HOOK-URL` or set the `hook_url.params` in `nextflow.config` to receive notifications from your pipeline in Microsoft Teams or Slack. Learn how to set up an incoming webhook in [MS Teams](https://learn.microsoft.com/en-us/microsoftteams/platform/webhooks-and-connectors/how-to/add-incoming-webhook?tabs=dotnet) and in [Slack](https://api.slack.com/messaging/webhooks).
- Include `-r <version-number>` when running your pipeline to specify a release version explicitly. That way, the same command will give identical results in future.
- Use `-resume` to restart pipelines that did not complete. This ensures that successful tasks from the previous run won't be re-executed.
- Use `nextflow log` to find names of all previous runs in your directory. These can be used with `-resume` to restart specific runs.
- Utilize multiple Nextflow configuration locations to your benefit. For example, use `-profile` for your cluster configuration, `~/.nextflow/config` for your personal configuration (with `params.email`, for example), and a working directory `nextflow.config` file for reproducible run-specific configuration.
- If you use Singularity, we recommend that you specify a cache directory with the [nextflow environment variable](https://www.nextflow.io/docs/latest/config.html#environment-variables) `NXF_SINGULARITY_CACHEDIR` in your `~./bash_profile` or `~/.bashrc` during the installation. This will store all your container images in one place, rather than downloading an image each time you run a pipeline. Only the base directory needs to be specified, Nextflow handles the folders and file names for you.

## Helper tools

To help you manage your nf-core pipelines and discover updates, we have written some command-line helper tools.
These allow you to list all available pipelines and versions, with information about the versions you're running locally.
There are also commands to download pipelines for offline use.

To find out more about these tools, see the [Tools](/tools) page.
