---
title: Getting started with nf-core
subtitle: Tutorial covering the basics of using nf-core pipelines.
---

> Material originally written for the Nextflow Camp 2019, Barcelona 2019-09-19: **_"Getting started with nf-core"_** _(see [programme](https://www.nextflow.io/nfcamp/2019/phil2.html))._
>
> Duration: **1hr**
>
> Updated for the nf-core Hackathon 2020, London 2020-03 _(see [event](https://nf-co.re/events#hackathon-francis-crick-2020))_.
>
> Updated for the Elixir workshop on November 2021 _(see [event](https://nf-co.re/events/2021/elixir-workflow-training-event))_.
>
> Updated during the March 2022 hackathon.

<!-- markdownlint-disable -->
<iframe src="https://www.slideshare.net/slideshow/embed_code/key/ewviLei9qkU7Ki?hostedIn=slideshare&page=upload" width="595" height="485" frameborder="0" marginwidth="0" marginheight="0" scrolling="no" style="border:1px solid #CCC; border-width:1px; margin-bottom:5px; max-width: 100%;" allowfullscreen> </iframe>
<!-- markdownlint-restore -->

**[Click here to download the slides associated with this tutorial.](/assets/markdown_assets/usage/nf_core_tutorial/nf-core-tutorial-usage.pdf)**

#### Overview

- [1 - Installation](#installing-the-nf-core-helper-tools)
- [2 - Listing pipelines](#listing-available-nf-core-pipelines)
- [3 - Running pipelines](#running-nf-core-pipelines)
- [4 - Troubleshooting pipelines](#troubleshooting-nf-core-pipelines)

## Abstract

The _nf-core_ community provides a range of tools to help new users get to grips with Nextflow - both by providing complete pipelines that can be used out of the box, and also by helping developers with best practices.
Companion tools can create a bare-bones pipeline from a template scattered with `TODO` pointers and CI with linting tools check code quality.
Guidelines and documentation help to get Nextflow newbies on their feet in no time.
Best of all, the _nf-core_ community is always on hand to help.

In this tutorial, we discuss the best-practice guidelines developed by the _nf-core_ community, why they're important and give insight into the best tips and tricks for budding Nextflow pipeline users. ✨

## Introduction

### What is _nf-core_

_nf-core_ is a community-led project to develop a set of best-practice pipelines built using Nextflow.
Pipelines are governed by a set of guidelines, enforced by community code reviews and automatic linting (code testing).
A suite of helper tools aim to help people run and develop pipelines.

### What this tutorial will cover

This tutorial attempts to give an overview of how _nf-core_ works:

- What are the most commonly used _nf-core_ tools.
- Listing pipelines available in the _nf-core_ project.
- How to run _nf-core_ pipelines.
- How to troubleshoot _nf-core_ pipelines.

### Where to get help

The beauty of _nf-core_ is that there is lots of help on offer!
The main place for this is Slack - an instant messaging service.

The _nf-core_ Slack can be found at [https://nfcore.slack.com](https://nfcore.slack.com/) (_NB: no hyphen in `nfcore`!_).
To join you will need an invite, which you can get at [https://nf-co.re/join/slack](https://nf-co.re/join/slack).

The _nf-core_ Slack organisation has channels dedicated for each pipeline, as well as specific topics (_eg._ [#help](https://nfcore.slack.com/channels/help), [#pipelines](https://nfcore.slack.com/channels/pipelines), [#tools](https://nfcore.slack.com/channels/tools), [#configs](https://nfcore.slack.com/channels/configs) and much more).

One additional tool which we like a lot is [TLDR](https://tldr.sh/) - it gives concise command line reference through example commands for most linux tools, including `nextflow`, `docker`, `singularity`, `conda`, `git` and more.
There are many clients, but [raylee/tldr](https://github.com/raylee/tldr) is arguably the simplest - just a single bash script.

## Installing the nf-core helper tools

Much of this tutorial will make use of the `nf-core` command line tool.
This has been developed to provide a range of additional functionality for the project such as pipeline creation, testing and more.

The `nf-core` tool is written in Python and is available from the [Python Package Index](https://pypi.org/project/nf-core/) and [Bioconda](https://bioconda.github.io/recipes/nf-core/README.html).
You can install the latest released version from PyPI as follows:

```bash
pip install nf-core
```

Or this command to install the `dev` version:

```bash
pip install --upgrade --force-reinstall git+https://github.com/nf-core/tools.git@dev
```

If using conda, first set up Bioconda as described in the [bioconda docs](https://bioconda.github.io/#usage) (especially setting the channel order), create and activate an environment and then install nf-core:

```bash
conda install nf-core
```

To update the package you can run the following command

```bash
conda update nf-core
```

The nf-core/tools source code is available at [https://github.com/nf-core/tools](https://github.com/nf-core/tools) - if you prefer, you can clone this repository and install the code locally:

```bash
git clone https://github.com/nf-core/tools.git nf-core-tools
cd nf-core-tools
python setup.py install
```

Once installed, you can check that everything is working by printing the help:

```bash
nf-core --help
```

You will also need to install [Prettier](https://prettier.io/) for formatting your code.
To do so, you can either use the following command with conda:

```bash
conda install prettier
```

Or use the Visual Studio Code extension [Prettier](https://marketplace.visualstudio.com/items?itemName=esbenp.prettier-vscode) also available in the pack of useful extension [NF-core](https://marketplace.visualstudio.com/items?itemName=nf-core.nf-core-extensionpack).

Besides, you can also add a comment with `@nf-core-bot fix linting` in your Pull Request and prettier will be used to apply the required fixes to your code.

### Exercise 1 (installation)

- Install nf-core/tools
- Use the help flag to list the available commands

## Listing available nf-core pipelines

As you saw from the `--help` output, the tool has a range of sub-commands.
The simplest is `nf-core list`, which lists all available _nf-core_ pipelines.
The output shows the latest version number, when that was released.
If the pipeline has been pulled locally using Nextflow, it tells you when that was and whether you have the latest version.

If you supply additional keywords after the command, the listed pipelines will be filtered.
Note that this searches more than just the displayed output, including keywords and description text.
The `--sort` flag allows you to sort the list (default is by most recently released) and `--json` returns the complete list, without any filtering, in JSON output for programmatic use.

The nf-core pipelines currently available and under development are also listed on the nf-core website, in the [pipelines page](https://nf-co.re/pipelines).

### Exercise 2 (listing pipelines)

- Use the help flag to print the list command usage
- List all available nf-core pipelines
- Sort pipelines alphabetically, then by popularity (stars)
- Fetch one of the pipelines using `nextflow pull`
- Use `nf-core list` to see if the pipeline you pulled is up to date
- Filter pipelines for those that work with RNA
- Save these pipeline details to a JSON file

## Running nf-core pipelines

### Software requirements for nf-core pipelines

In order to run _nf-core_ pipelines, you will need to have Nextflow installed ([https://www.nextflow.io](https://www.nextflow.io/)).
The only other requirement is a software packaging tool: [Conda](https://docs.conda.io/en/latest/miniconda.html), [Docker](https://www.docker.com) or [Singularity](https://sylabs.io/singularity/).
In theory it is possible to run the pipelines with software installed by other methods (_e.g._ environment modules, or manual installation), but this is not recommended.
Most people find either Docker or Singularity containers the best options, as conda environments cannot guarantee 100% reproducibility.

### Fetching pipeline code

Unless you are actively developing pipeline code, we recommend using the Nextflow [built-in functionality](https://www.nextflow.io/docs/latest/sharing.html) to fetch _nf-core_ pipelines.
Nextflow will automatically fetch the pipeline code when you run `nextflow run nf-core/PIPELINE`.
For the best reproducibility, it is good to explicitly reference the pipeline version number that you wish to use with the `-revision`/`-r` flag.
For example:

```bash
nextflow run nf-core/rnaseq -revision 3.4
```

If not specified, Nextflow will fetch the default branch.
For pipelines with a stable release this the default branch is `master` - this branch contains code from the latest release.
For pipelines in early development that don't have any releases, the default branch is `dev`.

If you would like to run the latest development code, use `-r dev`.

Note that once pulled, Nextflow will use the local cached version for subsequent runs.
Use the `-latest` flag when running the pipeline to always fetch the latest version.
Alternatively, you can force Nextflow to pull a pipeline again using the `nextflow pull` command:

```bash
nextflow pull nf-core/rnaseq -revision 3.4
```

### Usage instructions and documentation

You can find general documentation and instructions for Nextflow and _nf-core_ on the _nf-core_ website: [https://nf-co.re/](https://nf-co.re/).
Pipeline-specific documentation is bundled with each pipeline in the `/docs` folder.
This can be read either locally, on GitHub, or on the _nf-core_ website.
Each pipeline has its own webpage at `https://nf-co.re/<pipeline_name>` (_e.g._ [nf-co.re/rnaseq](https://nf-co.re/rnaseq)), including `Usage` documentation, `Output` documentation and `Parameter` documentation.

In addition to this documentation, each pipeline comes with basic command line reference.
This can be seen by running the pipeline with the `--help` flag, for example:

```bash
nextflow run nf-core/rnaseq --help
```

Example results of a pipeline run on full-sized test data can be browsed on the pipeline page, under the `aws` results tab.

### Config profiles

Nextflow can load pipeline configurations from multiple locations.
To make it easy to apply a group of options on the command line, Nextflow uses the concept of [config profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles).
_nf-core_ pipelines load configuration in the following order:

1. Pipeline: Default 'base' config
   - Always loaded. Contains pipeline-specific parameters and "sensible defaults" for things like computational requirements
   - Does _not_ specify any method for software packaging. If nothing else is specified, Nextflow will expect all software to be available on the command line.
2. Pipeline: Core config profiles
   - All _nf-core_ pipelines come with some generic config profiles. The most commonly used ones are for software packaging: `docker`, `singularity` and `conda`. To ensure reproducibility across different compute infrastructures, it is recommended to use containers instead of conda environments.
   - Other core profiles are `debug` and `test`
3. [nf-core/configs](https://github.com/nf-core/configs): Server profiles
   - At run time, _nf-core_ pipelines fetch configuration profiles from the [configs](https://github.com/nf-core/configs) remote repository. The profiles here are specific to clusters at different institutions.
   - Because this is loaded at run time, anyone can add a profile here for their system and it will be immediately available for all _nf-core_ pipelines.
4. Personal configuration under `~/.nextflow/config`.
5. Local config files given to Nextflow with the `-c` flag.
6. Command line configuration.

Multiple comma-separate config profiles can be specified in one go, so the following commands are perfectly valid:

```bash
nextflow run nf-core/rnaseq -profile test,docker
nextflow run nf-core/rnaseq -profile singularity,debug
```

Note that the order in which config profiles are specified matters.
Their priority increases from left to right.

> Our tip: Be clever with multiple Nextflow configuration locations. For example, use `-profile` for your cluster configuration, `~/.nextflow/config` for your personal config such as `params.email` and a working directory `config` (e.g. `custom.config` provided to the run with `-c custom.config`) file for reproducible run-specific configuration.

To know more about Nextflow configurations you can check the [pipeline configuration tutorial](https://nf-co.re/docs/usage/configuration).

### Running pipelines with test data

The `test` config profile is a bit of a special case.
Whereas all other config profiles tell Nextflow how to run on different computational systems, the `test` profile configures each `nf-core` pipeline to run without _any_ other command line flags.
It specifies URLs for test data and all required parameters.
Because of this, you can test any _nf-core_ pipeline with the following command:

```bash
nextflow run nf-core/<pipeline_name> -profile test --outdir <OUTDIR>
```

> Note that you will typically still need to combine this with a configuration profile for your system - _e.g._ `-profile test,docker`.
> Running with the test profile is a great way to confirm that you have Nextflow configured properly for your system before attempting to run with real data.

### The nf-core launch command

Most _nf-core_ pipelines have a number of flags that need to be passed on the command line: some mandatory, some optional.
To make it easier to launch pipelines, these parameters are described in a JSON file bundled with the pipeline.
The `nf-core launch` command uses this to build an interactive command-line wizard which walks through the different options with descriptions of each, showing the default value and prompting for values.

Once all prompts have been answered, non-default values are saved to a `params.json` file which can be supplied to Nextflow to run the pipeline. Optionally, the Nextflow command can be launched there and then.

To use the launch feature, just specify the pipeline name:

```bash
nf-core launch <pipeline_name>
```

### Using nf-core pipelines offline

Many of the techniques and resources described above require an active internet connection at run time - pipeline files, configuration profiles and software containers are all dynamically fetched when the pipeline is launched.
This can be a problem for people using secure computing resources that do not have connections to the internet.

To help with this, the `nf-core download` command automates the fetching of required files for running _nf-core_ pipelines offline.
The command can download a specific release of a pipeline with `-r`/`--release` and fetch the singularity container if `--singularity` is passed (this needs Singularity to be installed).
All files are saved to a single directory, ready to be transferred to the cluster where the pipeline will be executed.

To know more about running pipelines offline you can check the [pipeline configuration tutorial](https://nf-co.re/docs/usage/offline).

### Exercise 3 (using pipelines)

- Install required dependencies (Nextflow, Docker)
- Print the command-line usage instructions for the _nf-core/rnaseq_ pipeline
- In a new directory, run the _nf-core/rnaseq_ pipeline with the provided test data
- Try launching the RNA pipeline using the `nf-core launch` command
- Download the _nf-core/rnaseq_ pipeline for offline use using the `nf-core download` command

## Troubleshooting nf-core pipelines

Not everything always runs smoothly and you might be getting some errors when running nf-core pipelines. Here are some step-by-step tips that can help you troubleshoot your errors.

1. Start small: each nf-core pipeline comes with small test data that are checked by continuous integration and for each pipeline release.
   - Start by running the pipeline tests as described [above](#running-pipelines-with-test-data). If these tests fail, there is a good chance that you are missing some of the components needed to run Nextflow pipelines.
   - Nextflow: check that you have the latest version installed.
   - Check that you have docker/singularity/conda installed and that you are using the right docker/singularity/conda/custom profile.
   - Check the [troubleshooting docs](https://nf-co.re/docs/usage/troubleshooting).
2. Categorize the type of error. Check the Nextflow low to figure out if the error occurs:
   - Before the first process
   - In the first process
   - During the pipeline run
   - Problems with the process output
3. Read the Nextflow log. Check the work directory for the `.command.err` or `.command.log` files for more information.
4. Search the nf-core slack, google. Ask for help in the corresponding nf-core slack channel.
5. Report a pipeline bug on the nf-core GitHub if none of the above steps helps.

Here is a bytesize talk with a step by step explanation on how to troubleshoot failing pipelines.

<!-- markdownlint-disable -->
<div class="ratio ratio-16x9">
<iframe src="https://www.youtube.com/embed/z9n2F4ByIkY" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
</div>
<!-- markdownlint-restore -->

## Conclusion

We hope that this _nf-core_ tutorial has been helpful!
Remember that there is more in-depth documentation on many of these topics available on the [nf-core website](https://nf-co.re).
If in doubt, please ask for help on Slack.

If you have any suggestions for how to improve this tutorial, or spot any mistakes, please create an issue or pull request on the [nf-core/nf-co.re repository](https://github.com/nf-core/nf-co.re).

> [Phil Ewels](https://github.com/ewels/), [Maxime Garcia](https://github.com/MaxUlysse), [Gisela Gabernet](https://github.com/ggabernet), [Friederike Hanssen](https://github.com/FriederikeHanssen) for _nf-core_, March 2022
