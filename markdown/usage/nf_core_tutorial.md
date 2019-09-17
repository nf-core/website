# _nf-core_ Tutorial

> Originally written for the Nextflow Camp 2019, Barcelona 2019-09-19: ***"Getting started with nf-core"*** _(see [programme](https://www.nextflow.io/nfcamp/2019/phil2.html))_
>
> Duration: **1hr 45**

> Last updated: September  2019

## Table of Contents

#### Sections
* [Abstract](#abstract)
* [Introduction](#introduction)
* [Installing the _nf-core_ helper tools](#installing-the-nf-core-helper-tools)
* [Listing available _nf-core_ pipelines](#listing-available-nf-core-pipelines)
* [Running _nf-core_ pipelines](#running-nf-core-pipelines)
* [Creating _nf-core_ pipelines](#creating-nf-core-pipelines)
* [Testing _nf-core_ pipelines](#testing-nf-core-pipelines)
* [Releasing _nf-core_ pipelines](#releasing-nf-core-pipelines)

#### Exercises
* [1 - installation](#exercise-1-installation)
* [2 - listing pipelines](#exercise-2-listing-pipelines)
* [3 - using pipelines](#exercise-3-using-pipelines)
* [4 - creating pipelines](#exercise-4-creating-pipelines)
* [5 - testing pipelines](#exercise-5-testing-pipelines)
* [6 - releasing pipelines](#exercise-6-releasing-pipelines)

## Abstract

The _nf-core_ community provides a range of tools to help new users get to grips with nextflow - both by providing complete pipelines that can be used out of the box, and also by helping developers with best practices. Companion tools can create a bare-bones pipeline from a template scattered with `TO-DO` pointers and CI with linting tools check code quality. Guidelines and documentation help to get nextflow newbies on their feet in no time. Best of all, the _nf-core_ community is always on hand to help.

In this tutorial we discuss the best-practice guidelines developed by the _nf-core_ community, why they're important and give insight into the best tips and tricks for budding nextflow pipeline developers. ✨

## Introduction
### What is _nf-core_
_nf-core_ is a community-led project to develop a set of best-practice pipelines built using Nextflow. Pipelines are governed by a set of guidelines, enforced by community code reviews and automatic linting (code testing). A suite of helper tools aim to help people run and develop pipelines.

### What this tutorial will cover
This tutorial attempts to give an overview of how _nf-core_ works: how to run _nf-core_ pipelines, how to make new pipelines using the _nf-core_ template and how _nf-core_ pipelines are reviewed and ultimately released.

### Where to get help
The beauty of _nf-core_ is that there is lots of help on offer! The main place for this is Slack - an instant messaging service. The _nf-core_ Slack organisation has channels dedicated for each pipeline, as well as specific topics (_eg._ `#new-pipeliens`, `#tools` and `#aws` ).

The _nf-core_  Slack can be found at [https://nfcore.slack.com](https://nfcore.slack.com/) (_NB: no hyphen in `nfcore`!_). To join you will need an invite, which you can get at [https://nf-co.re/join/slack](https://nf-co.re/join/slack).

One additional tool which this author swears by is [TLDR](https://tldr.sh/) - it gives concise command line reference through example commands for most linux tools, including `nextflow`, `docker`, `singularity`,  `conda`, `git` and more. There are many clients, but [raylee/tldr](https://github.com/raylee/tldr) is arguably the simplest - just a single bash script.

## Installing the _nf-core_ helper tools
Much of this tutorial will make use of the `nf-core` command line tool. This has been developed to provide a range of additional functionality for the project such as pipeline creation, testing and more.

The `nf-core` tool is written in Python and is available from the [Python Package Index](https://pypi.org/project/nf-core/) and [Bioconda](https://bioconda.github.io/recipes/nf-core/README.html). You can install it from either of these channels as follows:

```bash
pip install nf-core
conda install -c bioconda nf-core
```

The source code is available at [https://github.com/nf-core/tools](https://github.com/nf-core/tools) - if you prefer, you can clone this repository and install the code locally:

```bash
git clone https://github.com/nf-core/tools.git nf-core-tools
cd nf-core-tools
python setup.py install
```

Once installed, you can check that everything is working by printing the help:

```console
$ nf-core --help
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

Usage: nf-core [OPTIONS] COMMAND [ARGS]...

Options:
  --version      Show the version and exit.
  -v, --verbose  Verbose output (print debug statements)
  --help         Show this message and exit.

Commands:
  list          List nf-core pipelines with local info
  launch        Run pipeline, interactive parameter prompts
  download      Download a pipeline and singularity container
  licences      List software licences for a given workflow
  create        Create a new pipeline using the template
  lint          Check pipeline against nf-core guidelines
  bump-version  Update nf-core pipeline version number
```

###  Exercise 1 (installation)
* Install nf-core/tools
* Use the help flag to list the available commands

## Listing available _nf-core_ pipelines
As you saw from the `--help` output, the tool has a range of subcommands. The simplest is `nf-core list`, which lists all available _nf-core_ pipelines. The output shows the latest version number, when that was released. If the pipeline has been pulled locally using Nextflow, it tells you when that was and whether you have the latest version.

If you supply additional keywords after the command, it will search the pipelines for that. Note that it searches more than just the displayed output, including keywords and description text. The `--sort` flag allows you to sort the list (default is by most recently released) and `--json` gives JSON  output for programmatic use.

###  Exercise 2 (listing pipelines)
* Use the help flag to print the list command usage
* List all pipelines
* Sort pipelines alphabetically, then by popularity (stars)
* Fetch one of the pipelines using `nextflow pull`
* Use `nf-core list` to see if the pipeline you pulled is up to date
* Filter pipelines for those that work with RNA
* Save these pipeline details to a JSON file

## Running _nf-core_ pipelines
### Software requirements for _nf-core_ pipelines
In order to run _nf-core_ pipelines, you will need to have Nextflow installed ([https://www.nextflow.io](https://www.nextflow.io/)). The only other requirement is a software packaging tool: [Conda](https://docs.conda.io/en/latest/miniconda.html), [Docker](https://www.docker.com) or [Singularity](https://sylabs.io/singularity/). In theory it is possible to run the pipelines with software installed by other methods (_e.g._ environment modules, or manual installation), but this is not recommended. Most people find either Docker or Singularity the best options.

### Fetching pipeline code
Unless you are actively developing pipeline code, we recommend using the Nextflow [built-in functionality](https://www.nextflow.io/docs/latest/sharing.html) to fetch _nf-core_ pipelines. Nextflow will automatically fetch the pipeline code when you run `nextflow run nf-core/PIPELINE-NAME`. For the best reproducibility, it is good to explicitly reference the pipeline version number that you wish to use with the `-revision`/`-r` flag. For example:

```bash
nextflow run nf-core/rnaseq -revision 1.3
```

If not specified, Nextflow will fetch the `master` branch - for _nf-core_ pipelines this will be the latest release. If you would like to run the latest development code, use `-r dev`.

Note that once pulled, Nextflow will use the local cached version for subsequent runs. Use the `-latest` flag when running the pipeline to always fetch the latest version. Alternatively, you can force Nextflow to pull a pipeline again using the `nextflow pull` command:

```bash
nextflow pull nf-core/rnaseq
```

### Usage instructions and documentation
You can find general documentation and instructions for Nextflow and _nf-core_ on the _nf-core_ website: [https://nf-co.re/](https://nf-co.re/). Pipeline-specific documentation is bundled with each pipeline in the `/docs` folder. This can be read either locally, on GitHub, or on the _nf-core_ website. Each pipeline has its own webpage at `https://nf-co.re/PIPELINE-NAME`.

In addition to this documentation, each pipeline comes with basic command line reference. This can be seen by running the pipeline with the `--help` flag, for example:

```bash
nextflow run nf-core/rnaseq --help
```

### Config profiles
Nextflow can load pipeline configurations from multiple locations. To make it easy to apply a group of options on the command line, Nextflow uses the concept of [config profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles).  _nf-core_ pipelines load configuration in the following order:

1. Pipeline: Default 'base' config
  * Always loaded. Contains pipeline-specific parameters and "sensible defaults" for things like computational requirements
  * Does _not_ specify any method for software packaging. If nothing else is specified, Nextflow will expect all software to be available on the command line.
2. Pipeline: Core config profiles
  * All _nf-core_ pipelines come with some generic config profiles. The most commonly used ones are for software packaging: `docker`, `singularity` and `conda`
  * Other core profiles are `awsbatch`, `debug` and `test`
3. [nf-core/configs](https://github.com/nf-core/configs): Server profiles
  * At run time, _nf-core_ pipelines fetch configuration profiles from the [configs](https://github.com/nf-core/configs) remote repository. The profiles here are specific to clusters at different institutions.
  * Because this is loaded at run time, anyone can add a profile here for their system and it will be immediately available for all _nf-core_ pipelines.
4. Local config files given to Nextflow with the `-c` flag
5. Command line configuration

Multiple comma-separate config profiles can be specified in one go, so the following commands are perfectly valid:

```bash
nextflow run nf-core/rnaseq -profile test,docker
nextflow run nf-core/hlatyping -profile singularity,debug
```

Note that the order in which config profiles are specified matters...

### Running pipelines with test data
The `test` config profile is a bit of a special case. Whereas all other config profiles tell Nextflow how to run on different computational systems, the `test` profile configures each `nf-core` pipeline to run without _any_ other command line flags. It specifies URLs for test data and all required parameters. Because of this, you can test any _nf-core_ pipeline with the following command:

```bash
nextflow run nf-core/PIPELINE -profile test
```

Note that you will typically still need to combine this with a configuration profile for your system - _e.g._ `-profile test,docker`.

### The _nf-core_ launch command
### Using _nf-core_ pipelines offline
### Exercise 3 (using pipelines)
* Install other dependencies (nextflow, docker)
* Print the command-line usage instructions for the nf-core/rnaseq pipeline
* In a new directory, run the nf-core/rnaseq pipeline with the provided test data
* Try launching the RNA pipeline using the `nf-core launch` command
* Download the nf-core/rnaseq pipeline for offline use using the `nf-core download` command
## Creating _nf-core_ pipelines
### Using the nf-core template
### `TODO` statements
### How _nf-core_ software packaging works
### Forks, branches and pull-requests
Common difficulties:
* Changelog must be updated
* PRs must be against `dev` branch
* PRs need to be reviewed before being merged _(after first release)_
### Setting up Docker and Travis CI
### Exercise 4 (creating pipelines)
* Make a new pipeline using the template
* Update the readme file to fill in the `TODO` statements
* Add a new process to the pipeline in `main.nf`
* Add the new software dependencies from this process in to `environment.yaml`
## Testing _nf-core_ pipelines
###  Linting _nf-core_ pipelines
### Choosing test data
### nf-core/test_datasets
### Travis CI configuration
### Exercise 5 (testing pipelines)
* Run `nf-core lint` on your pipeline and make note of any test warnings / failures
* Read up on one or two of the linting rules on the nf-core website and see if you can fix some.
* Take a look at `conf/test.config` and switch the test data for another dataset on nf-core/test_data
## Releasing _nf-core_ pipelines
### Forking to _nf-core_
### Initial community review
### Making the first release
### Template updates
### Pipeline version numbers
### Review releases
### Exercise 6 (releasing pipelines)
* Use `nf-core bump-version` to update the required version of Nextflow in your pipeline
* Bump your pipeline's version to 1.0, ready for its first release!
* Make sure that you're signed up to the _nf-core_ slack (get an invite on [nf-co.re](https://nf-cor.re)) and drop us a line about your latest and greatest pipeline plans!
* Ask to be a member of the _nf-core_ GitHub organisation by commenting on [this GitHub issue](https://github.com/nf-core/nf-co.re/issues/3)
* If you're a twitter user, make sure to follow the [@nf_core](https://twitter.com/nf_core) account
