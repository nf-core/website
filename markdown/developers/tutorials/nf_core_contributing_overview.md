---
title: nf-core tutorial
subtitle: Tutorial covering the basics of contributing to nf-core.
menu:
  main:
    weight: 10
---

> Material originally written for the Nextflow Camp 2019, Barcelona 2019-09-19: **_"Getting started with nf-core"_** _(see [programme](https://www.nextflow.io/nfcamp/2019/phil2.html))._
>
> Duration: **1hr 30min**
>
> Updated for the nf-core Hackathon 2020, London 2020-03 _(see [event](https://nf-co.re/events#hackathon-francis-crick-2020))_.
>
> Updated for the Elixir workshop on November 2021 _(see [event](https://nf-co.re/events/2021/elixir-workflow-training-event))_.
>
> Updated during the March 2022 hackathon.

<iframe src="https://www.slideshare.net/slideshow/embed_code/key/v7hOzFwZTgD7o0?hostedIn=slideshare&page=upload" width="595" height="485" frameborder="0" marginwidth="0" marginheight="0" scrolling="no" style="border:1px solid #CCC; border-width:1px; margin-bottom:5px; max-width: 100%;" allowfullscreen> </iframe>

**[Click here to download the slides associated with this tutorial.](/assets/markdown_assets/usage/nf_core_tutorial/nf-core-tutorial-contributing.pdf)**

#### Overview

- [1 - Installation](#installing-the-nf-core-helper-tools)
- [2 - Creating pipelines](#creating-nf-core-pipelines)
- [3 - Testing pipelines](#testing-nf-core-pipelines)
- [4 - nf-core modules](#nf-core-modules)
- [5 - Releasing pipelines](#releasing-pipelines)

## Abstract

The _nf-core_ community provides a range of tools to help new users get to grips with Nextflow - both by providing complete pipelines that can be used out of the box, and also by helping developers with best practices.
Companion tools can create a bare-bones pipeline from a template scattered with `TODO` pointers and CI with linting tools check code quality.
Guidelines and documentation help to get Nextflow newbies on their feet in no time.
Best of all, the _nf-core_ community is always on hand to help.

In this tutorial we discuss the best-practice guidelines developed by the _nf-core_ community, why they're important and give insight into the best tips and tricks for budding Nextflow pipeline developers. ✨

## Introduction

### What is _nf-core_

_nf-core_ is a community-led project to develop a set of best-practice pipelines built using Nextflow.

The aim of _nf-core_ is to develop with the community over GitHub and other communication channels to build a pipeline to solve a particular analysis task.
We therefore encourage cooperation and adding new features to existing pipelines.
Before considering starting a new pipeline, have a look if a pipeline already exists performing a similar task and consider contributing to that one instead.
If there is no pipeline for the analysis task at hand, let us know about your new pipeline plans on Slack in the [`#new-pipelines`](https://nfcore.slack.com/channels/new-pipelines) channel.

Pipelines are governed by a set of guidelines, enforced by community code reviews and automatic linting (code testing).
A suite of helper tools aim to help people run and develop pipelines.

### What this tutorial will cover

This tutorial attempts to give an overview of:

- creating new pipelines using the _nf-core_ template.
- what nf-core shared modules are.
- adding nf-core shared modules to a pipeline.
- creating new nf-core modules using the _nf-core_ module template.
- reviewing and releasing _nf-core_ pipelines.

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

If using conda, first set up Bioconda as described in the [bioconda docs](https://bioconda.github.io/#usage) (especially setting the channel order) and then install nf-core:

```bash
conda install nf-core
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

### Exercise 1 (installation)

- Install nf-core/tools
- Use the help flag to list the available commands

## Creating nf-core pipelines

### Using the nf-core template

The heart of _nf-core_ is the standardisation of pipeline code structure.
To achieve this, all pipelines adhere to a generalised pipeline template.
The best way to build an _nf-core_ pipeline is to start by using this template via the `nf-core create` command.
This launches an interactive prompt on the command line which asks for things such as pipeline name, a short description and the author's name.
These values are then propagated throughout the template files automatically.

> Contribution guidelines: one of the main ideas of nf-core is to develop with the community to build together best-practise analysis pipelines.
> We encourage cooperation rather than duplication, and contributing to and extending existing pipelines that might be performing similar tasks.
> For more details about this, please check out the [contribution guidelines](https://nf-co.re/docs/contributing/guidelines).

### TODO statements

Not everything can be completed with a template and all new pipelines will need to edit and add to the resulting pipeline files in a similar set of locations.
To make it easier to find these, the _nf-core_ template files have numerous comment lines beginning with `TODO nf-core:`, followed by a description of what should be changed or added.
These comment lines can be deleted once the required change has been made.

Most code editors have tools to automatically discover such `TODO` lines and the `nf-core lint` command will flag these.
This makes it simple to systematically work through the new pipeline, editing all files where required.

### Forks, branches and pull-requests

All _nf-core_ pipelines use GitHub as their code repository, and git as their version control system.
For newcomers to this world, it is helpful to know some of the basic terminology used:

- A _repository_ contains everything for a given project
- _Commits_ are code checkpoints.
- A _branch_ is a linear string of commits - multiple parallel branches can be created in a repository
- Commits from one branch can be _merged_ into another
- Repositories can be _forked_ from one GitHub user to another
- Branches from different forks can be merged via a _Pull Request_ (PR) on [github.com](https://github.com)

Typically, people will start developing a new pipeline under their own personal account on GitHub.
When it is ready for its first release and has been discussed on Slack, this repository is forked to the _nf-core_ organisation.
All developers then maintain their own forks of this repository, contributing new code back to the _nf-core_ fork via pull requests.

All _nf-core_ pipelines must have the following three branches:

1. `master` - commits from stable releases _only_. Should always have code from the most recent release.
2. `dev` - current development code. Merged into `master` for releases.
3. `TEMPLATE` - used for template automation by the [@nf-core-bot](https://github.com/nf-core-bot) GitHub account. Should only contain commits with unmodified template code.

Pull requests to the _nf-core_ fork have a number of automated steps that must pass before the PR can be merged.
A few points to remember are:

- The pipeline `CHANGELOG.md` must be updated
- PRs must not be against the `master` branch (typically you want `dev`)
- PRs should be reviewed by someone else before being merged (help can be asked on the [#request-review](https://nfcore.slack.com/channels/request-review) Slack channel)

### Exercise 2 (creating pipelines)

- Make a new pipeline using the template
- Update the readme file to fill in the `TODO` statements

## Testing nf-core pipelines

### Linting nf-core pipelines

Manually checking that a pipeline adheres to all _nf-core_ guidelines and requirements is a difficult job.
Wherever possible, we automate such code checks with a [code linter](<https://en.wikipedia.org/wiki/Lint_(software)>).
This runs through a series of tests and reports failures, warnings and passed tests.

The linting code is closely tied to the _nf-core_ template and both change over time.
When we change something in the template, we often add a test to the linter to make sure that pipelines do not use the old method.

Each lint test is documented on the [nf-core tools documentation website](https://nf-co.re/tools/docs/latest/pipeline_lint_tests/index.html).
When warnings and failures are reported on the command line, a short description is printed along with a link to the documentation for that specific test on the website.

Code linting is run automatically every time you push commits to GitHub, open a pull request or make a release.
You can also run these tests yourself locally with the following command:

```bash
nf-core lint /path/to/pipeline
```

When merging PRs from `dev` to `master`, the `lint` command will be run with the `--release` flag which includes a few additional tests.

### nf-core/test-datasets

When adding a new pipeline, you must also set up the `test` config profile.
To do this, we use the [nf-core/test-datasets](https://github.com/nf-core/test-datasets) repository.
Each pipeline has its own branch on this repository, meaning that the data can be cloned without having to fetch all test data for all pipelines:

```bash
git clone --single-branch --branch PIPELINE https://github.com/nf-core/test-datasets.git
```

To set up the test profile, make a new branch on the [nf-core/test-datasets](https://github.com/nf-core/test-datasets) repo through the web page (see [instructions](https://help.github.com/en/articles/creating-and-deleting-branches-within-your-repository#creating-a-branch)).
Fork the repository to your user and open a PR to your new branch with a really (really!) tiny dataset.
Once merged, set up the `conf/test.config` file in your pipeline to refer to the URLs for your test data.

These test datasets are used by the automated continuous integration tests.
The systems that run these tests are _extremely_ limited in the resources that they have available.
Typically, the pipeline should be able to complete in around 10 minutes and use no more than 6-7 GB memory.
To achieve this, input files and reference genomes need to be very tiny.
If possible, a good approach can be to use PhiX or Yeast as a reference genome.
Alternatively, a single small chromosome (or part of a chromosome) can be used.
If you are struggling to get the tests to run, ask for help on Slack.

When writing `conf/test.config` remember to define all required parameters so that the pipeline will run with only `-profile test`
(with the exception of `params.outdir`).
Note that remote URLs cannot be traversed like a regular file system - so glob file expansions such as `*.fa` will not work.

### GitHub Actions workflows

The automated tests with GitHub Actions are configured in the `.github/workflows` folder that is generated by the template.
Each file (`branch.yml`, `ci.yml` and `linting.yml`) defines several tests: branch protection, running the pipeline with the test data, linting the code with `nf-core lint`, linting the syntax of all Markdown documentation and linting the syntax of all YAML files.

#### Branch protection

The `branch.yml` workflow sets the branch protection for the nf-core repository.
It is already set up and does not need to be edited.
This test will check that pull-requests going to the nf-core repo `master` branch only come from the `dev` or `patch` branches (for releases).
In case you want to add branch protection for a repository outside _nf-core_, you can add an extra step to the workflow with that repository name.
You can leave the _nf-core_ branch protection in place so that the `nf-core lint` command does not throw an error - it only runs on nf-core repositories so it should be ignored.

#### Linting

The `linting.yml` workflow will run the nf-core linting tests, as well as other linting tests that ensure that general code styling is maintained.

This workflow is already set up and does not need to be edited.
It will check:

- That markdown / YAML / JSON etc files are properly formatted, using Prettier.
- That any Python is correctly formatted using Black.
- That general code structure is maintained (e.g. no trailing whitespaces, final newline is present in all files) with [EditorConfig](https://editorconfig.org/).
- That the pipeline code follows the nf-core linting rules.
  - You can run this manually with the `nf-core lint` command.
  - Some of the linting tests can be fixed automatically with the `--fix` subcommand, for example any `files_unchanged` type of error: `nf-core lint --fix files_unchanged`.

See the documentation about [code formatting](/developers/code_formatting.md) for more information on these tools and how to use them.

#### Prettier

Prettier is a tool used to make the code more consistent and readable.
It's a optionated code formatter that enforces a consistent style by parsing the code and re-printing it in its own rules.

Every code and file of the nf-core community needs to be checked and formatted by Prettier before the upload.

The command to check if the style of all files in your directory and its subdirectories is valid for Prettier is the following:

```bash
prettier -c .
```

While the command to automatically format the style of all the files in your directory with Prettier is the following:

```bash
prettier -w .
```

See pre-commit informations here: https://nf-co.re/docs/contributing/modules#new-module-workflow

#### Continuous integration tests

The `ci.yml` workflow will run tests for your pipeline with the small test dataset.
It should mostly work out of the box, but it may need some editing from you to customise it for your pipeline.

The matrix `nxf_ver` variable sets the `NXF_VER` environment variable twice.
This tells GitHub Actions to run the tests twice in parallel - once with the latest version of Nextflow (`NXF_VER=''`) and once with the minimum version supported by the pipeline.
Do not edit this version number manually - it appears in multiple locations through the pipeline code, so it's better to use `nf-core bump-version --nextflow` instead.

The provided tests run your pipeline with the `-profile test,docker --outdir <OUTDIR>` flags. This may be sufficient for your pipeline.
However, if it is possible to run the pipeline with significantly different options (for example, different alignment tools), then it is good to test all of these.
You can do this by adding additional tests in the `jobs` block.
Do not try to add a run for every possible combination of parameters, as this would take too long to run.

### Exercise 3 (testing pipelines)

- Run `nf-core lint` on your pipeline and make note of any test warnings / failures.
- Take a look at `conf/test.config` and edit the test data (you can use another dataset on nf-core/test-datasets).
- Check that running the pipeline with the test profile passes

## nf-core modules

The Nextflow DSL2 syntax allows the modularizing of Nextflow pipelines, so workflows, subworkflows and modules can be defined and imported into a pipeline.
This allows for the sharing of pipeline processes (modules, and also routine subworkflows) among nf-core pipelines.

Shared modules are stored in the [nf-core/modules](https://github.com/nf-core/modules) repository.
Modules on this repository are as atomic as possible, in general calling each one tool only.
If a tool consists of several subtools (e.g. `bwa index` and `bwa mem`), these will be stored in individual modules with the naming convention `tool/subtool`.
Each module defines the input and output channels, the process script, as well as the software packaging for a specific process. Conda environments, docker or singularity containers are defined within each module. We mostly rely on the [biocontainers](https://biocontainers.pro/) project for providing single-tool containers for each module.

The nf-core/modules repository also includes defined tests for each module which run on tiny test data on the nf-core/test-datasets repository (modules branch). The modules tests run in a similar way as pipeline tests on GitHub actions and ensure that modules are always functional and produce the desired results.

nf-core tools have a series of subcommands to work with nf-core modules.

```bash
nf-core modules --help
```

You can list all the modules that are available in the nf-core/modules repository with the following command:

```bash
nf-core modules list remote
```

You can filter the search by typing the name of a tool or part of it.

### Adding nf-core modules to a pipeline

Adding nf-core modules to a pipeline, if the modules already exist in the nf-core modules repository, can be done with the following command (executing it in the main pipeline directory):

```bash
nf-core modules install <module name>
```

The modules files will be added under the `modules/nf-core` directory. To be able to call the module inside the main pipeline workflow (such as `workflows/<pipeline-name>.nf`) or a subworkflow, an include statement needs to be added in the corresponding Nextflow file:

```bash
include { TOOL_SUBTOOL } from '../modules/nf-core/modules/<tool/subtool>/main'
```

Tool options or other options that should be passed to the module can be defined in the `modules.config` configuration file.

If the module defines additional parameters, these parameters are directly available to the module, and their default should be specified in the `conf/modules.config` file (a warning `Access to undefined parameter` will otherwise be thrown). These parameters should also be included in the `nextflow_schema.json` parameter definition file, which can be done interactively by using the nf-core tool:

```bash
nf-core schema build
```

### Adding modules to nf-core

In the same way as pipelines, nf-core modules are created from a template using the `nf-core modules create` command.

A step-by-step tutorial [adding new modules](https://nf-co.re/docs/contributing/tutorials/dsl2_modules_tutorial) has been created that explains the procedure in detail.

### Exercise 4 (add a module to a pipeline)

- Use the pipeline you created in Exercise 2 and add an nf-core module (e.g. trimgalore).
- Connect the module to the main pipeline workflow (under `workflows/<pipeline-name>.nf`).

#### Subworkflow

Subworkflows are chains of modules that offer a higher-level of functionality within the context of a pipeline. A subworkflow is able to run multiple tools on the same input file.
Subworkflows should be shipped with the pipeline implementation and if required they should be shared amongst different pipelines directly from there. Shareable nf-core subworkflow files are available in the 'subworkflow/' directory of nf-core/modules along with the required documentation and tests.

## Releasing nf-core pipelines

Your pipeline is written and ready to go!
Before you can release it with _nf-core_ there are a few steps that need to be done.
First, tell everyone about it on Slack in the [`#new-pipelines`](https://nfcore.slack.com/channels/new-pipelines) channel.
Hopefully you've already done this before you spent lots of time on your pipeline, to check that there aren't other similar efforts happening elsewhere.
Next, you need to be a member of the [nf-core GitHub organisation](https://github.com/nf-core/).
You can find instructions for how to do this at [https://nf-co.re/join](https://nf-co.re/join).

### Forking to nf-core

Once you're ready to go, you can fork your repository to _nf-core_.
A lot of stuff happens automatically when you do this: the website will update itself to include your new pipeline, complete with rendered documentation pages and usage statistics.
Your pipeline will also appear in the `nf-core list` command output and in various other locations.

Unfortunately, at the time of writing, the Zenodo (automated DOI assignment for releases) services are not created automatically.
These can only be set up by _nf-core_ administrators, so please ask someone to do this for you on Slack.

### Initial community review

Once everything is set up and all tests are passing on the `dev` branch, let us know on Slack and we will do a large community review.
This is a one-off process that is done before the first release for all pipelines.
In order to give a nice interface to review _all_ pipeline code, we create a _"pseudo pull request"_ comparing `dev` against the first commit in the pipeline (hopefully the template creation).
This PR will never be merged, but gives the GitHub review web pages where people can comment on specific lines in the code.

These first community reviews can take quite a long time and typically result in a lot of comments and suggestions (nf-core/deepvariant famously had 156 comments before it was approved).
Try not to be intimidated - this is the main step where the community attempts to standardise and suggest improvements for your code.
Your pipeline will come out the other side stronger than ever!

### Making the first release

Once the pseudo-PR is approved, you're ready to make the release.
To do this, first bump the pipeline version to a stable tag using `nf-core bump-version`, then open a pull-request from the `dev` branch to `master`.
Once tests are passing and two _nf-core_ members have approved this PR, it can be merged to `master`.
Then a GitHub release is made, using the contents of the changelog as a description.

Pipeline version numbers (release tags) should be numerical only, using [semantic versioning](https://semver.org/spec/v2.0.0.html).
For example, with a release version `1.4.3`, bumping to `2.0` would correspond to the _major_ release where results would no longer be backwards compatible.
Bumping to `1.5` would be a minor release, for example adding some new features.
Bumping to `1.4.4` would be a patch release for minor things such as fixing bugs.

We have compiled a detailed [release checklist](https://nf-co.re/docs/contributing/release_checklist) for your reference.

### Template updates

Over time, new versions of nf-core/tools will be released with changes to the template.
In order to keep all _nf-core_ pipelines in sync, we have developed an automated synchronisation procedure.
A GitHub bot account, [@nf-core-bot](https://github.com/nf-core-bot) is scripted on a new tools release to use `nf-core create` with the new template using the input values you used on your pipeline.
This is committed to the `TEMPLATE` branch and a pull-request created to incorporate these changes into `dev`.

Note that these PRs can sometimes create git _merge conflicts_ which will need to be resolved manually.
There are plugins for most code editors to help with this process.
Once resolved and checked, this PR can be merged and a new pipeline release created.

### Exercise 5 (releasing pipelines)

- Use `nf-core bump-version --nextflow` to update the required version of Nextflow in your pipeline
- Bump your pipeline's version to 1.0, ready for its first release!
- Make sure that you're signed up to the _nf-core_ slack (get an invite on [nf-co.re](https://nf-co.re)) and drop us a line about your latest and greatest pipeline plans!
- Ask to be a member of the _nf-core_ GitHub organisation by commenting on [this GitHub issue](https://github.com/nf-core/nf-co.re/issues/3)
- If you're a twitter user, make sure to follow the [@nf_core](https://twitter.com/nf_core) account

## Conclusion

We hope that this _nf-core_ tutorial has been helpful!
Remember that there is more in-depth documentation on many of these topics available on the [nf-core website](https://nf-co.re).
If in doubt, please ask for help on Slack.

If you have any suggestions for how to improve this tutorial, or spot any mistakes, please create an issue or pull request on the [nf-core/nf-co.re repository](https://github.com/nf-core/nf-co.re).

> [Phil Ewels](https://github.com/ewels/), [Maxime Garcia](https://github.com/MaxUlysse), [Gisela Gabernet](https://github.com/ggabernet), [Friederike Hanssen](https://github.com/FriederikeHanssen) for _nf-core_, March 2022
