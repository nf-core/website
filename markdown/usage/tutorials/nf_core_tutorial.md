---
title: nf-core tutorial
subtitle: Tutorial covering the basics of using and creating nf-core pipelines.
---

> Material originally written for the Nextflow Camp 2019, Barcelona 2019-09-19: ***"Getting started with nf-core"*** _(see [programme](https://www.nextflow.io/nfcamp/2019/phil2.html))._
>
> Duration: **1hr 45**
>
> Updated for the nf-core Hackathon 2020, London 2020-03 _(see [event](https://nf-co.re/events#hackathon-francis-crick-2020))._

<!-- markdownlint-disable -->
<iframe src="//www.slideshare.net/slideshow/embed_code/key/sToqg2FJcGUJZ2" width="595" height="485" frameborder="0" marginwidth="0" marginheight="0" scrolling="no" style="border:1px solid #CCC; border-width:1px; margin-bottom:5px; max-width: 100%;" allowfullscreen> </iframe>
<!-- markdownlint-restore -->

**[Click here to download the slides associated with this tutorial.](/assets/markdown_assets/usage/nf_core_tutorial/nf-core_tutorial.pdf)**

#### Exercises

* [1 - Installation](#exercise-1-installation)
* [2 - Listing pipelines](#exercise-2-listing-pipelines)
* [3 - Using pipelines](#exercise-3-using-pipelines)
* [4 - Creating pipelines](#exercise-4-creating-pipelines)
* [5 - Testing pipelines](#exercise-5-testing-pipelines)
* [6 - Releasing pipelines](#exercise-6-releasing-pipelines)

## Abstract

The _nf-core_ community provides a range of tools to help new users get to grips with Nextflow - both by providing complete pipelines that can be used out of the box, and also by helping developers with best practices.
Companion tools can create a bare-bones pipeline from a template scattered with `TODO` pointers and CI with linting tools check code quality.
Guidelines and documentation help to get Nextflow newbies on their feet in no time.
Best of all, the _nf-core_ community is always on hand to help.

In this tutorial we discuss the best-practice guidelines developed by the _nf-core_ community, why they're important and give insight into the best tips and tricks for budding Nextflow pipeline developers. ✨

## Introduction

### What is _nf-core_

_nf-core_ is a community-led project to develop a set of best-practice pipelines built using Nextflow.
Pipelines are governed by a set of guidelines, enforced by community code reviews and automatic linting (code testing).
A suite of helper tools aim to help people run and develop pipelines.

### What this tutorial will cover

This tutorial attempts to give an overview of how _nf-core_ works:

* how to run _nf-core_ pipelines
* how to make new pipelines using the _nf-core_ template
* how _nf-core_ pipelines are reviewed and ultimately released.

### Where to get help

The beauty of _nf-core_ is that there is lots of help on offer!
The main place for this is Slack - an instant messaging service.

The _nf-core_  Slack can be found at [https://nfcore.slack.com](https://nfcore.slack.com/) (_NB: no hyphen in `nfcore`!_).
To join you will need an invite, which you can get at [https://nf-co.re/join/slack](https://nf-co.re/join/slack).

The _nf-core_ Slack organisation has channels dedicated for each pipeline, as well as specific topics (_eg._ [#help](https://nfcore.slack.com/channels/help), [#pipelines](https://nfcore.slack.com/channels/pipelines), [#tools](https://nfcore.slack.com/channels/tools), [#configs](https://nfcore.slack.com/channels/configs) and much more).

One additional tool which we like a lot is [TLDR](https://tldr.sh/) - it gives concise command line reference through example commands for most linux tools, including `nextflow`, `docker`, `singularity`,  `conda`, `git` and more.
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

If using conda, first set up Bioconda as described in the [bioconda docs](https://bioconda.github.io/user/install.html) (especially setting the channel order) and then install nf-core:

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

* Install nf-core/tools
* Use the help flag to list the available commands

## Listing available nf-core pipelines

As you saw from the `--help` output, the tool has a range of sub-commands.
The simplest is `nf-core list`, which lists all available _nf-core_ pipelines.
The output shows the latest version number, when that was released.
If the pipeline has been pulled locally using Nextflow, it tells you when that was and whether you have the latest version.

If you supply additional keywords after the command, the listed pipeline will be filtered.
Note that this searches more than just the displayed output, including keywords and description text.
The `--sort` flag allows you to sort the list (default is by most recently released) and `--json` gives JSON  output for programmatic use.

### Exercise 2 (listing pipelines)

* Use the help flag to print the list command usage
* List all pipelines
* Sort pipelines alphabetically, then by popularity (stars)
* Fetch one of the pipelines using `nextflow pull`
* Use `nf-core list` to see if the pipeline you pulled is up to date
* Filter pipelines for those that work with RNA
* Save these pipeline details to a JSON file

## Running nf-core pipelines

### Software requirements for nf-core pipelines

In order to run _nf-core_ pipelines, you will need to have Nextflow installed ([https://www.nextflow.io](https://www.nextflow.io/)).
The only other requirement is a software packaging tool: [Conda](https://docs.conda.io/en/latest/miniconda.html), [Docker](https://www.docker.com) or [Singularity](https://sylabs.io/singularity/).
In theory it is possible to run the pipelines with software installed by other methods (_e.g._ environment modules, or manual installation), but this is not recommended.
Most people find either Docker or Singularity the best options.

### Fetching pipeline code

Unless you are actively developing pipeline code, we recommend using the Nextflow [built-in functionality](https://www.nextflow.io/docs/latest/sharing.html) to fetch _nf-core_ pipelines.
Nextflow will automatically fetch the pipeline code when you run `nextflow run nf-core/PIPELINE`.
For the best reproducibility, it is good to explicitly reference the pipeline version number that you wish to use with the `-revision`/`-r` flag.
For example:

```bash
nextflow run nf-core/rnaseq -revision 1.3
```

If not specified, Nextflow will fetch the default branch.
For pipelines with a stable release this the default branch is `master` - this branch contains code from the latest release.
For pipelines in early development that don't have any releases, the default branch is `dev`.

If you would like to run the latest development code, use `-r dev`.

Note that once pulled, Nextflow will use the local cached version for subsequent runs.
Use the `-latest` flag when running the pipeline to always fetch the latest version.
Alternatively, you can force Nextflow to pull a pipeline again using the `nextflow pull` command:

```bash
nextflow pull nf-core/rnaseq
```

### Usage instructions and documentation

You can find general documentation and instructions for Nextflow and _nf-core_ on the _nf-core_ website: [https://nf-co.re/](https://nf-co.re/).
Pipeline-specific documentation is bundled with each pipeline in the `/docs` folder.
This can be read either locally, on GitHub, or on the _nf-core_ website.
Each pipeline has its own webpage at `https://nf-co.re/<pipeline_name>` (_e.g._ [nf-co.re/rnaseq](https://nf-co.re/rnaseq))

In addition to this documentation, each pipeline comes with basic command line reference.
This can be seen by running the pipeline with the `--help` flag, for example:

```bash
nextflow run nf-core/rnaseq --help
```

### Config profiles

Nextflow can load pipeline configurations from multiple locations.
To make it easy to apply a group of options on the command line, Nextflow uses the concept of [config profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles).
 _nf-core_ pipelines load configuration in the following order:

1. Pipeline: Default 'base' config
    * Always loaded. Contains pipeline-specific parameters and "sensible defaults" for things like computational requirements
    * Does _not_ specify any method for software packaging. If nothing else is specified, Nextflow will expect all software to be available on the command line.
2. Pipeline: Core config profiles
    * All _nf-core_ pipelines come with some generic config profiles. The most commonly used ones are for software packaging: `docker`, `singularity` and `conda`
    * Other core profiles are `debug` and `test`
3. [nf-core/configs](https://github.com/nf-core/configs): Server profiles
    * At run time, _nf-core_ pipelines fetch configuration profiles from the [configs](https://github.com/nf-core/configs) remote repository. The profiles here are specific to clusters at different institutions.
    * Because this is loaded at run time, anyone can add a profile here for their system and it will be immediately available for all _nf-core_ pipelines.
4. Local config files given to Nextflow with the `-c` flag
5. Command line configuration

Multiple comma-separate config profiles can be specified in one go, so the following commands are perfectly valid:

```bash
nextflow run nf-core/rnaseq -profile test,docker
nextflow run nf-core/rnaseq -profile singularity,debug
```

Note that the order in which config profiles are specified matters.
Their priority increases from left to right.

Our tip: Be clever with multiple Nextflow configuration locations. For example, use `-profile` for your cluster configuration, `~/.nextflow/config` for your personal config such as `params.email` and a working directory `nextflow.config` file for reproducible run-specific configuration.

### Running pipelines with test data

The `test` config profile is a bit of a special case.
Whereas all other config profiles tell Nextflow how to run on different computational systems, the `test` profile configures each `nf-core` pipeline to run without _any_ other command line flags.
It specifies URLs for test data and all required parameters.
Because of this, you can test any _nf-core_ pipeline with the following command:

```bash
nextflow run nf-core/<pipeline_name> -profile test
```

Note that you will typically still need to combine this with a configuration profile for your system - _e.g._ `-profile test,docker`.
Running with the test profile is a great way to confirm that you have Nextflow configured properly for your system before attempting to run with real data.

### The nf-core launch command

Most _nf-core_ pipelines have a number of flags that need to be passed on the command line: some mandatory, some optional.
To make it easier to launch pipelines, these parameters are described in a JSON file bundled with the pipeline.
The `nf-core launch` command uses this to build an interactive command-line wizard which walks through the different options with descriptions of each, showing the default value and prompting for values.

> _NB: This is an experimental feature - JSON file and rich descriptions of parameters is not yet available for all pipelines._

Once all prompts have been answered, non-default values are saved to a `params.json` file which can be supplied to Nextflow to run the pipeline. Optionally, the Nextflow command can be launched there and then.

To use the launch feature, just specify the pipeline name:

```bash
nf-core launch <pipeline_name>
```

### Using nf-core pipelines offline

Many of the techniques and resources described above require an active internet connection at run time - pipeline files, configuration profiles and software containers are all dynamically fetched when the pipeline is launched.
This can be a problem for people using secure computing resources that do not have connections to the internet.

To help with this, the `nf-core download` command automates the fetching of required files for running _nf-core_  pipelines offline.
The command can download a specific release of a pipeline with `-r`/`--release` and fetch the singularity container if `--singularity` is passed (this needs Singularity to be installed).
All files are saved to a single directory, ready to be transferred to the cluster where the pipeline will be executed.

### Exercise 3 (using pipelines)

* Install required dependencies (Nextflow, Docker)
* Print the command-line usage instructions for the _nf-core/rnaseq_ pipeline
* In a new directory, run the _nf-core/rnaseq_ pipeline with the provided test data
* Try launching the RNA pipeline using the `nf-core launch` command
* Download the _nf-core/rnaseq_ pipeline for offline use using the `nf-core download` command

## Creating nf-core pipelines

### Using the nf-core template

The heart of _nf-core_ is the standardisation of pipeline code structure.
To achieve this, all pipelines adhere to a generalised pipeline template.
The best way to build an _nf-core_ pipeline is to start by using this template via the `nf-core create` command.
This launches an interactive prompt on the command line which asks for things such as pipeline name, a short description and the author's name.
These values are then propagated throughout the template files automatically.

### TODO statements

Not everything can be completed with a template and all new pipelines will need to edit and add to the resulting pipeline files in a similar set of locations.
To make it easier to find these, the _nf-core_ template files have numerous comment lines beginning with `TODO nf-core:`, followed by a description of what should be changed or added.
These comment lines can be deleted once the required change has been made.

Most code editors have tools to automatically discover such `TODO` lines and the `nf-core lint` command will flag these.
This makes it simple to systematically work through the new pipeline, editing all files where required.

### How nf-core software packaging works

The only hard requirement for all _nf-core_ pipelines is that software must be available in Docker images. However, it is recommended that pipelines use the following methodology where possible:

1. Software requirements are defined for Conda in `environment.yml`
2. Docker images are automatically built on Docker Hub, using Conda
3. Singularity images are generated from Docker Hub at run time for end users

This approach has the following merits:

* A single file contains a list of all required software, making it easy to maintain
* Identical (or as close as is possible) software is available for users using Conda, Docker or Singularity
* Having a single container image for the pipeline uses disk space efficiently for Singularity images, and is simple to manage and transfer.

The reason that the above approach is not a hard requirement is that some issues can prevent it from working, such as:

* It may not be possible to package software on conda due to software licensing limitations
* Different packages may have dependency conflicts which are impossible to resolve

Alternative approaches are then decided upon on a case-by-case basis.
We encourage you to discuss this on Slack early on as we have been able to resolve some such issues in the past.

### Building environment.yml

The _nf-core_ template will create a simple `environment.yml` file for you with an environment name, conda channels and one or two dependencies.
You can then add additional required software to this file.
Note that all software packages must have a specific version number pinned - the format is a single equals sign, _e.g_ `package=version`.

Where software packages are not already available on Bioconda or Conda-forge, we encourage developers to add them.
This benefits the wider community, as well as just users of the _nf-core_  pipeline.

### Running with Docker locally

You can use Docker for testing by building the image locally.
The pipeline expects a container with a specific name, so you must _tag_ the Docker image with this.
You can build and tag an image in a single step with the following command:

```bash
docker build -t nfcore/PIPELINE:dev .
```

Note that it is `nfcore` without a hyphen (Docker Hub doesn't allow any punctuation).
The `.` refers to the current working directory - if run in the root pipeline folder this will tell Docker to use the `Dockerfile` recipe found there.

### Forks, branches and pull-requests

All _nf-core_ pipelines use GitHub as their code repository, and git as their version control system.
For newcomers to this world, it is helpful to know some of the basic terminology used:

* A _repository_ contains everything for a given project
* _Commits_ are code checkpoints.
* A _branch_ is a linear string of commits - multiple parallel branches can be created in a repository
* Commits from one branch can be _merged_ into another
* Repositories can be _forked_ from one GitHub user to another
* Branches from different forks can be merged via a _Pull Request_ (PR) on [github.com](https://github.com)

Typically, people will start developing a new pipeline under their own personal account on GitHub.
When it is ready for its first release and has been discussed on Slack, this repository is forked to the _nf-core_ organisation.
All developers then maintain their own forks of this repository, contributing new code back to the _nf-core_ fork via pull requests.

All _nf-core_ pipelines must have the following three branches:

1. `master` - commits from stable releases _only_. Should always have code from the most recent release.
2. `dev` - current development code. Merged into `master` for releases.
3. `TEMPLATE` - used for template automation by the [@nf-core-bot](https://github.com/nf-core-bot) GitHub account. Should only contain commits with unmodified template code.

Pull requests to the _nf-core_ fork have a number of automated steps that must pass before the PR can be merged.
A few points to remember are:

* The pipeline `CHANGELOG.md` must be updated
* PRs must not be against the `master` branch (typically you want `dev`)
* PRs should be reviewed by someone else before being merged (help can be asked on the [#request-review](https://nfcore.slack.com/channels/request-review) Slack channel)

### Setting up Docker Hub

When you fork your pipeline repository to the _nf-core_ organisation, one of the core team will set up Docker Hub (automated Docker image creation) for you.
However, it can be helpful to set it up on your personal fork as well.
That way, you can be confident that everything will work when you fork or open a PR on the _nf-core_ organisation.

This service is free to use.
To set it up, visit [https://hub.docker.com](https://hub.docker.com) and link your personal GitHub repository.

To set up an automated docker container build for the master and the dev branches, as well as for all releases, go to the builds tab of your docker hub repository. There, select build configurations. The source repository should be set to your GitHub repository. Three build rules should be set:

* Master branch build rule (comes by default). Set the "Source Type" to "branch", "Source" to "master" and "Docker tag" to "latest".
* dev branch build rule. Set _Source Type_ to `branch`, _Source_ to `dev` and _Docker tag_ to `dev`.
* Tag build rule: this will build a container for your pipeline releases, and is similar to the scenario _"Match versions"_. Set the _Source Type_ to `Tag`, _Source_ to `/^[0-9.]+$/` and _Docker Tag_ to `{sourceref}`.
* For all build rules, the _Dockerfile location_ should be set to `Dockerfile` and the _Build context_ set to `/`. _Autobuild_ and _Build Caching_ should be on.

### Exercise 4 (creating pipelines)

* Make a new pipeline using the template
* Update the readme file to fill in the `TODO` statements
* Add a new process to the pipeline in `main.nf`
* Add the new software dependencies from this process in to `environment.yml`

## Testing nf-core pipelines

### Linting nf-core pipelines

Manually checking that a pipeline adheres to all _nf-core_ guidelines and requirements is a difficult job.
Wherever possible, we automate such code checks with a [code linter](https://en.wikipedia.org/wiki/Lint_(software)).
This runs through a series of tests and reports failures, warnings and passed tests.

The linting code is closely tied to the _nf-core_ template and both change over time.
When we change something in the template, we often add a test to the linter to make sure that pipelines do not use the old method.

Each lint test has a number and is documented on the [nf-core website](https://nf-co.re/tools-docs).
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

To set up the test profile, make a new branch on the  [nf-core/test-datasets](https://github.com/nf-core/test-datasets) repo through the web page (see [instructions](https://help.github.com/en/articles/creating-and-deleting-branches-within-your-repository#creating-a-branch)).
Fork the repository to your user and open a PR to your new branch with a really (really!) tiny dataset.
Once merged, set up the `conf/test.config` file in your pipeline to refer to the URLs for your test data.

These test datasets are used by the automated continuous integration tests.
The systems that run these tests are _extremely_ limited in the resources that they have available.
Typically, the pipeline should be able to complete in around 10 minutes and use no more than 6-7 GB memory.
To achieve this, input files and reference genomes need to be very tiny.
If possible, a good approach can be to use PhiX or Yeast as a reference genome.
Alternatively, a single small chromosome (or part of a chromosome) can be used.
If you are struggling to get the tests to run, ask for help on Slack.

When writing `conf/test.config` remember to define all required parameters so that the pipeline will run with only `-profile test`.
Note that remote URLs cannot be traversed like a regular file system - so glob file expansions such as `*.fa` will not work.

### GitHub Actions configuration

The automated tests with GitHub Actions are configured in the `.github/worfklows` folder that is generated by the template.
Each file Each file (`branch.yml`, `ci.yml` and `linting.yml`) defines several tests: branch protection, running the pipeline with the test data, linting the code with `nf-core lint`, linting the syntax of all Markdown documentation and linting the syntax of all YAML files.

#### `branch.yml`

This is already set up and does not need to be edited.
It will check that pull-requests going to the nf-core repo `master` branch only come from the `dev` or `patch` branches (for releases).
In case you want to add branch protection for a repository out of _nf-core_, you can add an extra step to the workflow with that repository name.
You can leave the _nf-core_ branch protection in place so that the `nf-core lint` command does not throw an error - it only runs on nf-core repositories so it should be ignored.

#### `linting.yml`

This is already set up and does not need to be edited.
It will check:

* That all Markdown documentation follow a proper syntax.
  * Many code editors have similar packages to help with markdown validation. For example, [markdownlint for atom](https://atom.io/packages/linter-markdownlint) and [markdownlint for Visual Studio Code](https://marketplace.visualstudio.com/items?itemName=DavidAnson.vscode-markdownlint))
  * You can run the markdown linting on the command line by [installing markdownlint](https://www.npmjs.com/package/markdownlint) and running the command `markdownlint . -c .github/markdownlint.yml`
* That all YAML file follow a proper syntax.
* That the pipeline code follows the nf-core linting rules.
  * You can run this manually with the `nf-core lint` command.

#### `ci.yml`

This file defines a GitHub Actions workflow that tests your pipeline with the small test dataset.
It should mostly work out of the box, but it may need some editing from you to customise it for your pipeline.

The matrix `nxf_ver` variable sets the `NXF_VER` environment variable twice.
This tells GitHub Actions to run the tests twice in parallel - once with the latest version of Nextflow (`NXF_VER=''`) and once with the minimum version supported by the pipeline.
Do not edit this version number manually - it appears in multiple locations through the pipeline code, so it's better to use `nf-core bump-version --nextflow` instead.

The provided tests run your pipeline with the `-profile test,docker` flags. This may be sufficient for your pipeline.
However, if it is possible to run the pipeline with significantly different options (for example, different alignment tools), then it is good to test all of these.
You can do this by adding additional tests in the `jobs` block.
Do not try add a run for every possible combination of parameters, as it will take too long to run.

### Exercise 5 (testing pipelines)

* Run `nf-core lint` on your pipeline and make note of any test warnings / failures
* Read up on one or two of the linting rules on the nf-core website and see if you can fix some.
* Take a look at `conf/test.config` and switch the test data for another dataset on nf-core/test_data.

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

Unfortunately, at the time of writing, Docker Hub and Zenodo (automated DOI assignment for releases) services are not created automatically.
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

### Template updates

Over time, new versions of nf-core/tools will be released with changes to the template.
In order to keep all _nf-core_ pipelines in sync, we have developed an automated synchronisation procedure.
A GitHub bot account, [@nf-core-bot](https://github.com/nf-core-bot) is scripted on a new tools release to use `nf-core create` with the new template using the input values you used on your pipeline.
This is committed to the `TEMPLATE` branch and a pull-request created to incorporate these changes into `dev`.

Note that these PRs can sometimes create git _merge conflicts_ which will need to be resolved manually.
There are plugins for most code editors to help with this process.
 Once resolved and checked this PR can be merged and a new pipeline release created.

### Exercise 6 (releasing pipelines)

* Use `nf-core bump-version --nextflow` to update the required version of Nextflow in your pipeline
* Bump your pipeline's version to 1.0, ready for its first release!
* Make sure that you're signed up to the _nf-core_ slack (get an invite on [nf-co.re](https://nf-co.re)) and drop us a line about your latest and greatest pipeline plans!
* Ask to be a member of the _nf-core_ GitHub organisation by commenting on [this GitHub issue](https://github.com/nf-core/nf-co.re/issues/3)
* If you're a twitter user, make sure to follow the [@nf_core](https://twitter.com/nf_core) account

## Conclusion

We hope that this _nf-core_ tutorial has been helpful!
Remember that there is more in-depth documentation on many of these topics available on the [nf-core website](https://nf-co.re).
If in doubt, please ask for help on Slack.

If you have any suggestions for how to improve this tutorial, or spot any mistakes, please create an issue or pull request on the [nf-core/nf-co.re repository](https://github.com/nf-core/nf-co.re).

> [Phil Ewels](https://github.com/ewels/), [Maxime Garcia](https://github.com/MaxUlysse) for _nf-core_, February 2020
