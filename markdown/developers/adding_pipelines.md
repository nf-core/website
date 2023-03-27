---
title: Adding a new pipeline
subtitle: Follow this walkthrough to add a new pipeline to nf-core.
---

## Before you start

So, you want to add a new pipeline to nf-core - brilliant!
Before you start typing, check that you're happy with the following points:

- You're familiar with nf-core and nextflow (see our [introduction docs](/docs/usage/introduction.md)).
- You're used to working with `git` and [GitHub](https://github.com)
  (see a [nice tutorial here](https://blog.scottlowe.org/2015/01/27/using-fork-branch-git-workflow/))
- The workflow you're thinking of meets the [nf-core guidelines](https://nf-co.re/docs/contributing/guidelines).

The main steps involved in adding a new nf-core pipeline covered below are:

1. [Joining the community](#join-the-community)
2. [Creating a pipeline](#create-a-pipeline-from-the-template)
3. [Adding test data](#add-some-test-data)
4. [Adding to the nf-core organisation](#adding-your-pipeline-to-the-nf-core-organisation)
5. [Making your first release](#making-the-first-release)
6. [Updates and new releases](#subsequent-releases)

## Join the community

At its heart, nf-core is a community - to add a pipeline you need to be part of that community!
Please request to join the [nf-core GitHub organisation](https://github.com/nf-core/nf-co.re/issues/3))
and introduce yourself on [Slack](https://nf-co.re/join/slack) or the
[mailing list](https://groups.google.com/forum/#!forum/nf-core).

**⚠️ It's good to introduce your idea early on so that it can be discussed, before you spend lots of time coding. ⚠️**

The [nf-core guidelines](/docs/contributing/guidelines) state that no two pipelines should overlap too much
in their purpose and results. There may be an existing pipeline that can be extended to give the
functionality that you are looking for, or there could be another group working on a similar to the
pipeline to the one you're planning.

To avoid problems at a later date, please come and discuss your plans with the nf-core community as early
as possible. Ideally before you make a start on your pipeline!

> Not all pipelines are suitable for inclusion in the main nf-core community (eg. bespoke or proprietary workflows). However, we hope that you may still wish to use the nf-core template and/or use components of nf-core code. All nf-core code is under a MIT license and where possible we have endeavoured to make the tools work with any Nextflow pipeline. If this is the case for you, please see the [unofficial pipelines tutorial](/docs/contributing/tutorials/unofficial_pipelines.md) for more details.

All nf-core discussion happens on the nf-core Slack, which you can join here:
[https://nf-co.re/join](https://nf-co.re/join)

These topics are specifically discussed in the `#new-pipelines` channel:
[https://nfcore.slack.com/channels/new-pipelines](https://nfcore.slack.com/channels/new-pipelines)

## Create a pipeline from the template

You'll start by making a new pipeline locally and working with it on your own GitHub account.
Only when it's ready do we move it to the nf-core GitHub organisation.

It's _highly_ recommended to use the nf-core template.
The guidelines for nf-core pipelines are pretty strict, but if you start your pipeline by using the
nf-core template (`nf-core create` - see [the docs](https://nf-co.re/tools#creating-a-new-pipeline))
then your life will be much easier.
This tool does lots of things for you: it gives you the correct file structure and boiler plate code
and also sets up the required `git` infrastructure for you to keep your pipeline in sync in the future.

Even if you already have a working pipeline, it may be easier in the long run to use this this template
and copy over your code in the relevant places.

If you really don't want to use the template it should possible to work without it.
Please see the [manual synchronisation](sync.md) documentation.

> Note that workflow names should be all lower-case and contain no punctuation.
> This is to allow consistent names between platforms.

### Push to GitHub

Create an empty repository on GitHub for your new pipeline under your personal account.
Do this by going to the GitHub website and clicking + then _New Repository_.
Make sure _not_ to initialise it with _any_ file, `README` or `LICENSE`: you just want an empty repository. It will be populated after you push the template you created locally with all you need to start working on your pipeline.

Once created, copy the URL and add this as a remote to your local git repository
and push your code:

```bash
## Add a remote called 'origin' - this is the default name for a primary remote
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPOSITORY.git
## Commit any new code changes since you ran the template
git add .
git commit -m "Starting to build my pipeline"
## Push to GitHub
git push --all origin
```

### Work on your pipeline

Ok, now you're all set with your own personal nf-core pipeline!
You can now start writing code for real.
Remember to run the `nf-core lint` command (see [docs](https://nf-co.re/tools#linting-a-workflow))
to make sure that your workflow passes all of the nf-core compatibility tests.
The automated tests on Github Actions also run this, so you should get a
notification from GitHub if something breaks.

## Running with test data

Whilst the linting tests are good, they're not sufficient by themselves.
It's also good actually run your pipeline on a minimal dataset.
We also automatically run tests with GitHub Actions anytime someone updates the pipeline code (see below).
Currently, we don't usually check the results that are produced, but it often catches
syntax errors and other serious problems that cause nextflow to exit with an error.

### Putting the test data on GitHub

To avoid bloating the workflow, we don't keep test data in the same repository as
nf-core workflows.
Instead, we use the dedicated [nf-core/test-datasets](https://github.com/nf-core/test-datasets/) repository.

To set this up, make a fork of that repo to your personal account.
Clone the repository and check out a new branch for your workflow:

```console
git clone https://github.com/YOUR_USERNAME/test-datasets.git
cd test-datasets
git checkout -b MY_WORKFLOW
```

Now add your test data files - note that they must be **very small**.
GitHub has quite a low file size limit, and the GitHub Actions will time out with anything
that's not tiny. We typically use PhiX / Yeast / part of a chromosome as a reference
and aggressively subsampled input data. I.e., as small as possible, as large as necessary.
We highly recommend that you ask in the [#test-data channel](https://nfcore.slack.com/channels/test-data) on the [nf-core slack](https://nf-co.re/join) for advice
before adding any test data!

Once added, push these new files to GitHub:

```console
git add .
git commit -m "Added test data for MY_WORKFLOW"
git push --set-upstream origin MY_WORKFLOW
```

Finally, make a pull-request against the main [nf-core/test-datasets](https://github.com/nf-core/test-datasets/) repository with your files.
You want this repo to also use a branch with the name of your workflow, so first go
to the [repository GitHub web page](https://github.com/nf-core/test-datasets/) and create
this new branch using the UI there.
Once created, you can open a pull request and select this as the target branch.

If in doubt, ask for help!
([Slack](https://nf-co.re/join/slack) or [mailing list](https://groups.google.com/forum/#!forum/nf-core))

### Setting up a test workflow

Now that your test data is hosted on the web, you can set up a `test` config profile in your
workflow that points to it.
In fact, the `test` profile should already exist if you've used the template.
Switch out the example URLs for the ones you added (view the files on GitHub and click 'Raw' to get the URL).

Add any other required parameters so that running the pipeline runs with as few extra
flags as possible. Note that the `test` profile can be combined with other profiles such as `docker`
or `conda`, so your config should not specify a hardware environment.

Have a go at running the pipeline and see if it works:

```console
nextflow run MY_WORKFLOW -profile test,docker --outdir <OUTDIR>
```

Note that if you do need to adjust this `nextflow run` command, you'll need to update it
in the `.github/workflows/` YAML files too.

## Adding your pipeline to the nf-core organisation

Ok, so you're essentially finished. Your pipeline is written, the tests pass and
you're ready to add your workflow to nf-core.

First, go to the settings of your repository. Under the General page, in the 'Danger Zone' you should have an option to Transfer Ownership. Transfer this to the nf-core organisation.

> You must make sure you are already a part of the nf-core organisation to allow transferring to nf-core. Alternatively you can add a core-team member to your repository, and ask them to do the transfer you.

Once transferred, go to the transferred repository on nf-core, and make a new fork back to your user name or organisation to continue development on a fork.

> We [prefer](https://nfcore.slack.com/archives/CQY2U5QU9/p1676366247726189?thread_ts=1676360232.837109&cid=CQY2U5QU9) transferring over _forking_ to nf-core. If we fork the original repository to nf-core whenever anyone opens a pull-request, it defaults to going to the original user's fork of the repository, not the nf-core repo. In this case the only way to fix to request manual detachment from GitHub support.

### Branch setup

All nf-core pipelines use branches called `dev` and `master`.
The `master` branch should contain the code from the latest stable release, `dev` should have the latest development code.
We want people to run the latest development code by default up until the first release.
To do this, we set `dev` as the default repository branch.
After an initial release is created, we set the default branch back to `master` so that the default
action is to run the latest stable release code.

Once you have forked the repository, create a new branch called `dev` for the active development.
In the repository settings, set `dev` to be the default branch.

### Repository setup

Remember to configure the repository on the GitHub website with the following:

- A description, the [https://nf-co.re](https://nf-co.re) URL and lots of keywords!
- Issues enabled, disable Wiki and Projects
- A protected `master` branch that requires review and passing tests
- Write permissions for `nf-core/all` and admin permissions for `nf-core/admin`

You can check that all of these settings are done correctly by referring to your pipeline
in the nf-core [Repository health web page](https://nf-co.re/pipeline_health).
This reports the status of various checks and also has the option of fixing errors for you via the GitHub API.

### Differences to your own fork

The main difference when working with the main nf-core fork of your workflow is
that tests for pull-requests against the `master` branch will fail. This is because
the `master` branch should only ever contain code from the last release.
Instead, use the `dev` branch for new work and always make pull-requests against
that. Then the tests should pass.

## Making the first release

### Reset the default branch

When the code is stable and ready for a release, set the `master` branch to be the default
branch again.

### Bump the version

At this point you should bump the version numbers on `dev`.

When developing the pipeline, the version numbers should be numeric with `dev` at the end.
Use the `nf-core bump-version` command to do this - there are quite a few locations in the
code that need updating and this ensures that it happens in the correct places.

When making a release, version numbers should all be numeric. Use `nf-core lint --release`
when ready - this will check that everything looks correct. Pipeline release numbers MUST
use [Semantic Versioning](https://semver.org/).

### Core pipeline review

Ok - now the tough bit - does your workflow stand up to the scrutiny of the nf-core team?!
Not to worry, we're a friendly bunch, just let us know about the new pipeline, when you're
ready, following the process below.

Make a pull-request from the `dev` branch to `master` on the nf-core fork. This is a
special case and the tests should pass, and once they do you can request a review from the
core team.

What happens next depends on the state of your master branch:

- If you have developed in such a way that your master branch is clean, .i.e. doesn't have
  any commits since the inital one, the PR created above will represent all changes
  associated with the proposed release, and the core team will use it for review and
  feedback.
- If your master branch already contains changes associated with the release, the core
  team may merge your PR and create a pseudo-PR against the first commit in the
  pipeline. This gives the PR review interface showing all code that you've written.

In either case we will go through everything and request changes that we think are
necessary until you're good to go. Common things that are flagged at this point are:

- A clear, short but descriptive readme
- Good documentation, especially describing the output files and all parameters
- Pipeline code

We typically tend to have two reviewers for most of the crucial code changes, e.g. adding
new major features to an existing pipeline or making an entirely new pipeline release. You
can also ping people from the nf-core core or maintainers team to review your pipeline
code by `@`ing them.

### Making the release

Once any requested changes have been made and the associated PR approved, you can go ahead
with releasing the pipeline. Put in a basic changelog entry describing the general
functionality at release. When you're ready, follow the instructions in the nf-core
[release checklist](release_checklist.md). We recommend you also explicitly tag contributors with their GitHub handles, so each release on GitHub will display their icons.

The nf-core website and helper tools will automatically detect the new release and be updated accordingly.

That's it, you're finished! Congratulations!

### Subsequent releases

Once you've made your first release you can continue to work on your fork and make pull-requests
against the `dev` branch on the nf-core repository. Now that we have a stable `master` branch,
there should be reviews of each PR against `dev` before merging.

When ready to make new releases, make sure that the version number is increased and create a
pull-request against `master`. If tests pass, it can be merged and a new release made.

The `master` branch should always have only the commit from the latest release. This is important
because the commit ID is used to reference whether the pipeline is up to date or not.

### Adding new pipeline features to existing pipelines

We are an open and inclusive community, welcoming any contributions to pipelines already present in nf-core. In many cases, the original developers might either not have experience with some new fancy method or simply doesn't have the time to implement everything themselves - so they might be really happy to see you actively contributing!

Basic rules for such contributions:

- Ask in the [Slack](https://nf-co.re/join/slack) channel for the specific pipeline whether there is an open issue on the respective pipeline's issue tracker for the feature you're planning to
- If not, create a new issue there, describing the purpose and ideas you have and wait for someone to comment/discuss
- If everyone is happy or there is some consensus in the community, start implementing the feature in your [fork](https://help.github.com/en/articles/fork-a-repo) of the respective pipeline
- Please do not write to multiple channels in the Slack community, rather collect all of the information in a single GitHub issue, which makes it also much easier to follow up on your proposal

### Adding new dependencies to an existing pipeline

Sometimes, especially when adding new features to a pipeline, the dependencies change as well. In such cases, you might want to have an updated Docker Container available before submitting a pull request, in order to have the GitHub Actions tests run through when testing your updated code. To achieve that, please follow these steps:

- Add _only_ the newly required dependencies to the `environment.yml` in the pipeline code
- List this new dependency as something new in the `CHANGELOG`
- Create a Pull Request including only these two changes against the `dev` branch of the pipeline you're working on

This way, a review process will be very fast and we can merge the changes into the `dev` branch, updating the Docker Image for that pipeline automatically. After ~30 Minutes, the Docker Image for that pipeline is then updated, and you can open your Pull Request containing your actual pipeline code changes.

## nf-core pipeline structure

### Files

You will find the following files in each nf-core pipeline. They are automatically generated, when running `nf-core create`.

- `main.nf`: This is the main nextflow file which will get executed if the pipeline is run. Typically, parameters are initialized and validated in this script before a workflow from the `workflow/` directory is called for execution.

* `nextflow.config`: The main nextflow configuration file. It contains the default pipeline parameters, nextflow configuration options and information like pipeline and minimum nextflow version, among others.
  The `nextflow.config` also defines different configuration profiles that can be used to run the pipeline. See the [Configuration docs](/docs/usage/configuration) for more information.

- `README.md`: Basic information about the pipeline and usage

- `nextflow_json.schema`: The JSON schema file is used for pipeline parameter specification. This is automatically created using the `nf-core schema build` command. It is used for printing command-line help, validating input parameters, building the website docs and for building pipeline launch interfaces (web and cli).

- `CHANGELOG.md`: Information about the changes made to the pipeline for each release.

- `LICENSE`: The license - should be MIT

- `CODE_OF_CONDUCT.md`: The nf-core code of conduct.

- `CITATIONS.md`: All citations needed when using the pipeline

- `.gitattributes`: Git settings, primarily getting the `.config` files to render with Nextflow syntax highlighting on <github.com>

- `.gitignore`: Files that should be ignored by git.

- `.editorconfig`: Editorconfig file that helps assuring consistent coding style

- `.markdownlint.yml`: Markdown lint configuration file to assure consistent markdown files

- `modules.json`: This file holds information (e.g. version) about all the modules in the pipeline that have been installed from `nf-core/modules`

### Directories

- `.github/`: Other GitHub specific files, e.g. for specifying templates and GitHub actions

- `assets/`: Any additional files needed for the pipeline

- `bin/`: Directory for scripts that must be directly accessible within a pipeline process. Anything in this directory can be directly called from within Nextflow processes.

- `conf/`: Configuration files, including a `base.config` file which is always loaded into `nextflow.config` and describes basic pipeline configurations, like CPU and memory usage for processes with low, medium and high requirements. Additionally, most pipelines also have a `igenomes.config` file which describes the locations of popular genomes that can be automatically downloaded for a pipeline run. Finally, the `test.config` and `test_full.config` files are test configurations that are loaded during test runs. Since DSL2, it also contains a `modules.config` file, which defines module-specific configurations and is explained further down in the "DSL2 and modules" section.

- `docs/`: Markdown files for documenting the pipeline

- `lib/`: The lib directory contains Groovy utility functions. These are called from within the nf-core pipeline to do common pipeline tasks (e.g. parameter schema validation) and to hold Groovy functions that may be useful in the pipeline context (e.g. to validate pipeline-specific parameters). Currently, the following files are included:

  - `NfcoreSchema.groovy` - Functions to validate input parameters using the pipeline JSON schema
  - `NfcoreTemplate.groovy` - Additional nf-core specific pipeline functions (sending emails, checking nf-core config profiles and more)
  - `Utils.groovy` - Additional generic pipeline functions (checking conda config and more)
  - `WorkflowMain.groovy` - Startup functions for the main pipeline (printing logs, custom params initialisation and more)
  - `WorkflowPipeline.groovy` - Functions for pipeline subworkflows
  - `nfcore_external_java_deps.jar` - Bundled Groovy dependencies so that pipelines work offline (mostly for JSON schema validation - see imports in `NfcoreSchema.groovy`)

- `modules/`: Contains pipeline-specific and common nf-core modules

- `workflows/`: Contains the main pipeline workflows to be executed in the `main.nf` file

- `subworkflows/`: Contains smaller subworkflows that typically consist out of a few modules chained together

## Continuous integration testing

To assure that nf-core pipelines don't break after some change is made to the code, we use automated continuous integration (CI) testing. This is done via GitHub actions, which are defined in the `.github/workflows` directory. Parameters and file paths are set in the `conf/test.config` and `conf/test_full.config`. Please see also [here](/docs/contributing/adding_pipelines#add-some-test-data) for how to set-up the test workflow for your pipeline.

## DSL2 and modules

Nextflow DSL2 allows for a more modularized design of pipelines and the reuse of components. Currently, most nf-core pipelines are still entirely written in DSL1, but in the near future all pipelines will be written in DSL2. The nf-core team has developed a set of design patterns on how to best implement DSL2 pipelines, which should be used by all nf-core pipelines in order to assure standardization and the reuse of components. The following is meant to help understand certain design choices and how a nf-core DSL2 pipeline should be build.

### Modules

Each pipeline has a `modules` directory which contains all the module code. A module here depicts a single process which involves - if possible - only a single tool/software. The `modules` directory is furthermore divided into `local`and `nf-core` sub-directories, which themselves each have a `process`/`software` and `subworkflow` directory. Modules contained in the `local` directory are specific to the pipeline, whereas `nf-core` modules are installed from the `nf-core/modules` repository. For instance, most pipelines that involve FastQ files will run the FastQC tool for quality control. The module needed for this can be easily reused from the `nf-core/modules` directory using the `nf-core/tools`package. The `process` directories contain modules which define single processes, which smaller workflows are contained in the `subworkflow` directories.

All modules load utility functions from a `functions.nf` script that must be contained in the `modules/local/process` directory. It contains simple functions to initialize the module options, get the software version, save files and get a path from a string. For further explanations of modules and how they should be structured in DSL2 pipelines, check out the [nf-core modules repo](https://github.com/nf-core/modules).

### Module parameters

One thing that might not be straightforward is how module parameters are handled in nf-core DSL2 pipelines. Every module and subworkflow, when loaded into the pipeline, has to be passed a groovy map containing module options. For single processes this is typically only a single `options` map, while for subworkflows these can be several maps that are then passed down to the correct processes within the subworkflows. These `options` maps are directly loaded from the `modules.config` file (contained in the pipeline `conf` directory), which is the place where all additional and optional parameters for modules are stored. Modules should be build in a way such that they are flexible with respect to the parameters, so that most command line parameters can be passed to them via the `modules.config`. This way, all command line parameters and other options can be modified within a single script, which makes it easy for users to adjust the pipeline and at the same time makes modules more reusable.

The `modules.config` file should contain a `params.modules` dictionary which lists every module used in the pipeline. For each module, the following fields can be specified:

- `args`: additional arguments appended to command in the module
- `args2`: Second set of arguments append to command in the module (multi-tool modules)
- `publish_dir`: Directory to publish the results
- `publish_by_id`: Publish results in separate folder by meta.id value
- `publish_files`: Groovy map where key = "file_ext" and value = "directory" to publish results for that file extension. The value of "directory" is appended to the standard "publish_dir" path as defined above. If publish_files == null (unspecified) all files are published. If publish_files == false no files are published.
- `suffix`: File name suffix for output files

### Sample meta information

In nf-core DSL2 pipelines, every channel that contains sample data in some way should also contain a `meta`variable, which must contain the fields `meta.id`, `meta.single_end` and `meta.strandedness`. The `meta` variable can be passed down to processes as a tuple of the channel containing the actual samples, e.g. FastQ files, and the `meta` variable. This meta information can easily be extracted from a samplesheet which specifies the input files. For an example process that reads a samplesheet, creates a `meta` variable and returns it along with the filepaths, have a look at the [rnaseq pipeline](https://github.com/nf-core/rnaseq/blob/master/modules/local/samplesheet_check.nf).
