---
title: Creating pipelines with nf-core
subtitle: Learn how to use nf-core to build your pipeline
menu:
  main:
    weight: 20
---

> Material originally written for the [Nextflow & nf-core training, October 2022](https://nf-co.re/events/2022/training-october-2022).

## Introduction

Welcome to the nf-core section of the Nextflow & nf-core training course!
You can follow the commands and examples covered within this training in the document below.

### Using GitPod

As with the rest of the Nextflow training up until this point, we will be using GitPod for this training.
We will use a different GitPod environment in order to get the very latest releases of the nf-core tools.
To launch GitPod, follow the link:

<p class="text-center">
  <a href="https://www.gitpod.io/#https://github.com/nf-core/tools" class="btn btn-lg btn-success" target="_blank">
    Launch GitPod
  </a>
</p>

#### Working in an empty directory

To work with a clean directory, you can do the following:

- In the terminal, make a new directory to work in:
  ```bash
  cd ~
  mkdir training
  cd training
  ```
- In the menu top left, select File > Open Folder (<kbd>^</kbd> <kbd>⇧</kbd> <kbd>O</kbd>)
- Enter `/home/gitpod/training`
- GitPod will probably reload the browser tab
- The file explorer on the left should now have an expandable box with the title `TRAINING`

> 💡 ♻️ As you create new files you should see the explorer populate with them. If you don't see a file you expect, try hovering over the `TRAINING` section title and clicking the <kbd>⟳</kbd> icon that appears to refresh the file list.

### Local installation

If you prefer, you can install the nf-core/tools package locally.
It is a Python package, available via the [Python Package Index](https://pypi.org/project/nf-core/) (`pip`) and [bioconda](https://bioconda.github.io/recipes/nf-core/README.html) (`conda` / `mamba`).
There is a docker container available, however for the purposes of this tutorial it not recommended to use this.

To install from PyPI:

```bash
pip install nf-core
```

To get the very latest development version of the code:

```bash
pip install --upgrade --force-reinstall git+https://github.com/nf-core/tools.git@dev
```

If using conda, first set up Bioconda as described in the [bioconda docs](https://bioconda.github.io/#usage) (especially setting the channel order) and then install nf-core:

```bash
conda install nf-core
```

### First run of nf-core

Whether running using GitPod or locally, you can confirm that nf-core is correctly installed by viewing the command line help:

```bash
nf-core --help
```

You should get something that looks like the following output:

![nf-core --help](/assets/markdown_assets/developers/creating_with_nf_core/nfcore_help.svg)

The first set of subcommands are typically useful for running pipelines, the second are for developing pipelines.

You can try out some commands, for example listing available nf-core pipelines:

![nf-core list](/assets/markdown_assets/developers/creating_with_nf_core/nfcore_list.svg)

In this tutorial we will focus on creating a pipeline, but please do look at the functionality that nf-core/tools provides to you as a user - especially `nf-core launch` and `nf-core download`.

## Creating a pipeline

To get started with your new pipeline, run the create command:

```bash
nf-core create
```

Although you can provide options on the command line, it's easiest to use the interactive prompts.

![nf-core create](/assets/markdown_assets/developers/creating_with_nf_core/nfcore_create.svg)

Follow the instructions and you should see a new pipeline appear in your file explorer.

Let's move into the new pipeline directory in the terminal:

```bash
# NB: The path will vary according to what you called your pipeline!
cd nf-core-demo/
```

### Pipeline git repo

The `nf-core create` command has made a fully fledged pipeline for you.
Before getting too carried away looking at all of the files, note that it has also initiated a git repository:

![git status](/assets/markdown_assets/developers/creating_with_nf_core/git_status.svg)

```console
$ git status
On branch master
nothing to commit, working tree clean
```

It's actually created three branches for you:

![git branch](/assets/markdown_assets/developers/creating_with_nf_core/git_branch.svg)

```console
$ git branch
  TEMPLATE
  dev
* master
```

Each have the same initial commit, with the vanilla template:

![git log](/assets/markdown_assets/developers/creating_with_nf_core/git_log.svg)

```console
$ git log --oneline
bb0db93 (HEAD -> master, dev, TEMPLATE) initial template build from nf-core/tools, version 2.5.1
```

This is important, because this shared git history with unmodified nf-core template in the `TEMPLATE` branch is how the [nf-core automated template synchronisation](https://nf-co.re/docs/contributing/sync) works (see the docs for more details).

The main thing to remember with this is that:

- When creating a new repository on <https://github.com> or equivalent, don't initialise it - leave it bare and push everything from your local clone
- Develop your code on either the `master` or `dev` branches and leave `TEMPLATE` alone.

### Testing the new pipeline

The new pipeline should run with Nextflow, right out of the box.
Let's try:

```bash
cd ../
nextflow run nf-core-demo/ -profile test,docker --outdir test_results
```

![nextflow run](/assets/markdown_assets/developers/creating_with_nf_core/nextflow_run_1.svg)

```console
$ cd ../

$ nextflow run nf-core-demo/ -profile test,docker --outdir test_results
N E X T F L O W  ~  version 22.04.0
Launching `nf-core-demo/main.nf` [ridiculous_colden] DSL2 - revision: 83bcba00f4


------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/demo v1.0dev
------------------------------------------------------
Core Nextflow options
  runName                   : ridiculous_colden
  containerEngine           : docker
  launchDir                 : /home/gitpod/training
  workDir                   : /home/gitpod/training/work
  projectDir                : /home/gitpod/training/nf-core-demo
  userName                  : gitpod
  profile                   : test,docker
  configFiles               : /home/gitpod/training/nf-core-demo/nextflow.config

Input/output options
  input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
  outdir                    : test_results

Reference genome options
  genome                    : R64-1-1
  fasta                     : s3://ngi-igenomes/igenomes/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/WholeGenomeFasta/genome.fa

Institutional config options
  config_profile_name       : Test profile
  config_profile_description: Minimal test dataset to check pipeline function

Max job request options
  max_cpus                  : 2
  max_memory                : 6.GB
  max_time                  : 6.h

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------
If you use nf-core/demo for your analysis please cite:

* The nf-core framework
  https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
  https://github.com/nf-core/demo/blob/master/CITATIONS.md
------------------------------------------------------
Downloading plugin nf-amazon@1.7.1
executor >  local (7)
[2f/f9d04a] process > NFCORE_DEMO:DEMO:INPUT_CHECK:SAMPLESHEET_CHECK (samplesheet_test_illumina_amplicon.csv) [100%] 1 of 1 ✔
[b0/f47909] process > NFCORE_DEMO:DEMO:FASTQC (SAMPLE1_PE_T1)                                                 [100%] 4 of 4 ✔
[f7/9c185c] process > NFCORE_DEMO:DEMO:CUSTOM_DUMPSOFTWAREVERSIONS (1)                                        [100%] 1 of 1 ✔
[cc/f88b11] process > NFCORE_DEMO:DEMO:MULTIQC                                                                [100%] 1 of 1 ✔
-[nf-core/demo] Pipeline completed successfully-
Completed at: 03-Oct-2022 21:38:18
Duration    : 1m 23s
CPU hours   : (a few seconds)
Succeeded   : 7
```

Sure enough, our new pipeline has run FastQC and MultiQC on a selection of test data.

### Customising the template

In many of the files generated by the nf-core template, you'll find code comments that look like this:

```groovy
// TODO nf-core: Do something here
```

These are markers to help you get started with customising the template code as you write your pipeline.
Editor tools such as [Todo tree](https://marketplace.visualstudio.com/items?itemName=Gruntfuggly.todo-tree) help you easily navigate these and work your way through them.

### Linting your pipeline

Customising the template is part of writing your new pipeline.
However, not _all_ files should be edited - indeed, nf-core strives to promote standardisation amongst pipelines.

To try to keep pipelines up to date and using the same code where possible, we have an automated code linting tool for nf-core pipelines.
Running `nf-core lint` will run a comprehensive test suite against your pipeline.

Linting tests can have one of four statuses: pass, ignore, warn or fail.
For example, you may see a large number of _warnings_ about `TODO` comments, letting you know that you haven't finished setting up your new pipeline.
Warnings are ok at this stage, but should be cleared up before a pipeline release.

Failures are more serious however, and will typically prevent pull-requests from being merged.
For example, if you edit `CODE_OF_CONDUCT.md`, which should match the template, you'll get a pipeline lint test failure:

```console
$ nf-core lint

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.5.1 - https://nf-co.re


INFO     Testing pipeline: .
╭─ [✗] 1 Pipeline Test Failed ─────────────────────────────────────────╮
│                                                                      │
│ files_unchanged: CODE_OF_CONDUCT.md does not match the template      │
│                                                                      │
╰──────────────────────────────────────────────────────────────────────╯
```

### Continuous integration testing

Whilst it's helpful to be able to run the nf-core lint tests locally, their real strength is the combination with CI (continuous integration) testing.
By default, nf-core pipelines are configured to run CI tests on GitHub every time you push to the repo or open a pull request.
The same `nf-core lint` command runs on your code on the automated machines.
If there are any failures, they will be reported with a ❌ and you will not be able to merge the pull-request until you push more commits that fix those failures.

These automated tests allow us to maintain code quality in a scalable and standardised way across the community.

### Configuring linting

In some special cases (especially using the nf-core template outside of the nf-core community) you may find that you want to ignore certain linting failures.

To do so, edit `.nf-core.yml` in the root of the repositry.
For example, to ignore the tests triggering warnings in the example above, you can add:

```yml
lint:
  pipeline_todos: false
  files_unchanged:
    - CODE_OF_CONDUCT.md
```

Please see the linting documentation for specific details of how to configure different tests.

## Nextflow Schema

All nf-core pipelines can be run with `--help` to see usage instructions.
We can try this with the demo pipeline that we just created:

```console
$ nextflow run nf-core-demo/ --help
N E X T F L O W  ~  version 22.04.0
Launching `nf-core-demo/main.nf` [compassionate_waddington] DSL2 - revision: 83bcba00f4


------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/demo v1.0dev
------------------------------------------------------
Typical pipeline command:

  nextflow run nf-core/demo --input samplesheet.csv --genome GRCh37 -profile docker

Input/output options
  --input                      [string]  Path to comma-separated file containing information about the samples in the experiment.
  --outdir                     [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud
                                         infrastructure.
  --email                      [string]  Email address for completion summary.
  --multiqc_title              [string]  MultiQC report title. Printed as page header, used for filename if not otherwise specified.

Reference genome options
  --genome                     [string]  Name of iGenomes reference.
  --fasta                      [string]  Path to FASTA genome file.

!! Hiding 22 params, use --show_hidden_params to show them !!
------------------------------------------------------
If you use nf-core/demo for your analysis please cite:

* The nf-core framework
  https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
  https://github.com/nf-core/demo/blob/master/CITATIONS.md
------------------------------------------------------
```

Here we get a rich help output, with command line parameter variable types, and help text.
However, the Nextflow syntax in `nextflow.config` does not allow for this kind of rich metadata natively.

To provide this, we created a standard for describing pipeline parameters in a file in the pipeline root called `nextflow_schema.json`.
These files are written using [JSON Schema](https://json-schema.org/), making them compatible with many other tools and validation libraries.

By describing our workflow parameters in this file we can do a lot of new things:

- Generate command line help text and web based documentation pages
- Validate pipeline inputs
- Create rich pipeline launch interfaces

Indeed, if you try to run the new pipeline without the required `--outdir` parameter, you will quickly get an error:

```console
$ nextflow run nf-core-demo/ -profile test,docker
Picked up JAVA_TOOL_OPTIONS:  -Xmx3489m
N E X T F L O W  ~  version 22.04.0
Launching `nf-core-demo/main.nf` [determined_lamarck] DSL2 - revision: 83bcba00f4

ERROR: Validation of pipeline parameters failed!


* Missing required parameter: --outdir

```

### Working with Nextflow schema

If you peek inside the `nextflow_schema.json` file you will see that it is quite an intimidating thing.
The file is large and complex, and very easy to break if edited manually.

Thankfully, we provide a user-friendly tool for editing this file: `nf-core schema build`.

To see this in action, let's add some new parameters to `nextflow.config`:

```groovy
params {
    demo                       = 'param-value-default'
    foo                        = null
    bar                        = false
    baz                        = 12
    // rest of the config file..
```

Then run `nf-core schema build` - it should prompt you to add each new:

```console
$ nf-core schema build

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.5.1 - https://nf-co.re


INFO     [✓] Default parameters match schema validation
INFO     [✓] Pipeline schema looks valid (found 28 params)
✨ Found 'params.demo' in the pipeline config, but not in the schema. Add to pipeline schema? [y/n]: y
✨ Found 'params.foo' in the pipeline config, but not in the schema. Add to pipeline schema? [y/n]: y
✨ Found 'params.bar' in the pipeline config, but not in the schema. Add to pipeline schema? [y/n]: y
✨ Found 'params.baz' in the pipeline config, but not in the schema. Add to pipeline schema? [y/n]: y
INFO     Writing schema with 32 params: './nextflow_schema.json'
🚀  Launch web builder for customisation and editing? [y/n]: y
```

Select `y` on the final prompt to launch a web browser to edit your schema graphically.

> ⚠️ When using GitPod, a slightly odd thing can happen at this point.
> GitPod may launch the schema editor with a terminal-based web browser called `lynx`.
> You'll see the terminal output suddenly fill with a text-based representation of the nf-core website.
> Press <kbd>Q</kbd> followed by <kbd>Y</kbd> to confirm to exit `lyx`.
> You should still see nf-core/tools running on the command line.
> Copy the URL that was printed to the terminal and open this in a new browser tab (or <kbd>alt</kbd> + click it).

Here in the schema editor you can edit:

- Description and help text
- Type (string / boolean / integer etc)
- Grouping of parameters
- Whether a parameter is required, or hidden from help by default
- Enumerated values (choose from a list)
- Min / max values for numeric types
- Regular expressions for validation
- Special formats for strings, such as `file-path`
- Additional fields for files such as `mime-type`

## nf-core modules
