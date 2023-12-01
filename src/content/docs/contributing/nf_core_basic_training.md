---
title: Basic training to create an nf-core pipeline
subtitle: A guide to create Nextflow pipelines using nf-core tools
---

# Introduction

This training course aims to demonstrate how to build an nf-core pipeline using the nf-core pipeline template and nf-core modules and subworkflows as well as custom, local modules. Be aware that we are not going to explain any fundamental Nextflow concepts, as such we advise anyone taking this course to have completed the [Basic Nextflow Training Workshop](https://training.nextflow.io/).

```md
During this course we are going to build a Simple RNA-Seq workflow.
This workflow is by no means ment to be a useful bioinformatics workflow,
but should only teach the objectives of the course, so please,
**DO NOT use this workflow to analyse RNA sequencing data**!
```

## Overview

### Layout of the pipeline

The course is going to build an (totally unscientific and useless) RNA seq pipeline that does the following:

1. Indexing of a transcriptome file
2. Quality control
3. Quantification of transcripts
4. [whatever the custom script does]
5. Generation of a MultiQC report

### Outline of the Course

The following sections will be handled in the course:

**1. Setting up the gitpod environment for the course**

The course is using gitpod in order to avoid the time expense for downloading and installing tools and data.

**2. Exploring the nf-core tools command**

A very basic walk-through of what can be done with nf-core tools

**3. Creating a new nf-core pipeline from the nf-core template**

**4. Exploring the nf-core template**

a) The git repository

b) running the pipeline

c) linting the pipeline

d) walk-through of the template files

**5. Building a nf-core pipeline using the template**

a) Adding a nf-core module to your pipeline

b) Adding a local custom module to your pipeline

c) Working with Nextflow schema

d) Linting your modules

## Prerequisites

## Follow the training videos

This training can be followed either based on this documentation alone, or via a training video hosted on youtube. You can find the youtube video in the Youtube playlist below:

(no such video yet)

# Using gitpod

For this tutorial we are going to use Gitpod, which is best for first-timers as this platform contains all the programs and data required.
Gitpod will contain a preconfigured Nextflow development environment and has the following requirements:

- A GitHub account
- Web browser (Google Chrome, Firefox)
- Internet connection

Simply click the link and log in using your GitHub account to start the tutorial:

<p class="text-center">
  <a href="https://www.gitpod.io/#https://github.com/nf-core/basic_training" class="btn btn-lg btn-success" target="_blank">
    Launch GitPod
  </a>
</p>

For more information about gitpod, including how to make your own gitpod environement, see the gitpod bytesize talk on youtube (link to the bytesize talk)

## Explore your Gitpod interface

You should now see something similar to the following:

(insert Gitpod welcome image)

- **The sidebar** allows you to customize your Gitpod environment and perform basic tasks (copy, paste, open files, search, git, etc.). Click the Explorer button to see which files are in this repository.
- **The terminal** allows you to run all the programs in the repository. For example, both `nextflow` and `docker` are installed and can be executed.
- **The main window** allows you to view and edit files. Clicking on a file in the explorer will open it within the main window. You should also see the nf-training material browser (<https://training.nextflow.io/>).

To test that the environment is working correctly, type the following into the terminal:

```bash
nextflow info
```

This should come up with the Nextflow version and runtime information:

```
Version: 23.10.0 build 5889
Created: 15-10-2023 15:07 UTC (15:07 GMT)
System: Linux 6.1.54-060154-generic
Runtime: Groovy 3.0.19 on OpenJDK 64-Bit Server VM 17.0.8-internal+0-adhoc..src
Encoding: UTF-8 (UTF-8)
```

## Reopening a Gitpod session

You can reopen an environment from <https://gitpod.io/workspaces>. Find your previous environment in the list, then select the ellipsis (three dots icon) and select Open.

If you have saved the URL for your previous Gitpod environment, you can simply open it in your browser.

Alternatively, you can start a new workspace by following the Gitpod URL: <https://gitpod.io/#https://github.com/nextflow-io/training>

If you have lost your environment, you can find the main scripts used in this tutorial in the `nf-training` directory.

## Saving files from Gitpod to your local machine

To save any file from the explorer panel, right-click the file and select Download.

# Explore nf-core/tools

The nf-core/tools package is already installed in the gitpod environment. Now you can check out which pipelines, subworkflows and modules are available via tools. To see all available commands of nf-core tools, run the following:

```bash
nf-core --help
```

We will touch on most of the commands for developers later throughout this tutorial.

# Create a pipeline from template

To get started with your new pipeline, run the create command:

```bash
nf-core create
```

This should open a command prompt similar to this:

```

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.10 - https://nf-co.re


? Workflow name demotest
? Description test pipeline for demo
? Author FBo
? Do you want to customize which parts of the template are used? No
INFO     Creating new nf-core pipeline: 'nf-core/demotest'
INFO     Initialising pipeline git repository
INFO     Done. Remember to add a remote and push to GitHub:
          cd /workspace/basic_training/nf-core-demotest
          git remote add origin git@github.com:USERNAME/REPO_NAME.git
          git push --all origin
INFO     This will also push your newly created dev branch and the TEMPLATE branch for syncing.
INFO     !!!!!! IMPORTANT !!!!!!
         If you are interested in adding your pipeline to the nf-core community,
         PLEASE COME AND TALK TO US IN THE NF-CORE SLACK BEFORE WRITING ANY CODE!

         Please read: https://nf-co.re/developers/adding_pipelines#join-the-community
```

Although you can provide options on the command line, it’s easiest to use the interactive prompts. For now we are assuming that we want to create a new nf-core pipeline, so we chose not to customize the template.
It is possible to use nf-core tools for non-nf-core pipelines, but the setup of such pipelines will be handled in a later chapter # ARE WE GOING TO DO THIS?

### Pipeline git repo

The nf-core create command has made a fully fledged pipeline for you. Before getting too carried away looking at all of the files, note that it has also initiated a git repository:

```bash
cd nf-core-demotest
git status
```

```
On branch master
nothing to commit, working tree clean
```

It’s actually created three branches for you:

```bash
git branch
```

```
  TEMPLATE
  dev
* master
```

Each have the same initial commit, with the vanilla template:

```bash
git log
```

```
commit 77e77783aab19e47ce6b2d736d766fbef2483de8 (HEAD -> master, dev, TEMPLATE)
Author: Franziska Bonath <franziska.bonath@scilifelab.se>
Date:   Tue Oct 31 15:21:33 2023 +0000

    initial template build from nf-core/tools, version 2.10
```

This is important, because this shared git history with unmodified nf-core template in the TEMPLATE branch is how the nf-core automated template synchronisation works (see the docs for more details).

The main thing to remember with this is that:

Never make changes to the TEMPLATE branch, otherwise it will interfere with the synchronisation of nf-core updates.

Ideally code should be developed on feature branches (i.e. a new branch made with `git checkout -b my_new_feature`), and when ready merged into the `dev` branch upon a successful code review. The `dev` branch is then merged to the `master` branch when a stable release of the workflow is ready to be made.

When creating a new repository on GitHub, create it as an empty repository without a README or any other file. Then push the repo with the template of your new pipeline from your local clone.

( OLD: When creating a new repository on https://github.com or equivalent, don’t initialise it - leave it bare and push everything from your local clone
Develop your code on either the master or dev branches and leave TEMPLATE alone.)

## Run the new pipeline

The new pipeline should run with Nextflow, right out of the box. Let’s try:

```bash
cd ../
nextflow run nf-core-demotest/ -profile test,docker --outdir test_results
```

This basic template pipeline contains already the FastQC and MultiQC modules, which do run on a selection of test data.

## Template code walk through

Now let us have a look at the files that were generated within the `nf-core-demotest` directory when we created this pipeline. You can see all files and directories either on the left hand side in the Explorer, or by running the command:

```bash
cd nf-core-demotest
tree
```

```
.
├── assets
│   ├── adaptivecard.json
│   ├── email_template.html
│   ├── email_template.txt
│   ├── methods_description_template.yml
│   ├── multiqc_config.yml
│   ├── nf-core-demotest_logo_light.png
│   ├── samplesheet.csv
│   ├── schema_input.json
│   ├── sendmail_template.txt
│   └── slackreport.json
├── bin
│   └── check_samplesheet.py
├── CHANGELOG.md
├── CITATIONS.md
├── CODE_OF_CONDUCT.md
├── conf
│   ├── base.config
│   ├── igenomes.config
│   ├── modules.config
│   ├── test.config
│   └── test_full.config
├── docs
│   ├── images
│   │   ├── mqc_fastqc_adapter.png
│   │   ├── mqc_fastqc_counts.png
│   │   ├── mqc_fastqc_quality.png
│   │   ├── nf-core-demotest_logo_dark.png
│   │   └── nf-core-demotest_logo_light.png
│   ├── output.md
│   ├── README.md
│   └── usage.md
├── lib
│   ├── nfcore_external_java_deps.jar
│   ├── NfcoreTemplate.groovy
│   ├── Utils.groovy
│   ├── WorkflowDemotest.groovy
│   └── WorkflowMain.groovy
├── LICENSE
├── main.nf
├── modules
│   ├── local
│   │   └── samplesheet_check.nf
│   └── nf-core
│       ├── custom
│       │   └── dumpsoftwareversions
│       │       ├── main.nf
│       │       ├── meta.yml
│       │       └── templates
│       │           └── dumpsoftwareversions.py
│       ├── fastqc
│       │   ├── main.nf
│       │   └── meta.yml
│       └── multiqc
│           ├── main.nf
│           └── meta.yml
├── modules.json
├── nextflow.config
├── nextflow_schema.json
├── pyproject.toml
├── README.md
├── subworkflows
│   └── local
│       └── input_check.nf
├── tower.yml
└── workflows
    └── demotest.nf
```

These are the files in detail:

1. **main.nf**

   This file contains the main nextflow pipeline code. Mostly this file is not touched.

2. **workflows/demotest.nf**

   This file is where the pipeline is going to be assembled. It connects the different modules and subworkflows.

3. **CHANGELOG.md, CODE_OF_CONDUCT.md, LICENSE, README.md, CITATIONS.md**

   These are standard files created for github repositories. As a default, this pipeline will be under an MIT licence. The CODE_OF_CONDUCT is specific for the nf-core community.

4. **assets/**

   This directory contains different templates such as email templates or the MultiQC config. In this course contents of this directory can be largely ignored.

5. **bin/**

   The `bin` directory contains custom executable scripts, and is automatically added to the `PATH` by Nextflow allowing these scripts to become findable by Nextflow. As such, they can be called by name without using their absolute or relative path. The python script `check_samplesheet.py` is part of the nf-core template, since typically, nf-core pipelines require a samplesheet as one of their inputs.

6. **conf/** and **nextflow.config**

   The `nextflow.config` file is the main config file. In addition to supplying default parameters, it imports all the configurations in `conf/`. Importantly, `conf/` contains the `test.config` file, which is used for pipeline testing. In this course we are not going to touch config files, but they have been extensively covered in the following bytesize talks: [How nf-core configs work (nf-core/bytesize #2)](https://www.youtube.com/watch?v=cXBYusdjrc0&list=PL3xpfTVZLcNiF4hgkW0yXeNzr0d35qlIB&index=8&pp=gAQBiAQB), [Making a new institutional config profile (nf-core/bytesize #10)](https://www.youtube.com/watch?v=Ym1s6sKGzkw&list=PL3xpfTVZLcNiF4hgkW0yXeNzr0d35qlIB&index=9&pp=gAQBiAQB), [nf-core/bytesize: Using nf-core configs in custom pipelines](https://www.youtube.com/watch?v=zgcrI_0SUgg&list=PL3xpfTVZLcNiF4hgkW0yXeNzr0d35qlIB&index=40&pp=gAQBiAQB)

7. **docs/**

   This directory contains additional information to the README file. The most important files are the `output.md` and the `usage.md` files. `usage.md` should describe what exactly is needed to run the pipeline and `output.md` should describe all outputs that can be expected. Importantly, for nf-core pipelines, the information from these two files will automatically be displayed on the nf-core website page of the pipeline.

8. **lib/**

   This directory contains groovy functions and classes that are imported into the `main.nf` file to provide additional functionality not native to Nextflow.

9. **modules/local**

   This is where all your custom non-nf-core modules go. We will cover when and how to make local modules later in the course.

10. **modules/nf-core**

    All nf-core modules that are installed using the nf-core tooling will automatically show up in this directory. Keep them here, it is important for automatic updates.

11. **modules.json**

This file keeps track of modules installed using nf-core tools from the nf-core/modules repository. This file should only be updated using nf-core tools, and never manually.

12. **nextflow_schema.json**

    This file hosts all the parameters for the pipeline. Any new parameter should be added to this file using the `nf-core schema build` command. Similar to `output.md` and `usage.md`, the contents of `nextflow_schema.json` will get displayed on the pipeline page of nf-core pipelines.

13. **pyproject.toml**
14. **subworkflows/local**
15. **tower.yml**
16. **hidden directories and files**

    a) _.devcontainer/devcontainer.json_

    b) _.github/_

    Files in here are used for Continuous integration tests (CI) with github actions as well as other github related defaults, such as a template for issues. We will not touch on these in the course.

    c) _.gitignore_

    d) _.editorconfig_

    e) _.gitpod.yml_

This file provides settings to create a Cloud development environment in your browser using Gitpod. It comes installed with the tools necessary to develop and test nf-core pipelines, modules, and subworkflows, allowing you to develop from anywhere without installing anything locally.

    f) _.nf-core.yml_

    g) _.pre-commit-config.yaml_

    h) _.prettierignore_

    i) _.prettierrc.yml_

## Customising the template

In many of the files generated by the nf-core template, you’ll find code comments that look like this:

```
// TODO nf-core: Do something here
```

These are markers to help you get started with customising the template code as you write your pipeline. Editor tools such as Todo tree help you easily navigate these and work your way through them.

## Linting your pipeline

Customising the template is part of writing your new pipeline. However, not all files should be edited - indeed, nf-core strives to promote standardisation amongst pipelines.

To try to keep pipelines up to date and using the same code where possible, we have an automated code linting tool for nf-core pipelines. Running nf-core lint will run a comprehensive test suite against your pipeline:

```bash
cd nf-core-demotest/
nf-core lint
```

The first time you run this command it will download some modules and then perform the linting tests. Linting tests can have one of four statuses: pass, ignore, warn or fail. For example, at first you will see a large number of warnings about TODO comments, letting you know that you haven’t finished setting up your new pipeline. Warnings are ok at this stage, but should be cleared up before a pipeline release.

```
                                         ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.10 - https://nf-co.re


INFO     Testing pipeline: .

╭─ [!] 24 Pipeline Test Warnings ──────────────────────────────────────────────────────────────────────────────────────────────────╮
│                                                                                                                                  │
│ readme: README contains the placeholder zenodo.XXXXXXX. This should be replaced with the zenodo doi (after the first release).   │
│ pipeline_todos: TODO string in README.md: TODO nf-core:                                                                          │
│ pipeline_todos: TODO string in README.md: Include a figure that guides the user through the major workflow steps. Many nf-core   │
│ pipeline_todos: TODO string in README.md: Fill in short bullet-pointed list of the default steps in the pipeline                 │
│ pipeline_todos: TODO string in README.md: Describe the minimum required steps to execute the pipeline, e.g. how to prepare       │
│ samplesheets.                                                                                                                    │
│ pipeline_todos: TODO string in README.md: update the following command to include all required parameters for a minimal example  │
│ pipeline_todos: TODO string in README.md: If applicable, make list of people who have also contributed                           │
│ pipeline_todos: TODO string in README.md: Add citation for pipeline after first release. Uncomment lines below and update Zenodo │
│ doi and badge at the top of this file.                                                                                           │
│ pipeline_todos: TODO string in README.md: Add bibliography of tools and data used in your pipeline                               │
│ pipeline_todos: TODO string in main.nf: Remove this line if you don't need a FASTA file                                          │
│ pipeline_todos: TODO string in nextflow.config: Specify your pipeline's command line flags                                       │
│ pipeline_todos: TODO string in awsfulltest.yml: You can customise AWS full pipeline tests as required                            │
│ pipeline_todos: TODO string in ci.yml: You can customise CI pipeline run tests as required                                       │
│ pipeline_todos: TODO string in methods_description_template.yml: #Update the HTML below to your preferred methods description,   │
│ e.g. add publication citation for this pipeline                                                                                  │
│ pipeline_todos: TODO string in base.config: Check the defaults for all processes                                                 │
│ pipeline_todos: TODO string in base.config: Customise requirements for specific processes.                                       │
│ pipeline_todos: TODO string in test.config: Specify the paths to your test data on nf-core/test-datasets                         │
│ pipeline_todos: TODO string in test.config: Give any required params for the test so that command line flags are not needed      │
│ pipeline_todos: TODO string in test_full.config: Specify the paths to your full test data ( on nf-core/test-datasets or directly │
│ in repositories, e.g. SRA)                                                                                                       │
│ pipeline_todos: TODO string in test_full.config: Give any required params for the test so that command line flags are not needed │
│ pipeline_todos: TODO string in output.md: Write this documentation describing your workflow's output                             │
│ pipeline_todos: TODO string in usage.md: Add documentation about anything specific to running your pipeline. For general topics, │
│ please point to (and add to) the main nf-core website.                                                                           │
│ pipeline_todos: TODO string in WorkflowDemotest.groovy: Optionally add in-text citation tools to this list.                      │
│ pipeline_todos: TODO string in WorkflowMain.groovy: Add Zenodo DOI for pipeline after first release                              │
│                                                                                                                                  │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

╭─ [!] 3 Module Test Warnings ─────────────────────────────────────────────────────────────────────────────────────────────────────╮
│                                           ╷                                          ╷                                           │
│ Module name                               │ File path                                │ Test message                              │
│╶──────────────────────────────────────────┼──────────────────────────────────────────┼──────────────────────────────────────────╴│
│ custom/dumpsoftwareversions               │ modules/nf-core/custom/dumpsoftwarevers… │ New version available                     │
│ fastqc                                    │ modules/nf-core/fastqc                   │ New version available                     │
│ multiqc                                   │ modules/nf-core/multiqc                  │ New version available                     │
│                                           ╵                                          ╵                                           │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭───────────────────────╮
│ LINT RESULTS SUMMARY  │
├───────────────────────┤
│ [✔] 183 Tests Passed  │
│ [?]   0 Tests Ignored │
│ [!]  27 Test Warnings │
│ [✗]   0 Tests Failed  │
╰───────────────────────╯
```

Failures are more serious however, and will typically prevent pull-requests from being merged. For example, if you edit CODE_OF_CONDUCT.md, which should match the template, you’ll get a pipeline lint test failure:

```bash
echo "Edited" >> CODE_OF_CONDUCT.md
nf-core lint
```

```
[...]
╭─ [✗] 1 Pipeline Test Failed ─────────────────────────────────────────────────────╮
│                                                                                  │
│ files_unchanged: CODE_OF_CONDUCT.md does not match the template                  │
│                                                                                  │
╰──────────────────────────────────────────────────────────────────────────────────╯
[...]
╭───────────────────────╮
│ LINT RESULTS SUMMARY  │
├───────────────────────┤
│ [✔] 182 Tests Passed  │
│ [?]   0 Tests Ignored │
│ [!]  27 Test Warnings │
│ [✗]   1 Test Failed   │
╰───────────────────────╯
[...]
```

# Building a pipeline from (existing) components

Nextflow pipelines can be build in a very modular fashion. In nf-core, we have simple building blocks available: nf-core/modules. They are wrappers around usually individual tools. In addition, we have subworkflows: smaller pre-build pipeline chunks. You can think about the modules as Lego bricks and subworkflows as pre-build chunks that can be added to various sets. These components are centrally available for all Nextflow pipelines. To make working with them easy, we have can use `nf-core/tools`

## Adding an existing nf-core module

### Identify available nf-core modules

The nf-core pipeline template comes with a few nf-core/modules pre-installed. You can list these with the command below:

```bash
nf-core modules list local
```

```
                                         ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.10 - https://nf-co.re


INFO     Modules installed in '.':

┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
┃ Module Name                 ┃ Repository                 ┃ Version SHA                 ┃ Message                    ┃ Date       ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
│ custom/dumpsoftwareversions │ https://github.com/nf-cor… │ 911696ea0b62df80e900ef244d… │ Remove quay from           │ 2023-05-04 │
│                             │                            │                             │ biocontainers (#3380)      │            │
│ fastqc                      │ https://github.com/nf-cor… │ bd8092b67b5103bdd52e300f75… │ Add singularity.registry = │ 2023-07-01 │
│                             │                            │                             │ 'quay.io' for tests        │            │
│                             │                            │                             │ (#3499)                    │            │
│ multiqc                     │ https://github.com/nf-cor… │ 911696ea0b62df80e900ef244d… │ Remove quay from           │ 2023-05-04 │
│                             │                            │                             │ biocontainers (#3380)      │            │
└─────────────────────────────┴────────────────────────────┴─────────────────────────────┴────────────────────────────┴────────────┘

```

These version hashes and repository information for the source of the modules are tracked in the modules.json file in the root of the repo. This file will automatically be updated by nf-core/tools when you create, remove or update modules.

Let’s see if all of our modules are up-to-date:

```bash
nf-core modules update
```

```
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.10 - https://nf-co.re


? Update all modules or a single named module? All modules
? Do you want to view diffs of the proposed changes? No previews, just update everything
INFO     Updating 'nf-core/custom/dumpsoftwareversions'
INFO     Updating 'nf-core/fastqc'
INFO     Updating 'nf-core/multiqc'
INFO     Updates complete ✨
```

You can list all of the modules available on nf-core/modules via the command below but we have added search functionality to the nf-core website to do this too!

```bash
nf-core modules list remote
```

In addition, all modules are listed on the website: [https://nf-co.re/modules](https://nf-co.re/modules)

### Install a remote nf-core module

To install a remote nf-core module, you can first get information about a tool, including the installation command by executing:

```bash
nf-core modules info salmon/index
```

```
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.10 - https://nf-co.re


╭─ Module: salmon/index  ───────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ 🌐 Repository: https://github.com/nf-core/modules.git                                                                             │
│ 🔧 Tools: salmon                                                                                                                  │
│ 📖 Description: Create index for salmon                                                                                           │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
                             ╷                                                                                              ╷
 📥 Inputs                   │Description                                                                                   │Pattern
╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━╸
  genome_fasta  (file)       │Fasta file of the reference genome                                                            │
╶────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────────────────┼───────╴
  transcriptome_fasta  (file)│Fasta file of the reference transcriptome                                                     │
                             ╵                                                                                              ╵
                    ╷                                                                                                  ╷
 📤 Outputs         │Description                                                                                       │     Pattern
╺━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
  index  (directory)│Folder containing the star index files                                                            │      salmon
╶───────────────────┼───────────────────────────────────────────────────────────────────────���──────────────────────────┼────────────╴
  versions  (file)  │File containing software versions                                                                 │versions.yml
                    ╵                                                                                                  ╵

 💻  Installation command: nf-core modules install salmon/index
```

The output from the info command will among other things give you the nf-core/tools installation command, lets see what it is doing:

```bash
nf-core modules install salmon/index
```

```

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.10 - https://nf-co.re


╭─ Module: salmon/index  ───────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ 🌐 Repository: https://github.com/nf-core/modules.git                                                                             │
│ 🔧 Tools: salmon                                                                                                                  │
│ 📖 Description: Create index for salmon                                                                                           │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
                             ╷                                                                                              ╷
 📥 Inputs                   │Description                                                                                   │Pattern
╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━╸
  genome_fasta  (file)       │Fasta file of the reference genome                                                            │
╶────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────────────────┼───────╴
  transcriptome_fasta  (file)│Fasta file of the reference transcriptome                                                     │
                             ╵                                                                                              ╵
                    ╷                                                                                                  ╷
 📤 Outputs         │Description                                                                                       │     Pattern
╺━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
  index  (directory)│Folder containing the star index files                                                            │      salmon
╶───────────────────┼───────────────────────────────────────────────────────────────────────���──────────────────────────┼────────────╴
  versions  (file)  │File containing software versions                                                                 │versions.yml
                    ╵                                                                                                  ╵

 💻  Installation command: nf-core modules install salmon/index

gitpod /workspace/basic_training/nf-core-demotest (master) $ nf-core modules install salmon/index

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.10 - https://nf-co.re


INFO     Installing 'salmon/index'
INFO     Use the following statement to include this module:

 include { SALMON_INDEX } from '../modules/nf-core/salmon/index/main'
```

The module is now installed into the folder `modules/nf-core`. Now open the file `workflow/demotest.nf`. You will find already several `include` statements there from the installed modules (`MultiQC` and `FastQC`):

```bash

include { FASTQC  } from '../modules/nf-core/fastqc/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
```

Now add the above line underneath it:

```bash

include { FASTQC  } from '../modules/nf-core/fastqc/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { SALMON_INDEX } from '../modules/nf-core/salmon/index/main'

```

This makes the module now available in the workflow script and it can be called with the right input data.

<!-- TODO/TODISCUSS here the user now needs to know about how to get their fasta. We could do this here or add a new point for this above -->



(lots of steps missing here)
exercise to add a different module would be nice! => salmon/quant!
comparison to simple nextflow pipeline from the basic Nextflow training would be nice!)

## Adding a remote module

If there is no nf-core module available for the software you want to include, you can add the module to the nf-core/modules repository. It will then become available to the wider Nextflow Community. See how to [here](https://nf-co.re/docs/contributing/tutorials/dsl2_modules_tutorial). If the module is very pipeline specific, you can also add a local module. The nf-core tools package can aid in the generation of a module template. To add a bare-bone local module run the following:

```
nf-core modules create
```

Open ./modules/local/demo/module.nf and start customising this to your needs whilst working your way through the extensive TODO comments! For further help and guidelines for the modules code, check out the [modules specific documentation](https://nf-co.re/docs/contributing/tutorials/dsl2_modules_tutorial).

### Making a remote module for a custom script

## Lint all modules

As well as the pipeline template you can lint individual or all modules with a single command:

```
nf-core modules lint --all
```

# Nextflow Schema

All nf-core pipelines can be run with --help to see usage instructions. We can try this with the demo pipeline that we just created:

```
cd ../
nextflow run nf-core-demo/ --help
```

## Working with Nextflow schema

If you peek inside the nextflow_schema.json file you will see that it is quite an intimidating thing. The file is large and complex, and very easy to break if edited manually.

Thankfully, we provide a user-friendly tool for editing this file: nf-core schema build.

To see this in action, let’s add some new parameters to nextflow.config:

```
params {
demo = 'param-value-default'
foo = null
bar = false
baz = 12
// rest of the config file..
```

Then run nf-core schema build:

```
cd nf-core-demo/
nf-core schema build
```

The CLI tool should then prompt you to add each new parameter.
Here in the schema editor you can edit:

- Description and help text
- Type (string / boolean / integer etc)
- Grouping of parameters
- Whether a parameter is required, or hidden from help by default
- Enumerated values (choose from a list)
- Min / max values for numeric types
- Regular expressions for validation
- Special formats for strings, such as file-path
- Additional fields for files such as mime-type

```

```
