---
title: Basic training to create an nf-core pipeline
subtitle: A guide to create Nextflow pipelines using nf-core tools
---

## Scope

- How do I create a pipeline using nf-core tools?
- How do I incorporate modules from nf-core modules?
- How can I use custom code in my pipeline?

:::note

### Learning objectives

- The learner will create a simple pipeline using the nf-core template.
- The learner will identify key files in the pipeline.
- The learner will lint their pipeline code to identify work to be done.
- The learner will incorporate modules from nf-core/modules into their pipeline.
- The learner will add custom code as a local module into their pipeline.
- The learner will build an nf-core schema to describe and validate pipeline parameters.

:::

This training course aims to demonstrate how to build an nf-core pipeline using the nf-core pipeline template and nf-core modules as well as custom, local modules. Be aware that we are not going to explain any fundamental Nextflow concepts, as such we advise anyone taking this course to have completed the [Basic Nextflow Training Workshop](https://training.nextflow.io/).

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

1. **Setting up the gitpod environment for the course**

The course is using gitpod in order to avoid the time expense for downloading and installing tools and data.

2. **Exploring the nf-core tools command**

A very basic walk-through of what can be done with nf-core tools

3. **Creating a new nf-core pipeline from the nf-core template**

4. **Exploring the nf-core template**

   a) The git repository

   b) running the pipeline

   c) linting the pipeline

   d) walk-through of the template files

5. **Building a nf-core pipeline using the template**

   a) Adding a nf-core module to your pipeline

   b) Adding a local custom module to your pipeline

   c) Working with Nextflow schema

   d) Linting your modules

## Preparation

### Prerequisites

- Familiarity with Nextflow syntax and configuration.

### Follow the training videos

This training can be followed either based on this documentation alone, or via a training video hosted on youtube. You can find the youtube video in the Youtube playlist below:

(no such video yet)

### Gitpod

For this tutorial we will use Gitpod, which runs in the learners web browser. The Gitpod environment contains a preconfigured Nextflow development environment
which includes a terminal, file editor, file browser, Nextflow, and nf-core tools. To use Gitpod, you will need:

- A GitHub account
- Web browser (Google Chrome, Firefox)
- Internet connection

Click the link and log in using your GitHub account to start the tutorial:

<p class="text-center">
  <a href="https://www.gitpod.io/#https://github.com/nf-core/basic_training" class="btn btn-lg btn-success" target="_blank">
    Launch GitPod
  </a>
</p>

For more information about Gitpod, including how to make your own Gitpod environement, see the Gitpod bytesize talk on youtube (link to the bytesize talk),
check the [nf-core Gitpod documentation](gitpod/index) or [Gitpod's own documentation](https://www.gitpod.io/docs).

<details>
<summary> Expand this section for instructions to explore your Gitpod environment</summary>

#### Explore your Gitpod interface

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

#### Reopening a Gitpod session

When a Gitpod session is not used for a while, i.e., goes idle, it will timeout and close the interface.
You can reopen the environment from <https://gitpod.io/workspaces>. Find your previous environment in the list, then select the ellipsis (three dots icon) and select Open.

If you have saved the URL for your previous Gitpod environment, you can simply open it in your browser.

Alternatively, you can start a new workspace by following the Gitpod URL: <https://gitpod.io/#https://github.com/nextflow-io/training>

If you have lost your environment, you can find the main scripts used in this tutorial in the `nf-training` directory.

#### Saving files from Gitpod to your local machine

To save any file locally from the explorer panel, right-click the file and select Download.

</details>

## Explore nf-core/tools

The nf-core/tools package is already installed in the gitpod environment. Now you can check out which pipelines, subworkflows and modules are available via tools. To see all available commands of nf-core tools, run the following:

```bash
nf-core --help
```

We will touch on most of the commands for developers later throughout this tutorial.

## Create a pipeline from template

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

:::tip{title="Exercise 1 - Getting around the git environment"}

1. Create and switch to a new git branch called `demo`.
   <details>
   <summary>solution 1</summary>

   ```bash
   git checkout -b demo
   ```

   </details>

2. Display all available git branches.
   <details>
   <summary>solution 2</summary>

   ```bash
   git branch
   ```

   </details>

3. Create a directory within the new pipeline directory called `results` and add it to the `.gitignore` file.
   <details>
   <summary>solution 3</summary>

   ```bash
   mkdir results
   ```

   ```groovy title=".gitignore"
   .nextflow*
    work/
    data/
    results/
    .DS_Store
    testing/
    testing*
    *.pyc
    results/
   ```

   </details>

4. Commit the changes you have made.
   <details>
   <summary>solution 4</summary>

   ```bash
   git add .
   git commit -m "creating results dir and adding it to gitignore"
   ```

   </details>
   :::

### Run the new pipeline

The new pipeline should run with Nextflow, right out of the box. Let’s try:

```bash
cd ../
nextflow run nf-core-demotest/ -profile test,docker --outdir test_results
```

This basic template pipeline contains already the FastQC and MultiQC modules, which do run on a selection of test data.

### Template code walk through

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

:::tip{title="Exercise 2 - Test your knowledge of the nf-core pipeline structure"}

1. In which directory can you find the main script of the nf-core module `fastqc`
   <details>
      <summary>solution 1</summary>

   ```
   modules/nf-core/fastqc/
   ```

      </details>

2. Which file contains the main workflow of your new pipeline?
   <details>
      <summary>solution 2</summary>

   ```
   workflows/demotest.nf
   ```

      </details>

3. `check_samplesheet.py` is a script that can be called by any module of your pipeline, where is it located?
   <details>
      <summary>solution 3</summary>

   ```
   bin/
   ```

   This directory can also contain a custom scripts that you may wish to call from within a custom module.

      </details>

[MORE QUESTIONS CAN BE ADDED HERE]
:::

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

:::tip{title="Exercise 3 - ToDos and linting"}

1. Add the following bullet point list to the README file, where the ToDo indicates to describe the default steps to execute the pipeline

   ```groovy title="pipeline overview"
   - Indexing of a transcriptome file
   - Quality control
   - Quantification of transcripts
   - [whatever the custom script does]
   - Generation of a MultiQC report
   ```

   <details>
      <summary>solution 1</summary>

   ```bash title="README.md"
   [...]

   ## Introduction

   **nf-core/a** is a bioinformatics pipeline that ...

   <!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You\'re giving an overview to someone new to nf-core here,
   in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction -->

   <!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.
     -->

   Default steps:
      - Indexing of a transcriptome file
    - Quality control
    - Quantification of transcripts
    - [whatever the custom script does]
    - Generation of a MultiQC report

   1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
   2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

   [...]

   ```

      </details>

2. Lint the changes you have made
   <details>
      <summary>solution 2</summary>

   ```bash
   nf-core lint
   ```

   You should see that we now get one less `warning` in our lint overview, since we removed one of the TODO items.

      </details>

3. Commit your changes
   <details>
      <summary>solution 3</summary>

   ```bash
   git add .
   git commit -m "adding pipeline overview to pipeline README"
   ```

      </details>

   :::

## Adding Modules to a pipeline

### Adding an existing nf-core module

#### Identify available nf-core modules

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

#### Install a remote nf-core module

To install a remote nf-core module from the website, you can first get information about a tool, including the installation command by executing:

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

╭─ Module: salmon/index ───────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ 🌐 Repository: https://github.com/nf-core/modules.git │
│ 🔧 Tools: salmon │
│ 📖 Description: Create index for salmon │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╷ ╷
📥 Inputs │Description │Pattern
╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━╸
genome_fasta (file) │Fasta file of the reference genome │
╶────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────────────────┼───────╴
transcriptome_fasta (file)│Fasta file of the reference transcriptome │
╵ ╵
╷ ╷
📤 Outputs │Description │ Pattern
╺━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
index (directory)│Folder containing the star index files │ salmon
╶───────────────────┼───────────────────────────────────────────────────────────────────────���──────────────────────────┼────────────╴
versions (file) │File containing software versions │versions.yml
╵ ╵

💻 Installation command: nf-core modules install salmon/index

```

:::tip{title="Exercise 4 - Identification of available nf-core modules"}

1. Get information abou the nf-core module `salmon/quant`.
   <details>
      <summary>solution 1</summary>

   ```
   nf-core modules info salmon/quant
   ```

    </details>

2. Is there any version of `salmon/quant` already installed locally?
   <details>
      <summary>solution 2</summary>

   ```
   nf-core modules list local
   ```

   If `salmon/quant` is not listed, there is no local version installed.

      </details>
   :::

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

(lots of steps missing here)
exercise to add a different module would be nice! => salmon/quant!
comparison to simple nextflow pipeline from the basic Nextflow training would be nice!)

:::tip{title="Exercise 5 - Installing a remote module from nf-core"}

1.  Install the nf-core module `adapterremoval`
    <details>
       <summary>solution 1</summary>

    ```bash
    nf-core modules install adapterremoval
    ```

       </details>

2.  Which file(s) were/are added and what does it / do they do?
    <details>
       <summary>solution 2</summary>

    ```
    Installation added the module directory `/workspace/basic_training/nf-core-demotest/modules/nf-core/adapterremoval`:
    .
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
       ├── main.nf.test
       ├── main.nf.test.snap
       └── tags.yml

    The `test` directory contains all information required to perform basic tests for the module, it rarely needs to be changed. `main.nf` is the main workflow file that contains the module code. All input and output variables of the module are described in the `meta.yml` file, whereas the `environment.yml` file contains the dependancies of the module.
    ```

       </details>

3.  Import the installed `adapterremoval` pipeline into your main workflow.
    <details>
       <summary>solution 3</summary>

    ```bash title="workflows/demotest.nf"
      [...]
       /*
       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       */

       include { FASTQC                 } from '../modules/nf-core/fastqc/main'
       include { MULTIQC                } from '../modules/nf-core/multiqc/main'
       include { paramsSummaryMap       } from 'plugin/nf-validation'
       include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
       include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
       include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_demotest_pipeline'
       include { ADAPTERREMOVAL         } from '../modules/nf-core/adapterremoval/main'

      [...]

    ```

       </details>

4.  Call the `ADAPTERREMOVAL` process in your workflow
    <details>
       <summary>solution 4</summary>

    ```bash title="workflows/demotest.nf"
    [...]
    FASTQC (
      ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: ADAPTERREMOVAL
    //
    ADAPTERREMOVAL(

    )
    [...]
    ```

       </details>

5.  Add required parameters for `adapterremoval`to the `ADAPTERREMOVAL` process
    <details>
       <summary>solution 5</summary>

    `adapterremoval` requires three input channels: `meta`, `reads` and `adapterlist`, as outlined in the `meta.yml` of the module. `meta` and `reads` are typically given in one channel as a metamap, whereas the `adapterlist` will be it's own channel for which we should give a path. See here:

    ```bash title="adapterremoval/main.nf"
       [...]
       input:
       tuple val(meta), path(reads)
       path(adapterlist)
       [...]
    ```

    The meta map containing the metadata and the reads can be taken directly from the samplesheet as is the case for FastQC, therefore we can give it the input channel `ch_samplesheet`. The `adapterlist` could either be a fixed path, or a parameter that is given on the command line. For now, we will just add a dummy channel called `adapterlist` assuming that it will be a parameter given in the command line. With this, the new module call for adapterremoval looks as follows:

    ```bash title="workflows/demotest.nf"
    [...]
    //
    // MODULE: ADAPTERREMOVAL
    //
    ADAPTERREMOVAL(
       ch_samplesheet
       params.adapterlist
    )
    [...]
    ```

       </details>

6.  Add the input parameter `adapterlist`
    <details>
       <summary>solution 7</summary>
      In order to use `params.adapterlist` we need to add the parameter to the `nextflow.config`.

    ```bash title="nextflow.config"
       /// Global default params, used in configs
       params {

       /// TODO nf-core: Specify your pipeline's command line flags
       /// Input options
       input                      = null
       adapterlist                = null

      [...]
    ```

    Then use the `nf-core schema build` tool to have the new parameter integrated into `nextflow_schema.json`. The output should look as follows.

    ```
    gitpod /workspace/basic_training/nf-core-demotest (master) $ nf-core schema build

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.13.1 - https://nf-co.re

    INFO [✓] Default parameters match schema validation
    INFO [✓] Pipeline schema looks valid (found 32 params)
    ✨ Found 'params.test' in the pipeline config, but not in the schema. Add to pipeline schema? [y/n]: y
    ```

    Select y on the final prompt to launch a web browser to edit your schema graphically.

    </details>

7.  Lint your pipeline
    <details>
       <summary>solution 7</summary>

    ```bash
    nf-core lint
    ```

       </details>

8.  Run the pipeline and inspect the results
    <details>
       <summary>solution 8</summary>

    To run the pipeline, be aware that we now need to specify a file containing the adapters. As such, we create a new file called "adapterlist.txt" and add the adapter sequence "[WE NEED AN ADAPTER SEQUENCE HERE]" to it. Then we can run the pipeline as follows:

    ```bash
    nextflow run nf-core-demotest/ -profile test,docker --outdir test_results --adapterlist /path/to/adapterlist.txt

    ```

       </details>

9.  Commit the changes
    <details>
       <summary>solution 9</summary>

    ```bash
    git add .
    git commit -m "add adapterremoval module"
    ```

       </details>

:::

### Adding a local module

If there is no nf-core module available for the software you want to include, the nf-core tools package can also aid in the generation of a local module that is specific for your pipeline. To add a local module run the following:

```

nf-core modules create

```

Open ./modules/local/demo/module.nf and start customising this to your needs whilst working your way through the extensive TODO comments!

### Making a local module for a custom script

To generate a module for a custom script you need to follow the same steps when adding a remote module.
Then, you can supply the command for your script in the `script` block but your script needs to be present
and _executable_ in the `bin`
folder of the pipeline.
In the nf-core pipelines,
this folder is in the main directory and you can see in [`rnaseq`](https://github.com/nf-core/rnaseq).
Let's look at an publicly available example in this pipeline,
for instance [`tximport.r`](https://github.com/nf-core/rnaseq/blob/master/bin/tximport.r).
This is an Rscript present in the [`bin`](https://github.com/nf-core/rnaseq/tree/master/bin) of the pipeline.
We can find the module that runs this script in
[`modules/local/tximport`](https://github.com/nf-core/rnaseq/blob/master/modules/local/tximport/main.nf).
As we can see the script is being called in the `script` block, note that `tximport.r` is
being executed as if it was called from the command line and therefore needs to be _executable_.

<blockquote style="border-left: 4px solid #F0AD4E; background-color: #FFF3CD; padding: 10px;">

<h4 style="margin-top: 0;">TL;DR</h4>

1. Write your script on any language (python, bash, R,
   ruby). E.g. `maf2bed.py`
2. If not there yet, move your script to `bin` folder of
   the pipeline and make it
   executable (`chmod +x <filename>`)
3. Create a module with a single process to call your script from within the workflow. E.g. `./modules/local/convert_maf2bed/main.nf`
4. Include your new module in your workflow with the command `include {CONVERT_MAF2BED} from './modules/local/convert_maf2bed/main'` that is written before the workflow call.
</blockquote>

_Tip: Try to follow best practices when writing a script for
reproducibility and maintenance purposes: add the
shebang (e.g. `#!/usr/bin/env python`), and a header
with description and type of license._

### 1. Write your script

Let's create a simple custom script that converts a MAF file to a BED file called `maf2bed.py` and place it in the bin directory of our nf-core-testpipeline::

```

#!/usr/bin/env python
"""bash title="maf2bed.py"
Author: Raquel Manzano - @RaqManzano
Script: Convert MAF to BED format keeping ref and alt info
License: MIT
"""
import argparse
import pandas as pd

def argparser():
parser = argparse.ArgumentParser(description="")
parser.add_argument("-maf", "--mafin", help="MAF input file", required=True)
parser.add_argument("-bed", "--bedout", help="BED input file", required=True)
parser.add_argument(
"--extra", help="Extra columns to keep (space separated list)", nargs="+", required=False, default=[]
)
return parser.parse_args()

def maf2bed(maf_file, bed_file, extra):
maf = pd.read_csv(maf_file, sep="\t", comment="#")
bed = maf[["Chromosome", "Start_Position", "End_Position"] + extra]
bed.to_csv(bed_file, sep="\t", index=False, header=False)

def main():
args = argparser()
maf2bed(maf_file=args.mafin, bed_file=args.bedout, extra=args.extra)

if **name** == "**main**":
main()

```

### 2. Make sure your script is in the right folder

Now, let's move it to the correct directory and make sure it is executable:

```bash
mv maf2bed.py /path/where/pipeline/is/bin/.
chmod +x /path/where/pipeline/is/bin/maf2bed.py
```

### 3. Create your custom module

Then, let's write our module. We will call the process
"CONVERT_MAF2BED" and add any tags or/and labels that
are appropriate (this is optional) and directives (via
conda and/or container) for
the definition of dependencies.

<blockquote style="border-left: 4px solid #F0AD4E; background-color: #FFF3CD; padding: 10px;">

<h4 style="margin-top: 0;">Some additional infos that might be of interest</h4>

<details>
<summary><span style="color: forestgreen; font-weight: bold;">More info on labels</span></summary>
A `label` will
annotate the processes with a reusable identifier of your
choice that can be used for configuring. E.g. we use the
`label` 'process_single', this looks as follows:

```

withLabel:process_single {
cpus = { check_max( 1 _ task.attempt, 'cpus' ) }
memory = { check_max( 1.GB _ task.attempt, 'memory') }
time = { check_max( 1.h \* task.attempt, 'time' ) }
}

```

</details>

<details>
<summary><span style="color: forestgreen; font-weight: bold;">More info on tags</span></summary>

A `tag` is simple a user provided identifier associated to
the task. In our process example, the input is a tuple
comprising a hash of metadata for the maf file called
`meta` and the path to the `maf` file. It may look
similar to: `[[id:'123', data_type:'maf'],
/path/to/file/example.maf]`. Hence, when nextflow makes
the call and `$meta.id` is `123` name of the job
will be "CONVERT_MAF2BED(123)". If `meta` does not have
`id` in its hash, then this will be literally `null`.

</details>

<details>
<summary><span style="color: forestgreen; font-weight: bold;">More info on conda/container directives</span></summary>

The `conda` directive allows for the definition of the
process dependencies using the [Conda package manager](https://docs.conda.io/en/latest/). Nextflow automatically sets up an environment for the given package names listed by in the conda directive. For example:

```

process foo {
conda 'bwa=0.7.15'

'''
your_command --here
'''
}

```

Multiple packages can be specified separating them with a blank space e.g. `bwa=0.7.15 samtools=1.15.1`. The name of the channel from where a specific package needs to be downloaded can be specified using the usual Conda notation i.e. prefixing the package with the channel name as shown here `bioconda::bwa=0.7.15`

```

process foo {
conda 'bioconda::bwa=0.7.15 bioconda::samtools=1.15.1'

'''
your_bwa_cmd --here
your_samtools_cmd --here
'''
}

```

Similarly, we can apply the `container` directive to execute the process script in a [Docker](http://docker.io/) or [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) container. When running Docker, it requires the Docker daemon to be running in machine where the pipeline is executed, i.e. the local machine when using the local executor or the cluster nodes when the pipeline is deployed through a grid executor.

```

process foo {
conda 'bioconda::bwa=0.7.15 bioconda::samtools=1.15.1'
container 'dockerbox:tag'

'''
your_bwa_cmd --here
your_samtools_cmd --here
'''
}

```

Additionally, the `container` directive allows for a more sophisticated choice of container and if it Docker or Singularity depending on the users choice of container engine. This practice is quite common on official nf-core modules.

```

process foo {
conda "bioconda::fastqc=0.11.9"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
'biocontainers/fastqc:0.11.9--0' }"

'''
your_fastqc_command --here
'''
}

```

</details>

</blockquote>

Since `maf2bed.py` is in the `bin` directory we can directory call it in the script block of our new module `CONVERT_MAF2BED`. You only have to be careful with how you call variables (some explanations on when to use `${variable}` vs. `$variable`):
A process may contain any of the following definition blocks: directives, inputs, outputs, when clause, and the process script. Here is how we write it:

```
process CONVERT_MAF2BED {
// HEADER
tag "$meta.id"
    label 'process_single'
    // DEPENDENCIES DIRECTIVES
    conda "anaconda::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
'quay.io/biocontainers/pandas:1.4.3' }"

// INPUT BLOCK
input:
tuple val(meta), path(maf)

// OUTPUT BLOCK
output:
tuple val(meta), path('\*.bed') , emit: bed
path "versions.yml" , emit: versions

// WHEN CLAUSE
when:
task.ext.when == null || task.ext.when

// SCRIPT BLOCK
script: // This script is bundled with the pipeline in bin
def args = task.ext.args ?: ''
def prefix = task.ext.prefix ?: "${meta.id}"

"""
maf2bed.py --mafin $maf --bedout ${prefix}.bed
"""
}
```

More on nextflow's process components in the [docs](https://www.nextflow.io/docs/latest/process.html).

### Include your module in the workflow

In general, we will call out nextflow module `main.nf` and save it in the `modules` folder under another folder called `conver_maf2bed`. If you believe your custom script could be useful for others and it is potentially reusable or calling a tool that is not yet present in nf-core modules you can start the process of making it official adding a `meta.yml` [explained above](#adding-modules-to-a-pipeline). In the `meta.yml` The overall tree for the pipeline skeleton will look as follows:

```

pipeline/
├── bin/
│ └── maf2bed.py
├── modules/
│ ├── local/
│ │ └── convert_maf2bed/
│ │ ├── main.nf
│ │ └── meta.yml
│ └── nf-core/
├── config/
│ ├── base.config
│ └── modules.config
...

```

To use our custom module located in `./modules/local/convert_maf2bed` within our workflow, we use a module inclusions command as follows (this has to be done before we invoke our workflow):

```bash title="workflows/demotest.nf"
include { CONVERT_MAF2BED } from './modules/local/convert_maf2bed/main'
workflow {
input_data = [[id:123, data_type='maf'], /path/to/maf/example.maf]
CONVERT_MAF2BED(input_data)
}
```

:::tip{title="Exercise 6 - Adding a custom module"}
In the directory `exercise_6` you will find the custom script `print_hello.py`, which will be used for this and the next exercise.

1.  Create a local module that runs the `print_hello.py` script
2.  Add the module to your main workflow
3.  Run the pipeline
4.  Lint the pipeline
5.  Commit your changes
    <details>
    <summary>solution 1</summary>

    ```

    ```

      </details>

:::

### Further reading and additional notes

#### What happens in I want to use containers but there is no image created with the packages I need?

No worries, this can be done fairly easy thanks to [BioContainers](https://biocontainers-edu.readthedocs.io/en/latest/what_is_biocontainers.html), see instructions [here](https://github.com/BioContainers/multi-package-containers). If you see the combination that you need in the repo, you can also use [this website](https://midnighter.github.io/mulled) to find out the "mulled" name of this container.

### I want to know more about software dependencies!

You are in luck, we have more documentation [here](https://nf-co.re/docs/contributing/modules#software-requirements)

#### I want to know more about modules!

See more info about modules in the nextflow docs [here](https://nf-co.re/docs/contributing/modules#software-requirements.)

## Lint all modules

As well as the pipeline template you can lint individual or all modules with a single command:

```

nf-core modules lint --all

```

## Nextflow Schema

All nf-core pipelines can be run with --help to see usage instructions. We can try this with the demo pipeline that we just created:

```

cd ../
nextflow run nf-core-demo/ --help

```

### Working with Nextflow schema

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

:::tip{title="Exercise 7 - Using nextflow schema to add command line parameters"}

1.  Feed a string to your custom script from exercise 6 from the command line. Use `nf-core schema build` to add the parameter to the `nextflow.config` file.

      </details>

:::

:::note

### Key points

- `nf-core create <pipeline>` creates a pipeline from the nf-core template.
- `nf-core lint` lints the pipeline code for things that must be completed.
- `nf-core modules list local` lists modules currently installed into your pipeline.
- `nf-core modules list remote` lists modules available to install into your pipeline.
- `nf-core modules install <tool/subtool>` installs the tool module into your pipeline.
- `nf-core modules create` creates a module locally to add custom code into your pipeline.
- `nf-core modules lint --all` lints your module code for things that must be completed.
- `nf-core schema build` opens an interface to allow you to describe your pipeline parameters and set default values, and which values are valid.

:::

```

```
