---
title: Basic training to create an nf-core pipeline
subtitle: Create a pipeline with `nf-core create`
---

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
It is possible to use nf-core tools for non-nf-core pipelines, but the setup of such pipelines will not be handled in this tutorial.

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

<p class="text-center">
  <a href="/docs/contributing/nf_core_basic_course/gitpod_environment/" class="btn btn-lg btn-success" style="font-size: 14px">
    < go to Chapter 1
  </a>
  <a href="/docs/contributing/nf_core_basic_course/template_walk_through/" class="btn btn-lg btn-success" style="font-size: 14px">
    go to Chapter 3 >
  </a>
</p>
