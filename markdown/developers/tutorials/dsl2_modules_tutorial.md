---
title: "Tutorial: Create a DSL2 Module"
subtitle: Creating a new module for the nf-core modules repository.
---

In this tutorial we will see how to create a new module for the nf-core modules repository. As an example, we will create a module to execute the FastqToBam function of the [FGBIO](http://fulcrumgenomics.github.io/fgbio/) suite of tools.

## Introduction

If you create a new module with the goal of contributing the code to *nf-core*, we recommend to familiarise with the community guidelines and use *nf-core tools* as explained below.

### Module guidelines

The nf-core community has agreed on a minimal set of [guidelines](https://nf-co.re/developers/modules#guidelines), intended to make module most suitable for general use, i.e. to be shared across a wide variety of community workflows.

### nf-core tools

Using [nf-core tools](https://nf-co.re/tools) is the best way to adhere to the guidelines, without worrying too much and writing things from scratch.
On the website you can find more details about [installation](https://nf-co.re/tools#installation), and all functionalities for [modules](https://nf-co.re/tools#modules).

## Before you start

Please check that the module you wish to add isn't already on [`nf-core/modules`](https://github.com/nf-core/modules/tree/master/modules):

- Use the [`nf-core modules list`](https://github.com/nf-core/tools#list-modules) command
- Check [open pull requests](https://github.com/nf-core/modules/pulls)
- Search [open issues](https://github.com/nf-core/modules/issues)

If the module doesn't exist on `nf-core/modules`:

- Please create a [new issue](https://github.com/nf-core/modules/issues/new?assignees=&labels=new%20module&template=new_nodule.md&title=new%20module:) before adding it
- Set an appropriate subject for the issue e.g. `new module: fastqc`
- Add yourself to the `Assignees` so we can track who is working on the module

### Test data

Even before beginning the development of a module, you should identify a small dataset you can use to test its functionality. Ideally, the dataset is existing already and can be shared with other test workflows for other modules.

> :recycle: This is in active development, keep an eye for available test data [here](https://github.com/nf-core/test-datasets/tree/modules/data) and how to access them using a config file (see this [change](https://github.com/nf-core/modules/pull/365)).

---

## Writing a new module

Check out [Tutorial: Create a DSL2 Module](https://nf-co.re/developers/tutorials/dsl2_modules_tutorial) for how to create a new DSL2 Module

### New module workflow

We have implemented a number of commands in the `nf-core/tools` package to make it incredibly easy for you to create and contribute your own modules to nf-core/modules.

1. Install the latest version of [`nf-core/tools`](https://github.com/nf-core/tools#installation) (`>=2.1`)

    > ⚠️ [2021-11-26] please use `dev` version of tools

2. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)
3. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Conda`](https://conda.io/miniconda.html)
4. [Fork and clone the nf-core/modules repo locally](#uploading-to-nf-coremodules)
5. Set up git on your computer by adding a new git remote of the main nf-core git repo called `upstream`

   ```bash
   git remote add upstream https://github.com/nf-core/modules.git
   ```

   Make a new branch for your module and check it out

   ```bash
   git checkout -b fastqc
   ```

6. Create a module using the [nf-core DSL2 module template](https://github.com/nf-core/tools/blob/master/nf_core/module-template/modules/main.nf):

    > ⚠️ [2021-11-26] please use `dev` version of tools

    ```console
    $ nf-core modules create fastqc --author @joebloggs --label process_low --meta

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.2.dev0

    INFO     Using Bioconda package: 'bioconda::fastqc=0.11.9'                      create.py:130
    INFO     Using Docker / Singularity container with tag: 'fastqc:0.11.9--0'      create.py:140
    INFO     Created / edited following files:                                      create.py:218
                ./modules/fastqc/main.nf
                ./modules/fastqc/meta.yml
                ./tests/modules/fastqc/main.nf
                ./tests/modules/fastqc/test.yml
                ./tests/modules/fastqc/nextflow.config
                ./tests/config/pytest_modules.yml
    ```

    All of the files required to add the module to `nf-core/modules` will be created/edited in the appropriate places. There are at most 5 files to modify:

    1. [`./modules/fastqc/main.nf`](https://github.com/nf-core/modules/blob/master/modules/fastqc/main.nf)

        This is the main script containing the `process` definition for the module. You will see an extensive number of `TODO` statements to help guide you to fill in the appropriate sections and to ensure that you adhere to the guidelines we have set for module submissions.

    2. [`./modules/fastqc/meta.yml`](https://github.com/nf-core/modules/blob/master/modules/fastqc/meta.yml)

        This file will be used to store general information about the module and author details - the majority of which will already be auto-filled. However, you will need to add a brief description of the files defined in the `input` and `output` section of the main script since these will be unique to each module.

    3. [`./tests/modules/fastqc/main.nf`](https://github.com/nf-core/modules/blob/master/tests/modules/fastqc/main.nf)

        Every module MUST have a test workflow. This file will define one or more Nextflow `workflow` definitions that will be used to unit test the output files created by the module. By default, one `workflow` definition will be added but please feel free to add as many as possible so we can ensure that the module works on different data types / parameters e.g. separate `workflow` for single-end and paired-end data.

        Minimal test data required for your module may already exist within the [nf-core/modules repository](https://github.com/nf-core/modules/blob/master/tests/config/test_data.config), in which case you may just have to change a couple of paths in this file - see the [Test data](#test-data) section for more info and guidelines for adding new standardised data if required.

    4. [`./tests/modules/fastqc/nextflow.config`](https://github.com/nf-core/modules/blob/master/tests/modules/amps/nextflow.config)

        Some modules MAY require additional parameters added to the test command to successfully run. These can be specified with an `ext.args` variable within the process scope of the `nextflow.config` file that exists alongside the test files themselves (and is automatically loaded when the test workflow `main.nf` is executed).

    5. [`./tests/modules/fastqc/test.yml`](https://github.com/nf-core/modules/blob/master/tests/modules/fastqc/test.yml)

        This file will contain all of the details required to unit test the main script in the point above using [pytest-workflow](https://pytest-workflow.readthedocs.io/). If possible, any outputs produced by the test workflow(s) MUST be included and listed in this file along with an appropriate check e.g. md5sum. The different test options are listed in the [pytest-workflow docs](https://pytest-workflow.readthedocs.io/en/stable/#test-options).

        As highlighted in the next point, we have added a command to make it much easier to test the workflow(s) defined for the module and to automatically create the `test.yml` with the md5sum hashes for all of the outputs generated by the module.

        `md5sum` checks are the preferable choice of test to determine file changes, however, this may not be possible for all outputs generated by some tools e.g. if they include time stamps or command-related headers. Please do your best to avoid just checking for the file being present e.g. it may still be possible to check that the file contains the appropriate text snippets.

7. Create a yaml file containing information required for module unit testing

    ```console
    $ nf-core modules create-test-yml

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.2.dev0


    INFO     Press enter to use default values (shown in brackets) or type your own responses                                             test_yml_builder.py:51
    ? Tool name: fastqc
    Test YAML output path (- for stdout) (tests/modules/fastqc/test.yml):
    INFO     Looking for test workflow entry points: 'tests/modules/fastqc/main.nf'                                                      test_yml_builder.py:116
    INFO     Building test meta for entry point 'test_fastqc_single_end'                                                                  test_yml_builder.py:150
    Test name (fastqc test_fastqc_single_end):
    Test command (nextflow run tests/modules/fastqc -entry test_fastqc_single_end -c tests/config/nextflow.config):
    Test tags (comma separated) (fastqc,fastqc_single_end):
    Test output folder with results (leave blank to run test):
    ? Choose software profile Singularity
    INFO     Setting env var '$PROFILE' to 'singularity'                                                                                  test_yml_builder.py:258
    INFO     Running 'fastqc' test with command:                                                                                          test_yml_builder.py:263
             nextflow run tests/modules/fastqc -entry test_fastqc_single_end -c tests/config/nextflow.config --outdir /tmp/tmpgbneftf5
    INFO     Test workflow finished!                                                                                                      test_yml_builder.py:276
    INFO     Writing to 'tests/modules/fastqc/test.yml'                                                                                  test_yml_builder.py:293
    ```

    > NB: See docs for [running tests manually](#running-tests-manually) if you would like to run the tests manually.

8. Check that the new module you've added follows the [new module guidelines](#new-module-guidelines-and-pr-review-checklist)

9. Lint the module locally to check that it adheres to nf-core guidelines before submission

    ```console
    $ nf-core modules lint fastqc --dir .

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.0

    INFO     Linting modules repo: .                                                lint.py:102
    INFO     Linting module: fastqc                                                 lint.py:106

    ╭────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
    │ [!] 3 Test Warnings                                                                                            │
    ╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
    ╭──────────────┬──────────────────────────────────────────────────────────────┬──────────────────────────────────╮
    │ Module name  │ Test message                                                 │ File path                        │
    ├──────────────┼──────────────────────────────────────────────────────────────┼──────────────────────────────────┤
    │ fastqc       │ TODO string in meta.yml: #Add a description of the module... │ modules/nf-core/modules/fastqc/  │
    │ fastqc       │ TODO string in meta.yml: #Add a description and other det... │ modules/nf-core/modules/fastqc/  │
    │ fastqc       │ TODO string in meta.yml: #Add a description of all of the... │ modules/nf-core/modules/fastqc/  │
    ╰──────────────┴──────────────────────────────────────────────────────────────┴──────────────────────────────────╯
    ╭────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
    │ [!] 1 Test Failed                                                                                              │
    ╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
    ╭──────────────┬──────────────────────────────────────────────────────────────┬──────────────────────────────────╮
    │ Module name  │ Test message                                                 │ File path                        │
    ├──────────────┼──────────────────────────────────────────────────────────────┼──────────────────────────────────┤
    │ fastqc       │ 'meta' map not emitted in output channel(s)                  │ modules/nf-core/modules/fastqc/  │
    ╰──────────────┴──────────────────────────────────────────────────────────────┴──────────────────────────────────╯
    ╭──────────────────────╮
    │ LINT RESULTS SUMMARY │
    ├──────────────────────┤
    │ [✔]  38 Tests Passed │
    │ [!]   3 Test Warning │
    │ [✗]   1 Test Failed  │
    ╰──────────────────────╯
    ```

10. Once ready, the code can be pushed and a pull request (PR) created

    On a regular basis you can pull upstream changes into this branch and it is recommended to do so before pushing and creating a pull request - see below. Rather than merging changes directly from upstream the rebase strategy is recommended so that your changes are applied on top of the latest master branch from the nf-core repo. This can be performed as follows

   ```bash
   git pull --rebase upstream master
   ```

   Once you are ready you can push the code and create a PR

   ```bash
   git push -u origin fastqc
   ```

   Once the PR has been accepted you should delete the branch and checkout master again.

   ```bash
   git checkout master
   git branch -d fastqc
   ```

   In case there are commits on the local branch that didn't make it into the PR (usually commits made after the PR), git will warn about this and not delete the branch. If you are sure you want to delete, use the following command

   ```bash
   git branch -D fastqc
   ```

### Test data

In order to test that each module added to `nf-core/modules` is actually working and to be able to track any changes to results files between module updates we have set-up a number of Github Actions CI tests to run each module on a minimal test dataset using Docker, Singularity and Conda.

- All test data for `nf-core/modules` MUST be added to the `modules` branch of [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) and organised by filename extension.

- In order to keep the size of the test data repository as minimal as possible, pre-existing files from [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) MUST be reused if at all possible.

- Test files MUST be kept as tiny as possible.

- If the appropriate test data doesn't exist in the `modules` branch of [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) please contact us on the [nf-core Slack `#modules` channel](https://nfcore.slack.com/channels/modules) (you can join with [this invite](https://nf-co.re/join/slack)) to discuss possible options.

- It may not be possible to add test data for some modules e.g. if the input data is too large or requires a local database. In these scenarios, it is recommended to use the Nextflow [`stub`](https://www.nextflow.io/docs/latest/process.html#stub) feature to test the module. Please refer to the [`gtdbtk/classify`](https://github.com/nf-core/modules/blob/79d38a306bdaf07000e0d6f300684d3ed38c8919/modules/gtdbtk/classifywf/main.nf#L66) module and its corresponding [test script](https://github.com/nf-core/modules/blob/79d38a306bdaf07000e0d6f300684d3ed38c8919/tests/modules/gtdbtk/classifywf/main.nf#L20) to understand how to use this feature for your module development.

### Running tests manually

As outlined in the [nf-core modules create](#nf-core-modules-create) section we have made it quite trivial to create an initial yaml file (via the `nf-core modules create-test-yml` command) containing a listing of all of the module output files and their associated md5sums. However, md5sum checks may not be appropriate for all output files if for example they contain timestamps. This is why it is a good idea to re-run the tests locally with `pytest-workflow` before you create your pull request adding the module. If your files do indeed have timestamps or other issues that prevent you from using the md5sum check, then you can edit the `test.yml` file to instead check that the file contains some specific content or as a last resort, if it exists. The different test options are listed in the [pytest-workflow docs](https://pytest-workflow.readthedocs.io/en/stable/#test-options).

Please follow the steps below to run the tests locally:

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Conda`](https://conda.io/miniconda.html)

3. Install [`pytest-workflow`](https://pytest-workflow.readthedocs.io/en/stable/#installation)

4. Start running your own tests using the appropriate [`tag`](https://github.com/nf-core/modules/blob/20d8250d9f39ddb05dfb437603aaf99b5c0b2b41/tests/modules/fastqc/test.yml) defined in the `test.yml`:

    - Typical command with Docker:

        ```console
        cd /path/to/git/clone/of/nf-core/modules/
        NF_CORE_MODULES_TEST=1 TMPDIR=~ PROFILE=docker pytest --tag fastqc --symlink --keep-workflow-wd
        ```

    - Typical command with Singularity:

        ```console
        cd /path/to/git/clone/of/nf-core/modules/
        NF_CORE_MODULES_TEST=1 TMPDIR=~ PROFILE=singularity pytest --tag fastqc --symlink --keep-workflow-wd
        ```

    - Typical command with Conda:

        ```console
        cd /path/to/git/clone/of/nf-core/modules/
        NF_CORE_MODULES_TEST=1 PROFILE=conda pytest --tag fastqc --symlink --keep-workflow-wd
        ```

    - See [docs on running pytest-workflow](https://pytest-workflow.readthedocs.io/en/stable/#running-pytest-workflow) for more info.

> 🛈 For docker/singularity`TMPDIR=~` is an example of a location the containers can mount (you can change this as you prefer). If you get test failures such as with Nextflow errors that end in `work doesn't exist in container`, check your container can mount your `TMPDIR`.
>
> :warning: if you have a module named `build` this can conflict with some pytest internal behaviour. This results in no tests being run (i.e. recieving a message of `collected 0 items`). In this case rename the `tests/<module>/build` directry to `tests/<module>/build_test`, and update the corresponding `test.yml` accordingly. An example can be seen with the [`bowtie2/build` module tests](https://github.com/nf-core/modules/tree/master/tests/modules/bowtie2/build_test).

### Uploading to `nf-core/modules`

[Fork](https://help.github.com/articles/fork-a-repo/) the `nf-core/modules` repository to your own GitHub account. Within the local clone of your fork add the module file to the `modules/` directory. Please try and keep PRs as atomic as possible to aid the reviewing process - ideally, one module addition/update per PR.

Commit and push these changes to your local clone on GitHub, and then [create a pull request](https://help.github.com/articles/creating-a-pull-request-from-a-fork/) on the `nf-core/modules` GitHub repo with the appropriate information.

When you are happy with your pull request, please <span class="x x-first x-last">select </span>the `Ready for Review` label on the GitHub PR tab, and providing that everything adheres to nf-core guidelines we will endeavour to approve your pull request as soon as possible. We also recommend to request reviews from the `nf-core/modules-team`<span class="x x-first x-last"> so </span>a core team of volunteers <span class="x x-first x-last">can try</span> to <span class="x x-first x-last">review </span>your <span class="x x-first x-last">PR</span> as fast as possible.

Once you<span class="x x-first x-last"> are </span>familiar with the module submission process, please consider joining the<span class="x x-first x-last"> reviewing</span> team by asking on the `#modules` slack channel.

## New module guidelines and PR review checklist

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Other tutorials

> ⚠️ these may include references to an older syntax, however the general idea remains the same

<div class="ratio ratio-16x9">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/xuNYATGFuw4" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
</div>

<div class="ratio ratio-16x9">
     <iframe src="https://widgets.figshare.com/articles/16825369/embed?show_title=1" width="568" height="351" allowfullscreen frameborder="0"></iframe>
</div>
