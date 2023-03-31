---
title: DSL2 Modules
subtitle: Guidelines and reference for DSL2 modules
---

If you decide to upload a module to `nf-core/modules` then this will ensure that it will become available to all nf-core pipelines, and to everyone within the Nextflow community! See [`modules/`](https://github.com/nf-core/modules/tree/master/modules) for examples.

See the [dsl2 modules tutorial](tutorials/dsl2_modules_tutorial) for a step by step guide for how to add a module!

## Terminology

The features offered by Nextflow DSL2 can be used in various ways depending on the granularity with which you would like to write pipelines. Please see the listing below for the hierarchy and associated terminology we have decided to use when referring to DSL2 components.

### Module

A `process` that can be used within different pipelines and is as atomic as possible i.e. cannot be split into another module. An example of this would be a module file containing the process definition for a single tool such as `FastQC`. Atomic nf-core module files are available in the [`modules/`](https://github.com/nf-core/modules/tree/master/modules) directory of nf-core/modules along with the required documentation and tests.

### Subworkflow

A chain of multiple modules that offer a higher-level of functionality within the context of a pipeline. For example, a subworkflow to run multiple QC tools with FastQ files as input. Subworkflows should be shipped with the pipeline implementation and if required they should be shared amongst different pipelines directly from there. Shareable nf-core subworkflow files are available in the [`subworkflow/`](https://github.com/nf-core/modules/tree/master/subworkflows) directory of nf-core/modules along with the required documentation and tests.

### Workflow

What DSL1 users would consider an end-to-end pipeline. For example, from one or more inputs to a series of outputs. This can either be implemented using a large monolithic script as with DSL1, or by using a combination of DSL2 modules and subworkflows. nf-core pipelines can have multiple workflows, such as processing different data types for the same ultimate purpose (such as in [nf-core/viralrecon](https://github.com/nf-core/viralrecon/tree/master/workflows))

## Writing a new module reference

See the [dsl2 modules tutorial](tutorials/dsl2_modules_tutorial) for a step by step guide for how to add a module!

### Before you start

Please check that the module you wish to add isn't already on [`nf-core/modules`](https://github.com/nf-core/modules/tree/master/modules):

- Use the [`nf-core modules list`](https://github.com/nf-core/tools#list-modules) command
- Check [open pull requests](https://github.com/nf-core/modules/pulls)
- Search [open issues](https://github.com/nf-core/modules/issues)

If the module doesn't exist on `nf-core/modules`:

- Please create a [new issue](https://github.com/nf-core/modules/issues/new?assignees=&labels=new%20module&template=new_nodule.md&title=new%20module:) before adding it
- Set an appropriate subject for the issue e.g. `new module: fastqc`
- Add yourself to the `Assignees` so we can track who is working on the module

### New module workflow

We have implemented a number of commands in the `nf-core/tools` package to make it incredibly easy for you to create and contribute your own modules to nf-core/modules.

1. Install the latest version of [`nf-core/tools`](https://github.com/nf-core/tools#installation) (`>=2.7`)
2. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)
3. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Conda`](https://conda.io/miniconda.html)
4. Setup up [pre-commit](https://pre-commit.com/) (comes packaged with [`nf-core/tools`](https://github.com/nf-core/tools#installation), watch the [pre-commit bytesize talk](https://www.youtube.com/watch?v=08d6zv6zvdM&t=215) if you want to know more about it) to ensure that your code is linted and formatted correctly before you commit it to the repository

   ```bash
   pre-commit install
   ```

5. [Fork and clone the nf-core/modules repo locally](#uploading-to-nf-coremodules)
6. Set up git on your computer by adding a new git remote of the main nf-core git repo called `upstream`

   ```bash
   git remote add upstream https://github.com/nf-core/modules.git
   ```

   Make a new branch for your module and check it out

   ```bash
   git checkout -b fastqc
   ```

7. Create a module using the [nf-core DSL2 module template](https://github.com/nf-core/tools/blob/master/nf_core/module-template/modules/main.nf):

   ```console
   $ nf-core modules create fastqc --author @joebloggs --label process_low --meta

                                         ,--./,-.
         ___     __   __   __   ___     /,-._.--~\
   |\ | |__  __ /  ` /  \ |__) |__         }  {
   | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                         `._,._,'

    nf-core/tools version 2.3dev0 - https://nf-co.re

    INFO     Using Bioconda package: 'bioconda::fastqc=0.11.9'                                                                                                           create.py:130
    INFO     Using Docker container: 'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1'                                                                                   create.py:190
    INFO     Using Singularity container: 'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--hdfd78af_1'                                                        create.py:191
    INFO     Created / edited following files:                                                                                                                           create.py:269
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

      This file will be used to store general information about the module and author details - the majority of which will already be auto-filled. However, you will need to add a brief description of the files defined in the `input` and `output` section of the main script since these will be unique to each module. We check it's formatting and validity based on a [JSON schema](https://github.com/nf-core/modules/blob/master/.yaml-schema.json) during linting (and in the pre-commit hook).

   3. [`./tests/modules/fastqc/main.nf`](https://github.com/nf-core/modules/blob/master/tests/modules/fastqc/main.nf)

      Every module MUST have a test workflow. This file will define one or more Nextflow `workflow` definitions that will be used to unit test the output files created by the module. By default, one `workflow` definition will be added but please feel free to add as many as possible so we can ensure that the module works on different data types / parameters e.g. separate `workflow` for single-end and paired-end data.

      When writing multiple tests, a common practice is to alias process names to differentiate them between tests. When using an alias, add a suffix to the process name so the CI tests can still find the output in the folder named after the tool, e.g.

      ```groovy
      include { FASTQC as FASTQC_POST } from '../../.....' // Good: Output folder is still 'fastqc'
      include { FASTQC as POST_FQC    } from '../../.....' // Bad: Generates problems with CI tests - Output folder is 'post'
      ```

      Minimal test data required for your module may already exist within the [nf-core/modules repository](https://github.com/nf-core/modules/blob/master/tests/config/test_data.config), in which case you may just have to change a couple of paths in this file - see the [Test data](#test-data) section for more info and guidelines for adding new standardised data if required.

   4. [`./tests/modules/fastqc/nextflow.config`](https://github.com/nf-core/modules/blob/master/tests/modules/amps/nextflow.config)

      Some modules MAY require additional parameters added to the test command to successfully run. These can be specified with an `ext.args` variable within the process scope of the `nextflow.config` file that exists alongside the test files themselves (and is automatically loaded when the test workflow `main.nf` is executed).

   5. [`./tests/modules/fastqc/test.yml`](https://github.com/nf-core/modules/blob/master/tests/modules/fastqc/test.yml)

      This file will contain all of the details required to unit test the main script in the point above using [pytest-workflow](https://pytest-workflow.readthedocs.io/). If possible, any outputs produced by the test workflow(s) MUST be included and listed in this file along with an appropriate check e.g. md5sum. The different test options are listed in the [pytest-workflow docs](https://pytest-workflow.readthedocs.io/en/stable/#test-options).

      As highlighted in the next point, we have added a command to make it much easier to test the workflow(s) defined for the module and to automatically create the `test.yml` with the md5sum hashes for all of the outputs generated by the module.

      `md5sum` checks are the preferable choice of test to determine file changes, however, this may not be possible for all outputs generated by some tools e.g. if they include time stamps or command-related headers. Please do your best to avoid just checking for the file being present e.g. it may still be possible to check that the file contains the appropriate text snippets.

8. Create a yaml file containing information required for module unit testing

   ```console
   $ nf-core modules create-test-yml

                                         ,--./,-.
         ___     __   __   __   ___     /,-._.--~\
   |\ | |__  __ /  ` /  \ |__) |__         }  {
   | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                         `._,._,'

    nf-core/tools version 2.3dev0 - https://nf-co.re


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

9. Check that the new module you've added follows the [new module guidelines](#new-module-guidelines-and-pr-review-checklist)

10. Lint the module locally to check that it adheres to nf-core guidelines before submission

    ```console
    $ nf-core modules lint fastqc --dir .

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

     nf-core/tools version 2.3dev0 - https://nf-co.re

     INFO     Linting modules repo: .                                                __init__.py:15
     INFO     Linting module: fastqc                                                 __init__.py:163

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

11. Once ready, the code can be pushed and a pull request (PR) created

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

   - Run the test with the helper tool `nf-core modules test` from the modules directory.

     ```console
     $ cd /path/to/git/clone/of/nf-core/modules/
     $ nf-core modules test fastqc
                                               ,--./,-.
               ___     __   __   __   ___     /,-._.--~\
         |\ | |__  __ /  ` /  \ |__) |__         }  {
         | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                               `._,._,'

         nf-core/tools version 2.4

     ? Choose software profile Docker
     INFO     Setting environment variable '$PROFILE' to 'docker'
     INFO     Running pytest for module 'fastqc'

     ========================================== test session starts ==========================================
     platform darwin -- Python 3.9.12, pytest-7.1.2, pluggy-1.0.0
     rootdir: ~/modules, configfile: pytest.ini
     plugins: workflow-1.6.0
     collecting ...
     collected 761 items

     fastqc single-end:
             command:   nextflow run ./tests/modules/fastqc/ -entry test_fastqc_single_end -c ./tests/config/nextflow.config -c ./tests/modules/fastqc/nextflow.config -c ./tests/modules/fastqc/nextflow.config
             directory: /var/folders/lt/b3cs9y610fg_13q14dckwcvm0000gn/T/pytest_workflow_ahvulf1v/fastqc_single-end
             stdout:    /var/folders/lt/b3cs9y610fg_13q14dckwcvm0000gn/T/pytest_workflow_ahvulf1v/fastqc_single-end/log.out
             stderr:    /var/folders/lt/b3cs9y610fg_13q14dckwcvm0000gn/T/pytest_workflow_ahvulf1v/fastqc_single-end/log.err
     'fastqc single-end' done.

     fastqc paired-end:
             command:   nextflow run ./tests/modules/fastqc/ -entry test_fastqc_paired_end -c ./tests/config/nextflow.config -c ./tests/modules/fastqc/nextflow.config -c ./tests/modules/fastqc/nextflow.config
             directory: /var/folders/lt/b3cs9y610fg_13q14dckwcvm0000gn/T/pytest_workflow_ahvulf1v/fastqc_paired-end
             stdout:    /var/folders/lt/b3cs9y610fg_13q14dckwcvm0000gn/T/pytest_workflow_ahvulf1v/fastqc_paired-end/log.out
             stderr:    /var/folders/lt/b3cs9y610fg_13q14dckwcvm0000gn/T/pytest_workflow_ahvulf1v/fastqc_paired-end/log.err
     'fastqc paired-end' done.

     tests/test_versions_yml.py ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss [ 17%]
     ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss [ 38%]
     ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss [ 59%]
     sssssssssssssssssssssssssssssssssssss..sssssssssssssssssssssssssssssssssssssssssssssssssssssssssss [ 80%]
     ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss     [ 98%]
     tests/modules/fastqc/test.yml ........
     Keeping temporary directories and logs. Use '--kwd' or '--keep-workflow-wd' to disable this behaviour.
     ============================= 10 passed, 751 skipped, 479 warnings in 50.76s =============================
     ```

   - See [docs on running pytest-workflow](https://pytest-workflow.readthedocs.io/en/stable/#running-pytest-workflow) for more info.

> 🛈 For docker/singularity, setting the environment variable `TMPDIR=~` is an example of a location the containers can mount (you can change this as you prefer). If you get test failures such as with Nextflow errors that end in `work doesn't exist in container`, check your container can mount your `TMPDIR`.
>
> :warning: if you have a module named `build` this can conflict with some pytest internal behaviour. This results in no tests being run (i.e. receiving a message of `collected 0 items`). In this case rename the `tests/<module>/build` directory to `tests/<module>/build_test`, and update the corresponding `test.yml` accordingly. An example can be seen with the [`bowtie2/build` module tests](https://github.com/nf-core/modules/tree/master/tests/modules/bowtie2/build_test).

### Uploading to `nf-core/modules`

[Fork](https://help.github.com/articles/fork-a-repo/) the `nf-core/modules` repository to your own GitHub account. Within the local clone of your fork add the module file to the `modules/` directory. Please try and keep PRs as atomic as possible to aid the reviewing process - ideally, one module addition/update per PR.

Commit and push these changes to your local clone on GitHub, and then [create a pull request](https://help.github.com/articles/creating-a-pull-request-from-a-fork/) on the `nf-core/modules` GitHub repo with the appropriate information.

When you are happy with your pull request, please <span class="x x-first x-last">select </span>the `Ready for Review` label on the GitHub PR tab, and providing that everything adheres to nf-core guidelines we will endeavour to approve your pull request as soon as possible. We also recommend to request reviews from the `nf-core/modules-team`<span class="x x-first x-last"> so </span>a core team of volunteers <span class="x x-first x-last">can try</span> to <span class="x x-first x-last">review </span>your <span class="x x-first x-last">PR</span> as fast as possible.

Once you<span class="x x-first x-last"> are </span>familiar with the module submission process, please consider joining the<span class="x x-first x-last"> reviewing</span> team by asking on the `#modules` slack channel.

### Talks

> ⚠️ these may include references to an older syntax, however the general idea remains the same

<div class="ratio ratio-16x9">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/xuNYATGFuw4" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
</div>

<div class="ratio ratio-16x9">
     <iframe src="https://widgets.figshare.com/articles/16825369/embed?show_title=1" width="568" height="351" allowfullscreen frameborder="0"></iframe>
</div>

## New module guidelines and PR review checklist

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

### General

1. All command-line tool non-file arguments MUST be provided as a string via the `$task.ext.args` variable, unless an argument is needed to modify the command (for example `lib_type` in [salmon/quant](https://github.com/nf-core/modules/blob/master/modules/nf-core/salmon/quant/main.nf)). The value of `task.ext.args` is supplied from the `modules.config` file by assigning a string value to `ext.args`.
   Mandatory command line arguments MUST be specified in long form where possible.

   `<module>.nf`:

   ```nextflow
   script:
   def args = task.ext.args ?: ''
   def prefix = task.ext.prefix ?: "${meta.id}"
   """
   fastqc \\
       $args \\
        <...>
   """
   ```

   `modules.config`:

   ```nextflow
   process {
       withName: <module> {
           ext.args = [                                                          // Assign either a string, or closure which returns a string
               '--quiet',
               params.fastqc_kmer_size ? "-k ${params.fastqc_kmer_size}" : ''    // Parameter dependent values can be provided like so
           ].join(' ')                                                           // Join converts the list here to a string.
           ext.prefix = { "${meta.id}" }                                         // A closure can be used to access variables defined in the script
       }
   }
   ```

   > ⚠️ Exceptions to non-file mandatory arguments may be acceptable in rare cases, however you must consult the community on Slack (#modules) in these cases.

2. Software that can be piped together SHOULD be added to separate module files
   unless there is a run-time, storage advantage in implementing in this way. For example,
   using a combination of `bwa` and `samtools` to output a BAM file instead of a SAM file:

   ```bash
   bwa mem | samtools view -B -T ref.fasta
   ```

3. Where applicable, the usage and generation of compressed files SHOULD be enforced as input and output, respectively:

   - `*.fastq.gz` and NOT `*.fastq`
   - `*.bam` and NOT `*.sam`

   If a tool does not support compressed input or output natively, we RECOMMEND passing the
   uncompressed data via unix pipes, such that it never gets written to disk, e.g.

   ```bash
   gzip -cdf $input | tool | gzip > $output
   ```

   The `-f` option makes `gzip` auto-detect if the input is compressed or not.

   If a tool cannot read from STDIN, or has multiple input files, it is possible to use
   named pipes:

   ```bash
   mkfifo input1_uncompressed input2_uncompressed
   gzip -cdf $input1 > input1_uncompressed &
   gzip -cdf $input2 > input2_uncompressed &
   tool input1_uncompressed input2_uncompressed > $output
   ```

   Only if a tool reads the input multiple times, it is required to uncompress the
   file before running the tool.

4. Where applicable, each module command MUST emit a file `versions.yml` containing the version number for each tool executed by the module, e.g.

```bash
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    samtools: \$( samtools --version 2>&1 | sed 's/^.*samtools //; s/Using.*\$// )
END_VERSION
```

resulting in, for instance,

```yaml
'FASTQC':
  fastqc: 0.11.9
  samtools: 1.12
```

All reported versions MUST be without a leading `v` or similar (i.e. must start with a numeric character), or for
unversioned software, a Git SHA commit id (40 character hexadecimal string).

We chose a [HEREDOC](https://tldp.org/LDP/abs/html/here-docs.html) over piping into the versions file
line-by-line as we believe the latter makes it easy to accidentally overwrite the file. Moreover, the exit status
of the sub-shells evaluated in within the HEREDOC is ignored, ensuring that a tool's version command does
not erroneously terminate the module.

If the software is unable to output a version number on the command-line then a variable called `VERSION` can be manually
specified to provide this information e.g. [homer/annotatepeaks module](https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf).
Please include the accompanying comments above the software packing directives and beside the version string.

```nextflow
process TOOL {
    ...

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::tool=0.9.1:" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tool:0.9.1--pl526hc9558a2_3' :
        'quay.io/biocontainers/tool:0.9.1--pl526hc9558a2_3' }"

    ...

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.9.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    ...

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool: $VERSION
    END_VERSIONS
    """

}
```

If the HEREDOC cannot be used because the script is not bash, the versions.yml must be written directly e.g. [ascat module](https://github.com/nf-core/modules/blob/master/modules/nf-core/ascat/main.nf).

5. The process definition MUST NOT change the `when` statement. `when` conditions can instead be supplied using the `process.ext.when` directive in a configuration file.

```groovy
process {
    withName: 'FOO' {
        ext.when = !params.skip_module
    }
    withName: 'BAR' {
        ext.when = { meta.single_end }
    }
}
```

### Naming conventions

1. The directory structure for the module name must be all lowercase e.g. [`modules/nf-core/bwa/mem/`](https://github.com/nf-core/modules/tree/master/modules/nf-core/bwa/mem/). The name of the software (i.e. `bwa`) and tool (i.e. `mem`) MUST be all one word.

2. The process name in the module file MUST be all uppercase e.g. `process BWA_MEM {`. The name of the software (i.e. `BWA`) and tool (i.e. `MEM`) MUST be all one word separated by an underscore.

3. All parameter names MUST follow the `snake_case` convention.

4. All function names MUST follow the `camelCase` convention.

5. Channel names MUST follow `snake_case` convention and be all lower case.

6. Output file (and/or directory) names SHOULD just consist of only `${prefix}` and the file-format suffix (e.g. `${prefix}.fq.gz` or `${prefix}.bam`).

   - This is primarily for re-usability so that other developers have complete flexibility to name their output files however they wish when using the same module.
   - As a result of using this syntax, if the module has the same named inputs and outputs then you can add a line in the `script` section like below (another example [here](https://github.com/nf-core/modules/blob/e20e57f90b6787ac9a010a980cf6ea98bd990046/modules/lima/main.nf#L37)) which will raise an error asking the developer to change the `args.prefix` variable to rename the output files so they don't clash.

     ```nextflow
     script:
     if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
     ```

### Input/output options

1. Input channel declarations MUST be defined for all _possible_ input files (i.e. both required and optional files).

   - Directly associated auxiliary files to an input file MAY be defined within the same input channel alongside the main input channel (e.g. [BAM and BAI](https://github.com/nf-core/modules/blob/e937c7950af70930d1f34bb961403d9d2aa81c7d/modules/samtools/flagstat/main.nf#L22)).
   - Other generic auxiliary files used across different input files (e.g. common reference sequences) MAY be defined using a dedicated input channel (e.g. [reference files](https://github.com/nf-core/modules/blob/3cabc95d0ed8a5a4e07b8f9b1d1f7ff9a70f61e1/modules/bwa/mem/main.nf#L21-L23)).

2. Named file extensions MUST be emitted for ALL output channels e.g. `path "*.txt", emit: txt`.

3. Optional inputs are not currently supported by Nextflow. However, passing an empty list (`[]`) instead of a file as a module parameter can be used to work around this issue.

For example, having a module (`MY_MODULE`) that can take a `cram` channel and an optional `fasta` channel as input, can be used in the following ways:

```nextflow
MY_MODULE(cram, [])     // fasta is optional, the module will run without the fasta present
MY_MODULE(cram, fasta)  // execution of the module will need an element in the fasta channel
```


4. Optional outputs SHOULD be marked as optional:

   ```nextflow
   tuple val(meta), path('*.tab'), emit: tab,  optional: true
   ```

5. Each output file SHOULD be emitted in its own channel (and no more than one), along with the `meta` map if provided ( the exception is the versions.yml ).

### Documentation

1. Each module MUST have a `meta.yaml` in the same directory as the `main.nf` of the module itself.

2. Keywords SHOULD be sufficient to make the module findable through research domain, data types, and tool function keywords

   - Keywords MUST NOT just be solely of the (sub)tool name

3. Keywords MUST be all lower case

4. Input and Output sections of the `meta.yaml` SHOULD only have entries of input and output channels

5. Input and output tuples MUST be split into separate entries

   - i.e., `meta` should be a separate entry to the `file` it is associated with

6. Input/output types MUST only be of the following categories: `map`, `file`, `directory`, `string`, `integer`, `float`

7. Input/output entries MUST match a corresponding channel in the module itself

   - There should be a one-to-one relationship between the module and the `meta.yaml`

   - Input/output entries MUST NOT combine multiple output channels

8. Input/output descriptions SHOULD be descriptive of the contents of file

   - i.e., not just 'A TSV file'

9. Input/output patterns (if present) MUST follow a [Java glob pattern](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob)

10. Input entries should be marked as Mandatory or Optional

### Module parameters

1. A module file SHOULD only define input and output files as command-line parameters to be executed within the process.

2. All `params` within the module MUST be initialised and used in the local context of the module. In other words, named `params` defined in the parent workflow MUST NOT be assumed to be passed to the module to allow developers to call their parameters whatever they want. In general, it may be more suitable to use additional `input` value channels to cater for such scenarios.

3. If the tool supports multi-threading then you MUST provide the appropriate parameter using the Nextflow `task` variable e.g. `--threads $task.cpus`.

4. Any parameters that need to be evaluated in the context of a particular sample e.g. single-end/paired-end data MUST also be defined within the process.

### Resource requirements

1. An appropriate resource `label` MUST be provided for the module as listed in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/conf/base.config#L29-L46) e.g. `process_single`, `process_low`, `process_medium` or `process_high`.

2. If the tool supports multi-threading then you MUST provide the appropriate parameter using the Nextflow `task` variable e.g. `--threads $task.cpus`. If the tool does not support multi-threading, consider `process_single` unless large amounts of RAM are required.

3. If a module contains _multiple_ tools that supports multi-threading (e.g. [piping output into a samtools command](https://github.com/nf-core/modules/blob/28b023e6f4d0d2745406d9dc6e38006882804e67/modules/bowtie2/align/main.nf#L32-L46)), you MUST assign cpus per tool such that the total number of used CPUs does not exceed `task.cpus`.
   - For example, combining two (or more) tools that both (all) have multi-threading, this can be assigned to the variable [`split_cpus`](https://github.com/nf-core/modules/blob/28b023e6f4d0d2745406d9dc6e38006882804e67/modules/bowtie2/align/main.nf#L32)
   - If one tool is multi-threaded and another uses a single thread, you can specify directly in the command itself e.g. with [`${task.cpus - 1}`](https://github.com/nf-core/modules/blob/6e68c1af9a514bb056c0513ebba6764efd6750fc/modules/bwa/sampe/main.nf#L42-L43)

### Software requirements

[BioContainers](https://biocontainers.pro/#/) is a registry of Docker and Singularity containers automatically created from all of the software packages on [Bioconda](https://bioconda.github.io/). Where possible we will use BioContainers to fetch pre-built software containers and Bioconda to install software using Conda.

1. Software requirements SHOULD be declared within the module file using the Nextflow `container` directive. For single-tool BioContainers, the `nf-core modules create` command will automatically fetch and fill-in the appropriate Conda / Docker / Singularity definitions by parsing the information provided in the first part of the module name:

```nextflow
conda "bioconda::fastqc=0.11.9"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
    'quay.io/biocontainers/fastqc:0.11.9--0' }"
```

2. If the software is available on Conda it MUST also be defined using the Nextflow `conda` directive. Using `bioconda::bwa=0.7.17` as an example, software MUST be pinned to the channel (i.e. `bioconda`) and version (i.e. `0.7.17`). Conda packages MUST not be pinned to a build because they can vary on different platforms.

3. If required, multi-tool containers may also be available on BioContainers e.g. [`bwa` and `samtools`](https://biocontainers.pro/#/tools/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40). You can install and use the [`galaxy-tool-util`](https://anaconda.org/bioconda/galaxy-tool-util) package to search for both single- and multi-tool containers available in Conda, Docker and Singularity format. e.g. to search for Docker (hosted on Quay.io) and Singularity multi-tool containers with both `bowtie` and `samtools` installed you can use the following command:

```console
mulled-search --destination quay singularity --channel bioconda --search bowtie samtools | grep "mulled"
```

> NB: Build information for all tools within a multi-tool container can be obtained in the `/usr/local/conda-meta/history` file within the container.

4. It is also possible for a new multi-tool container to be built and added to BioContainers by submitting a pull request on their [`multi-package-containers`](https://github.com/BioContainers/multi-package-containers) repository.

   - Fork the [multi-package-containers repository](https://github.com/BioContainers/multi-package-containers)
   - Make a change to the `hash.tsv` file in the `combinations` directory see [here](https://github.com/aunderwo/multi-package-containers/blob/master/combinations/hash.tsv#L124) for an example where `pysam=0.16.0.1,biopython=1.78` was added.
   - Commit the code and then make a pull request to the original repo, for [example](https://github.com/BioContainers/multi-package-containers/pull/1661)
   - Once the PR has been accepted a container will get built and you can find it using a search tool in the `galaxy-tool-util conda` package

   ```console
   mulled-search --destination quay singularity conda  --search pysam biopython  | grep "mulled"
   quay         mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f  185a25ca79923df85b58f42deb48f5ac4481e91f-0  docker pull quay.io/biocontainers/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0
   singularity  mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f  185a25ca79923df85b58f42deb48f5ac4481e91f-0  wget https://depot.galaxyproject.org/singularity/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0
   ```

   - You can copy and paste the `mulled-*` path into the relevant Docker and Singularity lines in the Nextflow `process` definition of your module
   - To confirm that this is correct. Spin up a temporary Docker container

     ```console
     docker run --rm -it quay.io/biocontainers/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0  /bin/sh
     ```

     And in the command prompt type

     ```console
     $ grep specs /usr/local/conda-meta/history
     # update specs: ['biopython=1.78', 'pysam=0.16.0.1']
     ```

     The packages should reflect those added to the multi-package-containers repo `hash.tsv` file

5. If the software is not available on Bioconda a `Dockerfile` MUST be provided within the module directory. We will use GitHub Actions to auto-build the containers on the [GitHub Packages registry](https://github.com/features/packages).

### Publishing results

Results are published using Nextflow's native [`publishDir`](https://www.nextflow.io/docs/latest/process.html#publishdir) directive defined in the `modules.config` of a workflow (see [here](https://github.com/nf-core/rnaseq/blob/f7702d5b76a1351e2e7796a5ed3f59943a139fbf/conf/modules.config#L100-L106) for an example.) Results were earlier published using a custom `publishDir` definition, using a Groovy Map defined by `params.modules`.

### Test data config file

If a new test dataset is added to [`tests/config/test_data.config`](https://github.com/nf-core/modules/blob/master/tests/config/test_data.config), check that the config name of the added file(s) follows the scheme of the entire file name with dots replaced with underscores.

For example: the nf-core/test-datasets file `genomics/sarscov2/genome/genome.fasta` labelled as `genome_fasta`, or `genomics/sarscov2/genome/genome.fasta.fai` as `genome_fasta_fai`.

### Using a stub test when required test data is too big

If the module absolute cannot run using tiny test data, there is a possibility to add [stub-run](https://www.nextflow.io/docs/edge/process.html#stub) to the test.yml. In this case it is required to test the module using larger scale data and document how this is done. In addition, an extra script-block labeled `stub:` must be added, and this block must create dummy versions of all expected output files as well as the `versions.yml`. An example is found in the [ascat module](https://github.com/nf-core/modules/blob/master/modules/nf-core/ascat/main.nf). In the `test.yml` the `-stub-run` argument is written as well as the md5sums for each of the files that are added in the stub-block. This causes the stub-code block to be activated when the unit test is run ([example](https://github.com/nf-core/modules/blob/master/tests/modules/nf-core/ascat/test.yml)):

```console
nextflow run tests/modules/<nameofmodule> -entry test_<nameofmodule> -c tests/config/nextflow.config -stub-run
```

## What is the `meta` map?

In nf-core DSL2 pipelines, to add sample-specific information and metadata that is carried throughout the pipeline, we use a meta variable. This avoids the need to create separate channels for each new characteristic.
The meta variable can be passed down to processes as a tuple of the channel containing the actual samples, e.g. FastQ files, and the meta variable.
The `meta map` is a [groovy map](https://www.tutorialspoint.com/groovy/groovy_maps.htm), which is like a python dictionary, as shown below:

```nextflow
[id: 'test', single_end: false]
```

Thus, the information can be accessed within processes and `module.conf` files with the key i.e. `meta.id`

The meta variable can be passed down to processes as a tuple of the channel containing the actual samples, e.g. FastQ files, and the meta variable.

```nextflow
input:
tuple val(meta), path(reads)
```

This pattern doesn't work out of the box with [fromFilePairs](https://www.nextflow.io/docs/edge/channel.html#fromfilepairs)

The difference between the two:

```nextflow
// fromFilePairs
filepairs = [
    SRR493366,
    [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]
]

// meta map
meta_map = [
    [id: 'test', single_end: false], // meta map
    [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]
]
```

As you can see the difference, they are both [groovy lists](https://www.tutorialspoint.com/groovy/groovy_lists.htm).
However, the filepairs just has a `val` that is a string, where as the `meta_map` the first value in the list, is a [groovy map](https://www.tutorialspoint.com/groovy/groovy_maps.htm), which is like a python dictionary.
The only required value is `meta.id` for most of the modules, however, they usually contain fields like `meta.single_end` and `meta.strandedness`

### Common patterns

The `meta map` is generated with [create_fastq_channel function in the input_check subworkflow](https://github.com/nf-core/rnaseq/blob/587c61b441c5e00bd3201317d48b95a82afe6aaa/subworkflows/local/input_check.nf#L23-L45) of most nf-core pipelines. Where the meta information is easily extracted from a samplesheet that contains the input file paths.

#### Generating a `meta map` from file pairs

Sometimes you want to use nf-core modules in small scripts. You don't want to make a samplesheet, or maintain a bunch of validation.
For instance, here's an example script to run fastqc

```nextflow
nextflow.enable.dsl = 2

params.input = "*.fastq.gz"

include { FASTQC } from "./modules/nf-core/modules/fastqc/main"

workflow {
    ch_fastq = Channel.fromFilePairs(params.input, size: -1)
        .map {
            meta, fastq ->
            def fmeta = [:]
            // Set meta.id
            fmeta.id = meta
            // Set meta.single_end
            if (fastq.size() == 1) {
                fmeta.single_end = true
            } else {
                fmeta.single_end = false
            }
            [ fmeta, fastq ]
        }

    FASTQC ( ch_fastq )
}
```

#### Sorting samples by groups

```nextflow
ch_genome_bam.map {
    meta, bam ->
    fmeta = meta.findAll { it.key != 'read_group' }
    fmeta.id = fmeta.id.split('_')[0..-2].join('_')
    [ fmeta, bam ] }
    .groupTuple(by: [0])
    .map { it ->  [ it[0], it[1].flatten() ] }
    .set { ch_sort_bam }
```

### Conclusion

As you can see the `meta map` is a quite flexible way for storing meta data in channels. Feel free to add whatever other key-value pairs your pipeline may need to it. We're looking to add [Custom objects](https://github.com/nf-core/modules/issues/1338) which will lock down the usage a bit more.

## What are the ext properties/keys?

Ext properties or keys are special process directives (See: [ext directive](https://www.nextflow.io/docs/latest/process.html#ext) ) that insert strings into the module scripts. For example, an nf-core module uses the string assigned to `ext.args` ( or `ext.args2`, `ext.args3`, ... ) to insert tool specific options in a module script:

Example:

The configuration

```nextflow
process {
  withName: 'TOOL_SUBTOOL' {
    ext.args = '-T -K'
  }
}
```

inserts the string `-T -K` as options to the module script:

```nextflow
process TOOL_SUBTOOL {
  input:
  tuple val(meta), path(bam)

  output:
  tuple val(meta), path("*.log"), emit: log
  path "versions.yml",            emit: versions

  script:
  def args   = task.ext.args ?: ''          // If ext.args is defined assign it to args
  def prefix = task.ext.prefix ?: meta.id   // If ext.prefix is defined assign it to prefix, otherwise assign meta.id value
  """
  tool subtool $args $bam > ${prefix}.log
  """
}
```

so the script becomes:

```bash
#! /usr/env/bin bash
tool subtool -T -K test.bam > test.log
```

The following table lists the available keys commonly used in nf-core modules.

| Key        | Description                                            |
| ---------- | ------------------------------------------------------ |
| ext.args   | Additional arguments appended to command in module.    |
| ext.args2  | Second set of arguments appended to command in module. |
| ext.args3  | Third set of arguments appended to command in module.  |
| ext.prefix | File name prefix for output files.                     |

To see some more advanced examples of these keys in use see:

- [Set ext.args based on parameter settings](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/modules.config#L222-L226)
- [Set ext.prefix based on task inputs](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/modules.config#L297)
- [Set ext.args based on both parameters and task inputs](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/modules.config#L377-L381)

## Help

For further information or help, don't hesitate to get in touch on [Slack `#modules` channel](https://nfcore.slack.com/channels/modules) (you can join with [this invite](https://nf-co.re/join/slack)).
