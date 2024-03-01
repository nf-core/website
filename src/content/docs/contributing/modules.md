---
title: DSL2 Modules
subtitle: Guidelines and reference for DSL2 modules
---

If you decide to upload a module to `nf-core/modules` then this will ensure that it will become available to all nf-core pipelines, and to everyone within the Nextflow community! See [`modules/`](https://github.com/nf-core/modules/tree/master/modules) for examples.

See the [dsl2 modules tutorial](tutorials/dsl2_modules_tutorial) for a step by step guide for how to add a module!

:::warning
DSL1 has reached it's end-of-life as of 2023. As of Nextflow version 22.04.x and 22.10.x it will not be possible to run DSL1 scripts.
:::

## Terminology

A domain-specific language (DSL) is a programming language that is developed for a specific application. Nextflow is based on a DSL, where DSL2 is the latest version. DSL2 allows data analysis pipelines to be scaled and modularised. The features offered by Nextflow DSL2 can be used in various ways depending on the granularity with which you would like to write pipelines. Please see the listing below for the hierarchy and associated terminology we have decided to use when referring to DSL2 components.

### Module

A `process` that can be used within different pipelines and is as atomic as possible i.e. cannot be split into another module. An example of this would be a module file containing the process definition for a single tool such as `FastQC`. Atomic nf-core module files are available in the [`modules/`](https://github.com/nf-core/modules/tree/master/modules) directory of nf-core/modules along with the required documentation and tests.

### Subworkflow

A chain of multiple modules that offer a higher-level of functionality within the context of a pipeline. For example, a subworkflow to run multiple QC tools with FastQ files as input. Subworkflows should be shipped with the pipeline implementation and if required they should be shared amongst different pipelines directly from there. Shareable nf-core subworkflow files are available in the [`subworkflow/`](https://github.com/nf-core/modules/tree/master/subworkflows) directory of nf-core/modules along with the required documentation and tests.

### Workflow

An end-to-end pipeline where one or more inputs produce a series of outputs. This can either be implemented using a large monolithic script or by using a combination of DSL2 modules and subworkflows. nf-core pipelines can have multiple workflows, such as processing different data types for the same ultimate purpose (such as in [nf-core/viralrecon](https://github.com/nf-core/viralrecon/tree/master/workflows))

## Writing a new module reference

See the [dsl2 modules tutorial](tutorials/dsl2_modules_tutorial/) for a step by step guide for how to add a module!

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
3. Install [`nf-test`](https://code.askimed.com/nf-test/installation/)
4. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Conda`](https://conda.io/miniconda.html)
5. Setup up [pre-commit](https://pre-commit.com/) (comes packaged with [`nf-core/tools`](https://github.com/nf-core/tools#installation), watch the [pre-commit bytesize talk](https://www.youtube.com/watch?v=08d6zv6zvdM&t=215) if you want to know more about it) to ensure that your code is linted and formatted correctly before you commit it to the repository

   ```bash
   pre-commit install
   ```

6. [Fork and clone the nf-core/modules repo locally](#uploading-to-nf-coremodules)
7. Set up git on your computer by adding a new git remote of the main nf-core git repo called `upstream`

   ```bash
   git remote add upstream https://github.com/nf-core/modules.git
   ```

   Make a new branch for your module and check it out

   ```bash
   git checkout -b fastqc
   ```

8. Create a module using the [nf-core DSL2 module template](https://github.com/nf-core/tools/blob/master/nf_core/module-template/modules/main.nf):

   ```console
   $ nf-core modules create fastqc --author @joebloggs --label process_low --meta

                                         ,--./,-.
         ___     __   __   __   ___     /,-._.--~\
   |\ | |__  __ /  ` /  \ |__) |__         }  {
   | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                         `._,._,'

    nf-core/tools version 2.11 - https://nf-co.re

    INFO     Using Bioconda package: 'bioconda::fastqc=0.12.1'
    INFO     Using Docker container: 'biocontainers/fastqc:0.12.1--hdfd78af_0'
    INFO     Using Singularity container: 'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0'
    INFO     Created component template: 'fastqc'
    INFO     Created following files:
             modules/nf-core/fastqc/main.nf
             modules/nf-core/fastqc/meta.yml
             modules/nf-core/fastqc/environment.yml
             modules/nf-core/fastqc/tests/tags.yml
             modules/nf-core/fastqc/tests/main.nf.test
   ```

   All of the files required to add the module to `nf-core/modules` will be created/edited in the appropriate places. There are at most 5 files to modify:

   1. [`./modules/fastqc/main.nf`](https://github.com/nf-core/modules/blob/master/modules/nf-core/fastqc/main.nf)

      This is the main script containing the `process` definition for the module. You will see an extensive number of `TODO` statements to help guide you to fill in the appropriate sections and to ensure that you adhere to the guidelines we have set for module submissions.

   2. [`./modules/fastqc/meta.yml`](https://github.com/nf-core/modules/blob/master/modules/nf-core/fastqc/meta.yml)

      This file will be used to store general information about the module and author details - the majority of which will already be auto-filled. However, you will need to add a brief description of the files defined in the `input` and `output` section of the main script since these will be unique to each module. We check it's formatting and validity based on a [JSON schema](https://github.com/nf-core/modules/blob/master/.yaml-schema.json) during linting (and in the pre-commit hook).

   3. [`./modules/nf-core/fastqc/tests/main.nf.test`](https://github.com/nf-core/modules/blob/master/modules/nf-core/fastqc/tests/main.nf.test)

      Every module MUST have a test workflow. This file will define one or more Nextflow `workflow` definitions that will be used to unit test the output files created by the module. By default, one `workflow` definition will be added but please feel free to add as many as possible so we can ensure that the module works on different data types / parameters e.g. separate `workflow` for single-end and paired-end data.

      Minimal test data required for your module may already exist within the [nf-core/modules repository](https://github.com/nf-core/modules/blob/master/tests/config/test_data.config), in which case you may just have to change a couple of paths in this file - see the [Test data](#test-data) section for more info and guidelines for adding new standardised data if required.

      Refer to the section [writing nf-test tests](#writing-nf-test-tests) for more information on how to write nf-tests

9. Create a snapshot of the tests

   ```console
   $ nf-core modules test fastqc

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.11 - https://nf-co.re


    INFO     Generating nf-test snapshot
    ╭────────────────────────────────────────────────────────── nf-test output ──────────────────────────────────────────────────────────╮
    │                                                                                                                                    │
    │ 🚀 nf-test 0.8.2                                                                                                                   │
    │ https://code.askimed.com/nf-test                                                                                                   │
    │ (c) 2021 - 2023 Lukas Forer and Sebastian Schoenherr                                                                               │
    │                                                                                                                                    │
    │ Found 1 files in test directory.                                                                                                   │
    │                                                                                                                                    │
    │ Test Workflow FASTQ_FASTQC_UMITOOLS_FASTP                                                                                          │
    │                                                                                                                                    │
    │   Test [f702dde6] 'sarscov2 paired-end [fastq]' PASSED (46.176s)                                                                   │
    │                                                                                                                                    │
    │ Test Process FASTQC                                                                                                                │
    │                                                                                                                                    │
    │   Test [e588093e] 'Single-Read' PASSED (26.528s)                                                                                   │
    │                                                                                                                                    │
    │                                                                                                                                    │
    │ SUCCESS: Executed 2 tests in 72.706s                                                                                               │
    │                                                                                                                                    │
    ╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
    INFO     Generating nf-test snapshot again to check stability
    ╭────────────────────────────────────────────────────────── nf-test output ──────────────────────────────────────────────────────────╮
    │                                                                                                                                    │
    │ 🚀 nf-test 0.8.2                                                                                                                   │
    │ https://code.askimed.com/nf-test                                                                                                   │
    │ (c) 2021 - 2023 Lukas Forer and Sebastian Schoenherr                                                                               │
    │                                                                                                                                    │
    │ Found 1 files in test directory.                                                                                                   │
    │                                                                                                                                    │
    │ Test Workflow FASTQ_FASTQC_UMITOOLS_FASTP                                                                                          │
    │                                                                                                                                    │
    │   Test [f702dde6] 'sarscov2 paired-end [fastq]' PASSED (39.268s)                                                                   │
    │                                                                                                                                    │
    │ Test Process FASTQC                                                                                                                │
    │                                                                                                                                    │
    │   Test [e588093e] 'Single-Read' PASSED (26.016s)                                                                                   │
    │                                                                                                                                    │
    │                                                                                                                                    │
    │ SUCCESS: Executed 2 tests in 65.286s                                                                                               │
    │                                                                                                                                    │
    ╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
    INFO     All tests passed!
   ```

   :::note
   See the [nf-test docs](https://code.askimed.com/nf-test/) if you would like to run the tests manually.
   :::

10. Check that the new module you've added follows the [new module guidelines](#new-module-guidelines-and-pr-review-checklist)

11. Run [`prettier`](/docs/contributing/code_formatting) on all edited and generated files:

    ```bash
    prettier -w .
    ```

12. Lint the module locally to check that it adheres to nf-core guidelines before submission

    ```console
    $ nf-core modules lint fastqc

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.8 - https://nf-co.re


    INFO     Linting modules repo: '.'
    INFO     Linting module: 'fastqc'

    ╭─ [!] 10 Module Test Warnings ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
    │                                           ╷                                       ╷                                                                                                                             │
    │ Module name                               │ File path                             │ Test message                                                                                                                │
    │╶──────────────────────────────────────────┼───────────────────────────────────────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╴│
    │ fastqc                                    │ modules/nf-core/fastqc/main.nf        │ Conda update: bioconda::fastqc 0.11.9 -> 0.12.1                                                                             │
    │ fastqc                                    │ tests/modules/nf-core/fastqc/test.yml │ TODO string in test.yml: - 'file md5sum was variable, please replace this text with a string found in the file instead '    │
    │ fastqc                                    │ tests/modules/nf-core/fastqc/test.yml │ TODO string in test.yml: - 'file md5sum was variable, please replace this text with a string found in the file instead '    │
    │ fastqc                                    │ tests/modules/nf-core/fastqc/test.yml │ TODO string in test.yml: - 'file md5sum was variable, please replace this text with a string found in the file instead '    │
    │ fastqc                                    │ tests/modules/nf-core/fastqc/test.yml │ TODO string in test.yml: - 'file md5sum was variable, please replace this text with a string found in the file instead '    │
    │ fastqc                                    │ tests/modules/nf-core/fastqc/test.yml │ TODO string in test.yml: - 'file md5sum was variable, please replace this text with a string found in the file instead '    │
    │ fastqc                                    │ tests/modules/nf-core/fastqc/test.yml │ TODO string in test.yml: - 'file md5sum was variable, please replace this text with a string found in the file instead '    │
    │ fastqc                                    │ tests/modules/nf-core/fastqc/test.yml │ TODO string in test.yml: - 'file md5sum was variable, please replace this text with a string found in the file instead '    │
    │ fastqc                                    │ tests/modules/nf-core/fastqc/test.yml │ TODO string in test.yml: - 'file md5sum was variable, please replace this text with a string found in the file instead '    │
    │ fastqc                                    │ tests/modules/nf-core/fastqc/test.yml │ TODO string in test.yml: - 'file md5sum was variable, please replace this text with a string found in the file instead '    │
    │                                           ╵                                       ╵                                                                                                                             │
    ╰─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
    ╭───────────────────────╮
    │ LINT RESULTS SUMMARY  │
    ├───────────────────────┤
    │ [✔]  23 Tests Passed  │
    │ [!]  10 Test Warnings │
    │ [✗]   0 Tests Failed  │
    ╰───────────────────────╯
    ```

13. Once ready, the code can be pushed and a pull request (PR) created

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

### Writing nf-test tests

We recently decided to use nf-test instead of pytest for testing modules. This is because nf-test is more flexible and allows us to test modules in a more realistic way. You can find more information at [nf-test official docs](https://code.askimed.com/nf-test/) and [in this bytesize talk](https://nf-co.re/events/2022/bytesize_nftest).

A simple example of a nf-test directory in nf-core/modules can be found [here](https://github.com/nf-core/modules/tree/master/modules/nf-core/fastqc/tests).

#### Philosophy of nf-tests

- Each module contains a `tests/` folder beside the `main.nf` of the module itself, containing the test files
- Test files come with a [snapshot](https://code.askimed.com/nf-test/docs/assertions/snapshots/) of module output channels

#### nf-test guidelines for a simple un-chained module

- Some modules MAY require additional parameters added to the test command to successfully run. These can be specified with an `ext.args` variable within the process scope of the `nextflow.config` file that exists alongside the test files themselves (and is automatically loaded when the test workflow `main.nf` is executed).

If your module requires a a `nextflow.config` file to run, create the file to the module's `tests/` directory and add the additional parameters there.

```bash
touch modules/nf-core/<tool>/<subtool>/tests/nextflow.config
```

Then add the path to the `main.nf.test` file.

```groovy title="main.nf.test"
process "MODULE"
config "./nextflow.config"
```

- When your test data is too big, the tests take too long or require too much resources, you can opt to run your tests in stub mode by adding the following option:

```groovy title="main.nf.test"
options "-stub"
```

:::note
this can be added at the top of `main.nf.test` to have all tests run in stub mode or this can also be added to a single test
:::

- You can find examples of different nf-tests assertions on [this tutorial](https://nf-co.re/docs/contributing/tutorials/nf-test_assertions).

#### nf-test guidelines for a chained module

- For modules that involve running more than one process to generate required test-data (aka chained modules), nf-test provides a [setup](https://code.askimed.com/nf-test/docs/testcases/setup/) method.

- For example, the module `abricate/summary` requires the process `abricate/run` to be run prior and takes its output as input. The `setup` method is to be declared before the primary `when` block in the test file as shown below:

```groovy title="main.nf.test"
setup {

            run("ABRICATE_RUN") {
                script "../../run/main.nf"
                process {
                    """
                    input[0] =  Channel.fromList([
                        tuple([ id:'test1', single_end:false ], // meta map
                            file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)),
                        tuple([ id:'test2', single_end:false ],
                            file(params.test_data['haemophilus_influenzae']['genome']['genome_fna_gz'], checkIfExists: true))
                    ])
                    """
                }
            }
        }
```

:::note
The setup method can run more than one process each enclosed in their own `run` block
:::

- Then, the output of setup process/es can be provided as input in the `process` section of `when` block

```groovy title="main.nf.test"
input[0] = ABRICATE_RUN.out.report.collect{ meta, report -> report }.map{ report -> [[ id: 'test_summary'], report]}
```

- Next, in the `then` block we can write our assertions that are used to verify the test. A test can have multiple assertions but, we recommend enclosing all assertions in a `assertAll()` block as shown below:

```groovy title="main.nf.test"
assertAll(
            { assert process.success },
            { assert snapshot(process.out).match() }
          )
```

- the `main.nf.test` file for chained modules will finally look as shown below:

```groovy title="main.nf.test"
nextflow_process {

    name "Test Process ABRICATE_SUMMARY"
    script "../main.nf"
    process "ABRICATE_SUMMARY"
    tag "modules"
    tag "modules_nfcore"
    tag "abricate"
    tag "abricate/summary"

    test("bacteroides_fragilis - genome_fna_gz") {

        setup {
            run("ABRICATE_RUN") {
                script "../../run/main.nf"
                process {
                """
                input[0] = Channel.fromList([
                                tuple([ id:'test1', single_end:false ], // meta map
                                    file(params.test_data['bacteroides_fragilis']['genome']['genome_fna_gz'], checkIfExists: true)),
                                tuple([ id:'test2', single_end:false ],
                                    file(params.test_data['haemophilus_influenzae']['genome']['genome_fna_gz'], checkIfExists: true))
                            ])
                """
            }
            }
        }

        when {
            process {
                """
                input[0] = ABRICATE_RUN.out.report.collect{ meta, report -> report }.map{ report -> [[ id: 'test_summary'], report]}
                """
            }
        }

        then {
            assertAll(
                { assert process.success },
                { assert snapshot(process.out).match() }
            )
        }
    }
}
```

### Migrating from pytest to nf-test

#### Steps for creating nf-test for a simple un-chained module

- Git checkout a new branch for your module tests.

```bash
git checkout -b <branch>
```

To create the necessary files for nf-test and ensure a smooth transition, we will use the template provided by nf-core/tools.

Here are the steps to follow:

- Use nf-core/tools to create a new module with the same name as the old one with the option `--migrate-pytest`.
  This command will rename the current module directory to `<module>_old` to avoid conflicts with the new module, create a new module, and copy the `main.nf`, `meta.yml` and `environment.yml` files over to preserve the original module code.

```bash
nf-core modules create <tool>/<subtool> --migrate-pytest
```

- (optional) If your module has a `nextflow.config` file to run (e.g. for `ext.args` specification), the command will also copy it to the module's `tests/` directory and the path will be added to the `main.nf.test` file.

```groovy title="main.nf.test"
process "MODULE"
config "./nextflow.config"
```

- When using the `--migrate-pytest` option you will be asked if you want to delete the old module directory and see the content of the old pytests in the terminal, or to keep the old module directory. For the following steps, use the information from the pytest tests to create the new nf-test tests.

- Provide a test name preferably indicating the test-data and file-format used. Example: `test("homo_sapiens - [bam, bai, bed] - fasta - fai")`

:::note
multiple tests are allowed in a single test file
:::

- If migrating an existing module, get the inputs from current pytest files `tests/modules/nf-core/module/main.nf` and provide as positional inputs `input[0]` in nf-test file

```groovy title="main.nf.test"
input[0] = [
            [id:"ref"],
            file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
           ]
```

- Next, in the `then` block we can write our assertions that are used to verify the test. A test can have multiple assertions but, we recommend enclosing all assertions in a `assertAll()` block as shown below:

```groovy title="main.nf.test"
assertAll(
            { assert process.success },
            { assert snapshot(process.out).match() }
          )
```

- Run the test to create a snapshot of your module test. This will create a `main.nf.test.snap` file

```bash
nf-core modules test <tool>/<subtool>
```

If you chose to not remove the old module directory with nf-core/tools:

- Remove the corresponding tags from `tests/config/pytest_modules.yml` so that py-tests for the module will be skipped on github CI.

- Remove the corresponding pytest files in `tests/modules/nf-core`

```bash
rm -r tests/modules/nf-core/<tool>/<subtool>
```

- Remove the old module

```bash
rm -r modules/nf-core/<tool>/<subtool>_old
```

- Check if everything is according to the nf-core guidelines with:

```bash
nf-core modules lint <tool>/<subtool>
```

- create a PR and add the `nf-test` label to it.

#### Steps for creating nf-test for chained modules

- Follow the steps listed above for simple modules for test generation, tags and test-name

- Refer to the section [nf-test guidelines for a chained module](#nf-test-guidelines-for-a-chained-module)

- Run the test to create a snapshot of your module test. This will create a `.nf.test.snap` file

```bash
nf-core modules test <tool>/<sub-tool>
```

- Add the corresponding module tag from `tests/config/pytest_modules.yml` to the `tags.yml` in `modules/nf-core/<module>/tests/`.

```yaml title="tags.yml"
<tool>/<sub-tool>:
  - modules/nf-core/<tool>/<sub-tool-1>/**
  - modules/nf-core/<tool>/<sub-tool-2>/**
```

:::note
Remove the corresponding tags from `tests/config/pytest_modules.yml` so that py-tests for the module will be skipped on github CI
:::

- create PR and add the `nf-test` label to it.

:::info
The implementation of nf-test in nf-core is still evolving. Things might still change and the information might here might be outdated. Please report any issues you encounter [on the nf-core/website repository](https://github.com/nf-core/website/issues/new?assignees=&labels=bug&projects=&template=bug_report.md) and the `nf-test` channel on nf-core slack.

<!-- NOTE: update when nf-core/tools gets nf-test support -->

:::

### Uploading to `nf-core/modules`

[Fork](https://help.github.com/articles/fork-a-repo/) the `nf-core/modules` repository to your own GitHub account. Within the local clone of your fork add the module file to the `modules/` directory. Please try and keep PRs as atomic as possible to aid the reviewing process - ideally, one module addition/update per PR.

Commit and push these changes to your local clone on GitHub, and then [create a pull request](https://help.github.com/articles/creating-a-pull-request-from-a-fork/) on the `nf-core/modules` GitHub repo with the appropriate information.

When you are happy with your pull request, please <span class="x x-first x-last">select </span>the `Ready for Review` label on the GitHub PR tab, and providing that everything adheres to nf-core guidelines we will endeavour to approve your pull request as soon as possible. We also recommend to request reviews from the `nf-core/modules-team`<span class="x x-first x-last"> so </span>a core team of volunteers <span class="x x-first x-last">can try</span> to <span class="x x-first x-last">review </span>your <span class="x x-first x-last">PR</span> as fast as possible.

Once you<span class="x x-first x-last"> are </span>familiar with the module submission process, please consider joining the<span class="x x-first x-last"> reviewing</span> team by asking on the `#modules` slack channel.

### Talks

:::warning
These may include references to an older syntax, however the general idea remains the same
:::

<div class="ratio ratio-16x9">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/xuNYATGFuw4" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
</div>

<div class="ratio ratio-16x9">
     <iframe src="https://widgets.figshare.com/articles/16825369/embed?show_title=1" width="568" height="351" allowfullscreen frameborder="0"></iframe>
</div>

## New module guidelines and PR review checklist

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

### General

1.  All mandatory and optional input files MUST be included in `input` channel definitions.

2.  Non-file mandatory arguments or arguments needed to modify the command to make the module run with no error, SHOULD be provided as value channels (for example `lib_type` in [salmon/quant](https://github.com/nf-core/modules/blob/master/modules/nf-core/salmon/quant/main.nf)) - see 'Input/output options' below.

3.  All _non-mandatory_ command-line tool _non-file_ arguments MUST be provided as a string via the `$task.ext.args` variable.

    - The value of `task.ext.args` is supplied from the `modules.config` file by assigning a closure that returns a string value to `ext.args`. The closure is necessary to update parameters supplied in a config with `-c`.

      ```groovy {2} title="<module>.nf"
      script:
      def args = task.ext.args ?: ''
      def prefix = task.ext.prefix ?: "${meta.id}"
      """
      fastqc \\
          $args \\
            <...>
      """
      ```

      ```groovy {2-5} title="modules.config"
          withName: <module> {
              ext.args = { [                                                        // Assign a closure which returns a string
                  '--quiet',
                  params.fastqc_kmer_size ? "-k ${params.fastqc_kmer_size}" : ''    // Parameter dependent values can be provided like so
              ].join(' ') }                                                         // Join converts the list here to a string.
              ext.prefix = { "${meta.id}" }                                         // A closure can be used to access variables defined in the script
          }
      }
      ```

    <details markdown="1">
      <summary>Rationale</summary>
      A disadvantage of passing arguments via ext.args is that it splits up how information is passed to a module, which can be difficult to understand where module inputs are defined.

    The justification behind using the `ext.args` is to provide more flexibility to users. As `ext.args` is derived from the configuration (e.g. `modules.config`), advanced users can overwrite the default `ext.args` and supply their own arguments to modify the behaviour of a module. This can increase the capabilities of a pipeline beyond what the original developers intended.

    Initially these were passed via the main workflow script using custom functions (e.g. `addParams`) and other additional nf-core custom methods, but this had a syntax overhead and other limitations that were found to be more difficult to use and understand by pipeline developers. Therefore using the 'native' `ext` functionality provided by Nextflow was easier to understand, maintain and use.

    Note that sample-specific parameters can still be provided to an instance of a process by storing these in `meta`, and providing these to the `ext.args` definition in `modules.config`. A closure is used to make Nextflow evaluate the code in the code in the string.

    ```nextflow
    ext.args = { "--id ${meta.id}" }
    ```

      </details>

4.  Software that can be piped together SHOULD be added to separate module files
    unless there is a run-time, storage advantage in implementing in this way. For example,
    using a combination of `bwa` and `samtools` to output a BAM file instead of a SAM file:

    ```bash
    bwa mem $args | samtools view $args2 -B -T ref.fasta
    ```

    :::info
    The addition of multi-tool modules to nf-core/modules adds increased burden on the nf-core
    maintainers. Where possible, if a multi-tool module is desired, it should be implemented as
    a local module in the nf-core pipeline. If another nf-core pipeline also desires to
    use this module, a PR can be made adding it to nf-core/modules.

    For guidelines regarding multi-tool modules, please search this page for the phrase `multi-tool`.

    Existing local multi-tool modules can be searched for using the Github search box, searching across
    the nf-core org for terms such as `args2` `samtools` `collate` `fastq`.

    ```
    org:nf-core args2 samtools collate fastq
    ```

    Modules intended to batch process files by parallelizing repeated calls to a tool, for example with
    `xargs` or `parallel`, also fall under the category of multi-tool modules. Multi-tool modules
    should chain tools in an explicit order given by the module name, e.g. `SAMTOOLS/COLLATEFASTQ`.
    :::

5.  Each tool in a multi-tool module MUST have an `$args` e.g.,

    ```bash
    bwa mem $args | samtools view $args2 -B -T ref.fasta | samtools sort $args3
    ```

    or

    ```bash
    <tool> \\
       <subcommand> \\
       $args

    gzip \\
        $args2
    ```

    The numbering of each `$args` variable MUST correspond to the order of the tools, and MUST be documented in `meta.yml`. E.g. in the first example, `bwa mem` is the first tool so is given `$args`, `samtools view` is the second tool so is `$args2`, etc.

6.  Modules MUST NOT use 'custom' hardcoded `meta` fields. The only accepted 'standard' meta fields are `meta.id` or `meta.single_end`. Proposals for other 'standard' fields for other disciplines must be discussed with the maintainers team.

    <details markdown="1">
      <summary>Rationale</summary>
      Modules should be written to allow as much flexibility to pipeline developers as possible.

    Hardcoding `meta` fields in a module will reduce the freedom of developers to use their own names for metadata, which would make more sense in that particular context.

    As all non-mandatory arguments must go via `$args`, pipeline developers can insert such `meta` information into `$args` with whatever name they wish.

    So, in the module code we DO NOT do:

    ```bash
    my_command -r ${meta.strand} input.txt output.txt
    ```

    ... but rather, in `modules.conf`

    ```nextflow
    ext.args = { "--r ${meta.<pipeline_authors_choice_of_name>}" }
    ```

    ... and then in the module code `main.nf`:

    ```bash
    my_command $args input.txt output.txt
    ```

      </details>

7.  Where applicable, the usage and generation of compressed files SHOULD be enforced as input and output, respectively:

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

8.  Where applicable, each module command MUST emit a file `versions.yml` containing the version number for each tool executed by the module, e.g.

    ```bash
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
        samtools: \$( samtools --version |& sed '1!d ; s/samtools //' )
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

    <details class="mb-3">
    <summary>Tips for extracting the version string</summary>

    `sed{:bash}` is a powerful stream editor that can be used to manipulate the input text into the desired output.
    Start by piping the output of the version command to `sed{:bash}` and try to select the line with the version number:

    ```bash
    tool --version | sed '1!d'
    ```

    - `sed '1!d'{:bash}` Extracts only line 1 of the output printed by `tools --version{:bash}`.
    - The line to process can also be selected using a pattern instead of a number: `sed '/pattern/!d'{:bash}`, e.g. `sed '/version:/!d'{:bash}`.
    - If the line extraction hasn't worked, then it's likely the version information is written to stderr, rather than stdout.
      In this case capture stderr using `|&{:bash}` which is shorthand for `2>&1 |{:bash}`.
    - `sed 's/pattern/replacement/'{:bash}` can be used to remove parts of a string. `.` matches any character, `+` matches 1 or more times.
    - You can separate `sed{:bash}` commands using `;`. Often the pattern : `sed 'filter line ; replace string'{:bash}` is enough to get the version number.
    - It is not necessary to use `echo`, `head`, `tail`, or `grep`.
    - Use `|| true` for tools that exit with a non-zero error code: `command --version || true{:bash}` or `command --version | sed ... || true{:bash}`.
    </details>

    We chose a [HEREDOC](https://tldp.org/LDP/abs/html/here-docs.html) over piping into the versions file
    line-by-line as we believe the latter makes it easy to accidentally overwrite the file. Moreover, the exit status
    of the sub-shells evaluated in within the HEREDOC is ignored, ensuring that a tool's version command does
    not erroneously terminate the module.

    If the software is unable to output a version number on the command-line then a variable called `VERSION` can be manually
    specified to provide this information e.g. [homer/annotatepeaks module](https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf).
    Please include the accompanying comments above the software packing directives and beside the version string.

    ```nextflow {4,15,21}
    process TOOL {
        ...

        // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
        conda (params.enable_conda ? "bioconda::tool=0.9.1:" : null)
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/tool:0.9.1--pl526hc9558a2_3' :
            'biocontainers/tool:0.9.1--pl526hc9558a2_3' }"

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

9.  The process definition MUST NOT change the `when` statement. `when` conditions can instead be supplied using the `process.ext.when` directive in
    a configuration file.

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

10. In some cases, STDOUT and STDERR need to be saved to file, for example for reporting purposes. Use the shell command `tee` to redirect the
    streams to both file and it's original stream. This allows for the streams to be captured by the job scheduler's stream logging capabilities
    and print them to screen when Nextflow encounters an error. In particular, when using `process.scratch`, the log files may not be preserved when
    the job scheduler relinquishes the job allocation.

    ```nextflow {7-8}
    script:
    """
    tool \\
      --input $input \\
      --threads $task.cpus \\
      --output_prefix $prefix \\
      2> >( tee ${prefix}.stderr.log >&2 ) \\
      | tee ${prefix}.stdout.log
    """
    ```

11. Occasionally, some tools do not exit with the expected exit code 0 upon successful use of the tool. In these cases one can use
    the `||` operator to run another useful command when the exit code is not 0 (for example, testing if a file is not size 0).

    ```nextflow {6}
    script:
    """
    tool \\
      --input $input \\
      --summary ${prefix}.summary \\
      || test -s ${prefix}.summary
    """
    ```

    See the [Bash manual on file operators](https://tldp.org/LDP/abs/html/fto.html) for examples of properties of files which could be tested.
    Alternate suggestions include using `grep -c` to search for a valid string match, or other tool which will appropriately
    error when the expected output is not successfully created.

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

1. Input channel `path` declarations MUST be defined for all _possible_ input files (i.e. both required and optional files).

   - Directly associated auxiliary files to an input file MAY be defined within the same input channel alongside the main input channel (e.g. [BAM and BAI](https://github.com/nf-core/modules/blob/e937c7950af70930d1f34bb961403d9d2aa81c7d/modules/samtools/flagstat/main.nf#L22)).
   - Other generic auxiliary files used across different input files (e.g. common reference sequences) MAY be defined using a dedicated input channel (e.g. [reference files](https://github.com/nf-core/modules/blob/3cabc95d0ed8a5a4e07b8f9b1d1f7ff9a70f61e1/modules/bwa/mem/main.nf#L21-L23)).

2. Input channel `val` declarations SHOULD be defined for all mandatory non-file inputs that are essential for the functioning of the tool (e.g. parameters, flags etc).

   - Mandatory non-file inputs are options that the tool MUST have to be able to be run.
   - These non-file inputs are typically booleans or strings, and must be documented as such in the corresponding entry in the `meta.yaml`.
   - Options, flags, parameters that are _not_ required by the tool to function should NOT be included - rather these can be passed via `ext.args`.

       <details markdown="1">
       <summary>Rationale</summary>
       It was decided by a [vote](https://nfcore.slack.com/archives/C043UU89KKQ/p1677581560661679) amongst interested parties within the 2023 Maintainers group on 2023-02-28 to allow non-file mandatory input channels.

     The reasoning behind this was that it is important to have documented (using the existing display on the website) the bare minimum information required for a module to run. It also allows module code to consume parameter values without parsing them out of the `ext.args` string and reduces possible risks of entire breakage of modules with future [expected config changes](https://github.com/nextflow-io/nextflow/issues/2723) at a Nextflow level.

     Downsides to this approach are readability (now multiple places must be checked on how to modify a module execution - modules.conf `ext.args`, the module invocation in pipeline code etc.), and reduced user freedom. However it was felt that it was more important for stability in and 'installation' and 'execution' of modules was preferred (e.g. for tools that require position arguments etc.)

       </details>

       <details markdown="1">
       <summary>Inputs particular cases</summary>
        When one and only one of multiple argument are required:

     - If they all are string argument : use 1 argument that will be equal to the string

     e.g. Parameter model of [glimpse2 chunk](https://nf-co.re/modules/glimpse2_chunk)

     - If some are files put them all in one channel and test if only one is present

     e.g. Grouping output parameters of [glimpse2 concordance](https://nf-co.re/modules/glimpse2_concordance)

     `if (((file1 ? 1:0) + (val1 ? 1:0) + (val2 ? 1:0)) != 1) error "One and only one argument required"`
       </details>

3. Named file extensions MUST be emitted for ALL output channels e.g. `path "*.txt", emit: txt`.

4. Optional inputs are not currently supported by Nextflow. However, passing an empty list (`[]`) instead of a file as a module parameter can be used to work around this issue.

For example, having a module (`MY_MODULE`) that can take a `cram` channel and an optional `fasta` channel as input, can be used in the following ways:

```nextflow
MY_MODULE(cram, [])     // fasta is optional, the module will run without the fasta present
MY_MODULE(cram, fasta)  // execution of the module will need an element in the fasta channel
```

5. Optional outputs SHOULD be marked as optional:

   ```nextflow
   tuple val(meta), path('*.tab'), emit: tab,  optional: true
   ```

6. Each output file SHOULD be emitted in its own channel (and no more than one), along with the `meta` map if provided ( the exception is the versions.yml ).

### Documentation

1. Each module MUST have a `meta.yaml` in the same directory as the `main.nf` of the module itself.

2. Keywords SHOULD be sufficient to make the module findable through research domain, data types, and tool function keywords

   - Keywords MUST NOT just be solely of the (sub)tool name

   :::info
   For multi-tool modules, please add the keyword `multi-tool`, as well as all the (sub)tools involved.
   :::

3. Keywords MUST be all lower case

4. The tools section MUST list every tool used in the module. For example

   ```yml
   tools:
     - bowtie2: <....>
     - samtools: <....>
   ```

5. The tools section MUST have a `args_id:` field for every tool in the module that describes which `$args` (`$args2`, `$args3`) variable is used for that specific module. A single tool module will only have `args_id: "$args"`.

   ```yml
   tools:
     - bowtie2:
         <...>
         args_id: "$args"
     - samtools:
         <...>
         args_id: "$args2"
   ```

6. Input and Output sections of the `meta.yaml` SHOULD only have entries of input and output channels

7. Input and output tuples MUST be split into separate entries

   - i.e., `meta` should be a separate entry to the `file` it is associated with

8. Input/output types MUST only be of the following categories: `map`, `file`, `directory`, `string`, `boolean`, `integer`, `float`, `boolean`, `list`

9. Input/output entries MUST match a corresponding channel in the module itself

   - There should be a one-to-one relationship between the module and the `meta.yaml`

   - Input/output entries MUST NOT combine multiple output channels

10. Input/output descriptions SHOULD be descriptive of the contents of file

- i.e., not just 'A TSV file'

11. Input/output patterns (if present) MUST follow a [Java glob pattern](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob)

12. Input entries should be marked as Mandatory or Optional

### Module parameters

1. A module file SHOULD only define input and output files as command-line parameters to be executed within the process.

2. All `params` within the module MUST be initialised and used in the local context of the module. In other words, named `params` defined in the parent workflow MUST NOT be assumed to be passed to the module to allow developers to call their parameters whatever they want. In general, it may be more suitable to use additional `input` value channels to cater for such scenarios.

3. If the tool supports multi-threading then you MUST provide the appropriate parameter using the Nextflow `task` variable e.g. `--threads $task.cpus`.

4. Any parameters that need to be evaluated in the context of a particular sample e.g. single-end/paired-end data MUST also be defined within the process.

### Resource requirements

1. An appropriate resource `label` MUST be provided for the module as listed in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/conf/base.config#L29-L46) e.g. `process_single`, `process_low`, `process_medium` or `process_high`.

2. If the tool supports multi-threading then you MUST provide the appropriate parameter using the Nextflow `task` variable e.g. `--threads $task.cpus`. If the tool does not support multi-threading, consider `process_single` unless large amounts of RAM are required.

3. If a module contains _multiple_ tools that supports multi-threading (e.g. [piping output into a samtools command](https://github.com/nf-core/modules/blob/c4cc1db284faba9fc4896f64bddf7703cedc7430/modules/nf-core/bowtie2/align/main.nf#L47-L54)), you can assign CPUs per tool.
   - Note that [`task.cpus`] is supplied unchanged when a process uses multiple cores
   - If one tool is multi-threaded and another uses a single thread, you can specify directly in the command itself e.g. with [`${task.cpus}`](https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/sampe/main.nf#L34)

### Software requirements

[BioContainers](https://biocontainers.pro/#/) is a registry of Docker and Singularity containers automatically created from all of the software packages on [Bioconda](https://bioconda.github.io/). Where possible we will use BioContainers to fetch pre-built software containers and Bioconda to install software using Conda.

1. Software requirements SHOULD be declared within the module file using the Nextflow `container` directive. For single-tool BioContainers, the `nf-core modules create` command will automatically fetch and fill-in the appropriate Conda / Docker / Singularity definitions by parsing the information provided in the first part of the module name:

```nextflow
conda "bioconda::fastqc=0.11.9"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
    'biocontainers/fastqc:0.11.9--0' }"
```

2. If the software is available on Conda it MUST also be defined using the Nextflow `conda` directive. Using `bioconda::bwa=0.7.17` as an example, software MUST be pinned to the channel (i.e. `bioconda`) and version (i.e. `0.7.17`). Conda packages MUST not be pinned to a build because they can vary on different platforms.

3. If required, multi-tool containers may also be available on BioContainers e.g. [`bwa` and `samtools`](https://biocontainers.pro/#/tools/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40). You can install and use the [`galaxy-tool-util`](https://anaconda.org/bioconda/galaxy-tool-util) package to search for both single- and multi-tool containers available in Conda, Docker and Singularity format. e.g. to search for Docker (hosted on Quay.io) and Singularity multi-tool containers with both `bowtie` and `samtools` installed you can use the following command:

```console
mulled-search --destination quay singularity --channel bioconda --search bowtie samtools | grep "mulled"
```

:::note
Build information for all tools within a multi-tool container can be obtained in the `/usr/local/conda-meta/history` file within the container.
:::

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

   - If the multi-tool container already exists and you want to obtain the `mulled-*` path, you can use (this)[https://midnighter.github.io/mulled] helper tool.

5. If the software is not available on Bioconda a `Dockerfile` MUST be provided within the module directory. We will use GitHub Actions to auto-build the containers on the [GitHub Packages registry](https://github.com/features/packages).

### Misc

1. All code must be aligned to follow the '[Harshil Alignment™️](#what-is-the-harshil-alignment)' format.

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

### PR Review Checklist

A PR review is the process of examining a new modules' submission or the changes proposed to a module. The reviewer provides constructive feedback on those changes before they are merged into the nf-core repository. The goal of a PR review is to ensure that the code meets the coding standards of the project, is consistent and of high-quality. While the team of [maintainers](https://github.com/orgs/nf-core/teams/maintainers/members) is responsible for overseeing the PR review process for modules, these guidelines can assist community members in reviewing PRs and ensure that the review process is consistent and effective. The following is a collection of community suggestions to have into account during the review process.

#### General reviews of submissions to modules:

- Ensure all checks pass, including linting, conda, singularity, and docker.
- Check that the module is suitable for offline running, without automatic database downloads assumed.
- If running docker containers, check that Nextflow changes the `--entrypoint` to `/bin/bash` and that environment variables used by certain programs (e.g., Busco, Merqury) are sourced again to use them in container settings.
- Check that it adheres to nf-core coding standards (e.g. use of meta map).
- Check that the code is readable and the formatting is correct (e.g. indenting, extra spaces).

#### In `modules/nf-core/modulename/main.nf`:

- Check that all optional parameters are in the `$args` section.
- Check that the software version extraction command is optimized, if required.
- Check if the bioconda version of the tool is the latest version.
- Ensure that temporary unzipped files are removed to avoid mitigating benefits and worsening problems.
- Ensure that large outputs are compressed with the correct tool (follow guidelines for gzip vs bzip2 vs other options).

#### In `../tests/modules/nf-core/modulename/main.nf` and `../tests/modules/nf-core/modulename/meta.yml`:

- Check that there are tests for all outputs, including optional ones.
- Check that the `meta.yml` file has correct documentation links and patterns of files.
- Run the tool help and check that important input (usually optional) has not been missed.
- Check that all outputs are captured by running pytest (e.g. on Gitpod).

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

### Generating a `meta map` from file pairs

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

### Sorting samples by groups

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

### Combining channel on meta subset

Sometimes it is necessary to combine multiple channels based on a subset of the meta maps.
Unfortunately this is not yet supported as the argument `by` isn't a closure in `.combine()` and `.join()` and it probably won't ([Nextflow issue #3175](https://github.com/nextflow-io/nextflow/issues/3175)).

To bypass this restriction one of the solution is to create a new map with only the necessary keys and make the junction on it. Here is an example:

```nextflow
ch_input = [[["id":"Ind1","ref":"RefA"],"file1"],[["id":"Ind2","ref":"RefB"],"file2"]]
ch_ref   = [[["ref":"RefA"],"fileA"],[["ref":"RefB"],"fileB"]]

ch_join  = ch_input
            .map{metaIR, file -> [metaIR.subMap(["ref"]), metaIR, file]}
            .combine(chr_ref)
            .map{metaR, metaIR, file, ref -> [metaIR, file, ref]}
```

### Modify the meta map

There is multiple ways to modify the meta map.
Here are some examples:

```nextflow
// Add to map - adding two maps makes a new Map object
ch.map { meta, files -> [ meta + [ single_end: files instanceof Path ], files ] }

// Remove certain keys (and their entries) from a map
ch.map { meta, files -> [ meta.subMap( ['id','rg'] ), files ] }
  // OR by specifying what not to include
ch.map { meta, files -> [ meta.findAll { ! it.key in ['single_end'] }, files ] }

// Split a map - use both methods of removing keys ( there is a split method for Maps, but the results are not Maps )
ch.map { meta, files -> def keyset = ['id', 'read_group']; [ meta.subMap(keyset), meta.findAll { ! it.key in keyset },  files ] }
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

:::note
The order of the numeric ID of `args` must match the order of the tools as used in the module.
:::

To see some more advanced examples of these keys in use see:

- [Set ext.args based on parameter settings](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/modules.config#L222-L226)
- [Set ext.prefix based on task inputs](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/modules.config#L297)
- [Set ext.args based on both parameters and task inputs](https://github.com/nf-core/rnaseq/blob/e049f51f0214b2aef7624b9dd496a404a7c34d14/conf/modules.config#L377-L381)

## Advanced pattern

### Multimaping

It is possible with `multiMap` to split a channel in to and to call them separately afterwards.

```nextflow
ch_input = reads.combine(db).multiMap{ it ->
   reads: it[0]
   db: it[1]
}
MODULE(ch_input.reads, ch_input.db)
```

### Adding additional information to the meta map

It is possible to combine a input channel with a set of parameters as follows:

```nextflow
ch_input.flatMap { meta, filetype ->
    [300, 500, 1000].collect {
      def new_meta = meta.clone()
      new_meta.window_size = it
      [ new_meta, filetype]
    }
}
```

You can also combine this technique with others for more processing:

```nextflow
workflow {

    input = [
        [
            [ patient: 'sample', sample: 'test', id: 'test' ],
            file ("chr21_23355001-46709983.bed")
        ],
        [
            [ patient: 'sample', sample: 'test', id: 'test' ],
            file ("chr21_2-23354000.bed")
        ],
        [
            [ patient: 'sample2', sample: 'test5', id: 'test' ],
            file ("chr21_23355001-46709983.bed")
        ],
        [
            [ patient: 'sample2', sample: 'test5', id: 'test' ],
            file ("chr21_2-23354000.bed")
        ]
    ]
    Channel.fromList ( input )
        .map { meta, intervals ->
            new_meta = meta.clone()
            new_meta.id = intervals.baseName != "no_intervals" ? new_meta.sample + "_" + intervals.baseName : new_meta.sample
            intervals = intervals.baseName != "no_intervals" ? intervals : []
            [new_meta, intervals]
        }.view { meta, intervals -> meta.id }
}
```

## What is the Harshil Alignment

The Harshil Alignment™️ format is the whitespace-happy code style that was introduced by a certain core member to get on everyone's nerves, but then make subsequently develop Stockholm Syndrome so that no-one in nf-core else now can look at Nextflow code without it.

The Harshil Alignment™️ format involves ensuring that common punctuation across multiple lines in a group are placed in the same location as each other.

There are many places where the format can be applied - it's not just code, it can also applies to comment formatting - however common examples are as follows:

### Curly Bracket Example

❌ Bad

```nextflow
include { SAMTOOLS_SORT } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'
```

✅ Good

```nextflow
include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'
```

### Equals Example

❌ Bad

```nextflow
stats = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
```

✅ Good

```nextflow
stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats    ] ]
flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
```

### Comma Example

❌ Bad

```nextflow
tuple val(meta), path("*.bam"), emit: bam, optional:true
tuple val(meta), path("*.log"), emit: log
tuple val(meta), path("*fastq.gz"), emit: fastq, optional:true
path  "versions.yml", emit: versions
```

✅ Good

```nextflow
tuple val(meta), path("*.bam")    , emit: bam     , optional:true
tuple val(meta), path("*.log")    , emit: log
tuple val(meta), path("*fastq.gz"), emit: fastq   , optional:true
path  "versions.yml"              , emit: versions
```

### Colon Example (Comments)

```nextflow
take:
print_version        // boolean: print version
dump_parameters      // boolean: dump parameters
outdir               // path: base directory used to publish pipeline results
check_conda_channels // boolean: check conda channels
```

✅ Good

```nextflow
take:
print_version        // boolean: print version
dump_parameters      // boolean: dump parameters
outdir               //    path: base directory used to publish pipeline results
check_conda_channels // boolean: check conda channels
```

## Help

For further information or help, don't hesitate to get in touch on [Slack `#modules` channel](https://nfcore.slack.com/channels/modules) (you can join with [this invite](https://nf-co.re/join/slack)).
