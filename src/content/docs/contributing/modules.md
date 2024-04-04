---
title: Guidelines for Modules
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

1. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Conda`](https://conda.io/miniconda.html)

:::tip{title="Single step conda installation" collapse}

If you use the conda package manager you can setup a new environment and install all dependencies for the new module workflow in one step with:

```bash
  conda create -n nf-core -c bioconda "nextflow>=21.04.0" "nf-core>=2.7" nf-test
  conda activate nf-core
```

and proceed with Step 5.
:::

2. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)
3. Install the latest version of [`nf-core/tools`](https://github.com/nf-core/tools#installation) (`>=2.7`)
4. Install [`nf-test`](https://code.askimed.com/nf-test/installation/)
5. [Fork and clone the nf-core/modules repo locally](#uploading-to-nf-coremodules)
6. Setup up [pre-commit](https://pre-commit.com/) (comes packaged with [`nf-core/tools`](https://github.com/nf-core/tools#installation), watch the [pre-commit bytesize talk](https://www.youtube.com/watch?v=08d6zv6zvdM&t=215) if you want to know more about it) to ensure that your code is linted and formatted correctly before you commit it to the repository

   ```bash
   pre-commit install
   ```

7. Set up git on your computer by adding a new git remote of the main nf-core git repo called `upstream`

   ```bash
   git remote add upstream https://github.com/nf-core/modules.git
   ```

   Make a new branch for your module and check it out

   ```bash
   git checkout -b fastqc
   ```

8. Create a module using the [nf-core DSL2 module template](https://github.com/nf-core/tools/blob/master/nf_core/module-template/main.nf):

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

   1. [`./modules/nf-core/fastqc/main.nf`](https://github.com/nf-core/modules/blob/master/modules/nf-core/fastqc/main.nf)

      This is the main script containing the `process` definition for the module. You will see an extensive number of `TODO` statements to help guide you to fill in the appropriate sections and to ensure that you adhere to the guidelines we have set for module submissions.

   2. [`./modules/nf-core/fastqc/meta.yml`](https://github.com/nf-core/modules/blob/master/modules/nf-core/fastqc/meta.yml)

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

- Please adhere to the [test-data specifications](/docs/contributing/component-specifications/test_data.md) when adding new test-data

- In order to keep the size of the test data repository as minimal as possible, pre-existing files from [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) MUST be reused if at all possible.

- Test files MUST be kept as tiny as possible.

- If the appropriate test data doesn't exist in the `modules` branch of [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) please contact us on the [nf-core Slack `#modules` channel](https://nfcore.slack.com/channels/modules) (you can join with [this invite](https://nf-co.re/join/slack)) to discuss possible options.

- It may not be possible to add test data for some modules e.g. if the input data is too large or requires a local database. In these scenarios, it is recommended to use the Nextflow [`stub`](https://www.nextflow.io/docs/latest/process.html#stub) feature to test the module. Please refer to the [`gtdbtk/classify`](https://github.com/nf-core/modules/blob/79d38a306bdaf07000e0d6f300684d3ed38c8919/modules/gtdbtk/classifywf/main.nf#L66) module and its corresponding [test script](https://github.com/nf-core/modules/blob/79d38a306bdaf07000e0d6f300684d3ed38c8919/tests/modules/gtdbtk/classifywf/main.nf#L20) to understand how to use this feature for your module development.

If a new test dataset is added to [`tests/config/test_data.config`](https://github.com/nf-core/modules/blob/master/tests/config/test_data.config), check that the config name of the added file(s) follows the scheme of the entire file name with dots replaced with underscores.

For example: the nf-core/test-datasets file `genomics/sarscov2/genome/genome.fasta` labelled as `genome_fasta`, or `genomics/sarscov2/genome/genome.fasta.fai` as `genome_fasta_fai`.

### Publishing results

Results are published using Nextflow's native [`publishDir`](https://www.nextflow.io/docs/latest/process.html#publishdir) directive defined in the `modules.config` of a workflow (see [here](https://github.com/nf-core/rnaseq/blob/f7702d5b76a1351e2e7796a5ed3f59943a139fbf/conf/modules.config#L100-L106) for an example.) Results were earlier published using a custom `publishDir` definition, using a Groovy Map defined by `params.modules`.

### Using a stub test when required test data is too big

If the module absolute cannot run using tiny test data, there is a possibility to add [stub-run](https://www.nextflow.io/docs/edge/process.html#stub) to the test.yml. In this case it is required to test the module using larger scale data and document how this is done. In addition, an extra script-block labeled `stub:` must be added, and this block must create dummy versions of all expected output files as well as the `versions.yml`. An example is found in the [ascat module](https://github.com/nf-core/modules/blob/master/modules/nf-core/ascat/main.nf). In the `test.yml` the `-stub-run` argument is written as well as the md5sums for each of the files that are added in the stub-block. This causes the stub-code block to be activated when the unit test is run ([example](https://github.com/nf-core/modules/blob/master/tests/modules/nf-core/ascat/test.yml)):

```console
nextflow run tests/modules/<nameofmodule> -entry test_<nameofmodule> -c tests/config/nextflow.config -stub-run
```

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

### PR Review Checklist

A PR review is the process of examining a new modules' submission or the changes proposed to a module. The reviewer provides constructive feedback on those changes before they are merged into the nf-core repository.The goal of a PR review is to ensure that the code meets the coding standards of the project, is consistent and of high-quality.

While the team of [maintainers](https://github.com/orgs/nf-core/teams/maintainers/members) is responsible for overseeing the PR review process for modules, these guidelines can assist community members in reviewing PRs and ensure that the review process is consistent and effective. The following is a collection of community suggestions to have into account during the review process.

#### General reviews of submissions to modules:

In general, the main purpose of the review is to ensure

- All modules adhere to the nf-core [module specifications](/docs/contributing/component-specifications/modules.md)
- Ensure all checks pass, including linting, conda, singularity, and docker.

Otherwise, you can cover most of the specifications by checking for the following:

- The module is suitable for offline running, without automatic database downloads assumed.
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
- Check that all outputs are captured by running nf-test or pytest (e.g. on Gitpod).

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

## Deprecating a module

Sometimes modules or subworkflows become outdated and need to be deprecated (available, but no longer recommended).
These modules or subworkflows should not be deleted as they could be used in private repositories, or used on other
platforms. The recommended procedure is, once the alternative is available on nf-core modules, add a message to the
top of the module code saying this module is deprecated, and an `assert` in the code body to print a deprecation
message like so:

```groovy title="main.nf"
def deprecation_message = """
WARNING: This module has been deprecated. Please use nf-core/modules/path/to/new/module

Reason:
This module is no longer fit for purpose because ...

"""

process OLD_MODULE {
  ...

  script:
  assert false: deprecation_message
}
```

The purpose of the `assert` is to introduce a mechanism which stops the pipeline and alerts the developer when
an automatic update of the module/subworkflow is performed.

## Help

For further information or help, don't hesitate to get in touch on [Slack `#modules` channel](https://nfcore.slack.com/channels/modules) (you can join with [this invite](https://nf-co.re/join/slack)).
