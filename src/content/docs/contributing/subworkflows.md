---
title: Guidelines for Subworkflows
subtitle: Guidelines and reference for DSL2 subworkflows
---

If you decide to upload a subworkflow to `nf-core/modules` then this will ensure that it will become available to all nf-core pipelines, and to everyone within the Nextflow community! See [`subworkflows/nf-core/`](https://github.com/nf-core/modules/tree/master/subworkflows/nf-core) for examples.

<!-- See the [dsl2 modules tutorial](tutorials/dsl2_modules_tutorial) for a step by step guide for how to add a module! -->

## Terminology

The features offered by Nextflow DSL2 can be used in various ways depending on the granularity with which you would like to write pipelines. Please see the listing below for the hierarchy and associated terminology we have decided to use when referring to DSL2 components.

### Module

A `process` that can be used within different pipelines and is as atomic as possible i.e. cannot be split into another module. An example of this would be a module file containing the process definition for a single tool such as `FastQC`. Atomic nf-core module files are available in the [`modules/`](https://github.com/nf-core/modules/tree/master/modules) directory of nf-core/modules along with the required documentation and tests.

### Subworkflow

A chain of multiple modules that offer a higher-level of functionality within the context of a pipeline. For example, a subworkflow to run multiple QC tools with FastQ files as input. Subworkflows should be shipped with the pipeline implementation and if required they should be shared amongst different pipelines directly from there. Shareable nf-core subworkflow files are available in the [`subworkflow/`](https://github.com/nf-core/modules/tree/master/subworkflows) directory of nf-core/modules along with the required documentation and tests.

### Workflow

What DSL1 users would consider an end-to-end pipeline. For example, from one or more inputs to a series of outputs. This can either be implemented using a large monolithic script as with DSL1, or by using a combination of DSL2 modules and subworkflows. nf-core pipelines can have multiple workflows, such as processing different data types for the same ultimate purpose (such as in [nf-core/viralrecon](https://github.com/nf-core/viralrecon/tree/master/workflows))

<!-- ## Writing a new subworkflow reference

See the [dsl2 subworkflows tutorial](tutorials/dsl2_subworkflows_tutorial) for a step by step guide for how to add a module! -->

### Before you start

Please check that the subworkflow you wish to add isn't already in the [`nf-core/modules` repository](https://github.com/nf-core/modules/tree/master/subworkflows/nf-core):

- Use the [`nf-core subworkflows list`](https://github.com/nf-core/tools#list-subworkflows) command
- Check [open pull requests](https://github.com/nf-core/modules/pulls)
- Search [open issues](https://github.com/nf-core/modules/issues)

If the subworkflow doesn't exist on `nf-core/modules`:

- Please create a [new issue](https://github.com/nf-core/modules/issues/new?assignees=&labels=new%20module&template=new_nodule.md&title=new%20module:) before adding it
- Set an appropriate subject for the issue e.g. `new subworkflow: bam_sort_stats_samtools`
- Add yourself to the `Assignees` so we can track who is working on the subworkflow

### Adding a new subworkflow

We have implemented a number of commands in the `nf-core/tools` package to make it incredibly easy for you to create and contribute your own subworkflow to nf-core/modules.

1. Install the latest version of [`nf-core/tools`](https://github.com/nf-core/tools#installation) (`>=2.7`)
2. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)
3. Install [`nf-test`](https://code.askimed.com/nf-test/installation/)
4. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Conda`](https://conda.io/miniconda.html)
5. [Fork and clone the nf-core/modules repo locally](#uploading-to-nf-coremodules)
6. Set up git on your computer by adding a new git remote of the main nf-core git repo called `upstream`

   ```bash
   git remote add upstream https://github.com/nf-core/modules.git
   ```

   Make a new branch for your subworkflow and check it out

   ```bash
   git checkout -b bam_sort_stats_samtools
   ```

7. Create a subworkflow using the [nf-core DSL2 subworkflow template](https://github.com/nf-core/tools/blob/master/nf_core/subworkflow-template/subworkflows/main.nf) in the root of the clone of the nf-core/modules repository:

   ```console
   $ nf-core subworkflows create bam_sort_stats_samtools --author @joebloggs

                                             ,--./,-.
             ___     __   __   __   ___     /,-._.--~\
       |\ | |__  __ /  ` /  \ |__) |__         }  {
       | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                             `._,._,'

       nf-core/tools version 2.11 - https://nf-co.re


   INFO     Repository type: modules
   INFO     Press enter to use default values (shown in brackets) or type your own responses. ctrl+click underlined text to open links.
   INFO     Created following files:
              subworkflows/nf-core/bam_sort_stats_samtools/main.nf
              subworkflows/nf-core/bam_sort_stats_samtools/meta.yml
              subworkflows/nf-core/bam_sort_stats_samtools/tests/tags.yml
              subworkflows/nf-core/bam_sort_stats_samtools/tests/main.nf.test
   ```

All of the files required to add the subworkflow to `nf-core/modules` will be created/edited in the appropriate places. There are at most 5 files to modify:

1. [`./subworkflows/nf-core/bam_sort_stats_samtools/main.nf`](https://github.com/nf-core/modules/blob/master/subworkflows/nf-core/bam_sort_stats_samtools/main.nf)

   This is the main script containing the `workflow` definition for the subworkflow. You will see an extensive number of `TODO` statements to help guide you to fill in the appropriate sections and to ensure that you adhere to the guidelines we have set for module submissions.

2. [`./subworkflows/nf-core/bam_sort_stats_samtools/meta.yml`](https://github.com/nf-core/modules/blob/master/subworkflows/nf-core/bam_sort_stats_samtools/meta.yml)

   This file will be used to store general information about the subworkflow and author details. You will need to add a brief description of the files defined in the `input` and `output` section of the main script since these will be unique to each subworkflow.

3. [`./subworkflows/nf-core/bam_sort_stats_samtools/tests/main.nf.test`](https://github.com/nf-core/modules/blob/master/subworkflows/nf-core/bam_sort_stats_samtools/tests/main.nf.test)

   Every subworkflow MUST have a test workflow. This file will define one or more Nextflow `workflow` definitions that will be used to unit test the output files created by the subworkflow. By default, one `workflow` definition will be added but please feel free to add as many as possible so we can ensure that the subworkflow works on different data types / parameters e.g. separate `workflow` for single-end and paired-end data.

   Minimal test data required for your subworkflow may already exist within the [nf-core/modules repository](https://github.com/nf-core/modules/blob/master/tests/config/test_data.config), in which case you may just have to change a couple of paths in this file - see the [Test data](#test-data) section for more info and guidelines for adding new standardised data if required.

   Refer to the section [writing nf-test tests](#writing-nf-test-tests) for more information on how to write nf-tests

4. Create a snapshot file for subworkflow testing

   ```console
   $ nf-core subworkflows test bam_sort_stats_samtools

                                             ,--./,-.
             ___     __   __   __   ___     /,-._.--~\
       |\ | |__  __ /  ` /  \ |__) |__         }  {
       | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                             `._,._,'

       nf-core/tools version 2.8 - https://nf-co.re


   INFO     Generating nf-test snapshot
   ╭────────────────────────────────────────────────────────── nf-test output ──────────────────────────────────────────────────────────╮
    │                                                                                                                                    │
    │ 🚀 nf-test 0.8.2                                                                                                                   │
    │ https://code.askimed.com/nf-test                                                                                                   │
    │ (c) 2021 - 2023 Lukas Forer and Sebastian Schoenherr                                                                               │
    │                                                                                                                                    │
    │ Found 1 files in test directory.                                                                                                   │
    │                                                                                                                                    │
    │ Test Workflow BAM_SORT_STATS_SAMTOOLS                                                                                              │
    │                                                                                                                                    │
    │   Test [4510cb97] 'test_bam_sort_stats_samtools_single_end' PASSED (469.169s)                                                      │
    │   Test [ea818516] 'test_bam_sort_stats_samtools_paired_end' PASSED (18.328s)                                                       │
    │                                                                                                                                    │
    │                                                                                                                                    │
    │ SUCCESS: Executed 2 tests in 487.499s                                                                                              │
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
    │ Test Workflow BAM_SORT_STATS_SAMTOOLS                                                                                              │
    │                                                                                                                                    │
    │   Test [4510cb97] 'test_bam_sort_stats_samtools_single_end' PASSED (17.28s)                                                        │
    │   Test [ea818516] 'test_bam_sort_stats_samtools_paired_end' PASSED (16.62s)                                                        │
    │                                                                                                                                    │
    │                                                                                                                                    │
    │ SUCCESS: Executed 2 tests in 33.902s                                                                                               │
    │                                                                                                                                    │
    ╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
    INFO     All tests passed!
   ```

   :::note
   See the [nf-test docs](https://code.askimed.com/nf-test/) if you would like to run the tests manually.
   :::

5. Run [`prettier`](https://nf-co.re/docs/contributing/code_formating) on all edited and generated files
   prettier -w .

6. Check that the new subworkflow you've added follows the [new subworkflow guidelines](#new-subworkflow-guidelines-and-pr-review-checklist)

7. Lint the subworkflow locally to check that it adheres to nf-core guidelines before submission

   ```console
   $ nf-core subworkflows lint bam_sort_stats_samtools

                                         ,--./,-.
         ___     __   __   __   ___     /,-._.--~\
   |\ | |__  __ /  ` /  \ |__) |__         }  {
   | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                         `._,._,'

   nf-core/tools version 2.11 - https://nf-co.re


   INFO     Linting modules repo: '.'
   INFO     Linting subworkflow: 'bam_sort_stats_samtools'

   ╭───────────────────────╮
   │ LINT RESULTS SUMMARY  │
   ├───────────────────────┤
   │ [✔]  40 Tests Passed  │
   │ [!]   0 Test Warnings │
   │ [✗]   0 Tests Failed  │
   ╰───────────────────────╯
   ```

8. Once ready, the code can be pushed and a pull request (PR) created

   On a regular basis you can pull upstream changes into this branch and it is recommended to do so before pushing and creating a pull request - see below. Rather than merging changes directly from upstream the rebase strategy is recommended so that your changes are applied on top of the latest master branch from the nf-core repo. This can be performed as follows:

   ```bash
   git pull --rebase upstream master
   ```

   Once you are ready you can push the code and create a PR

   ```bash
   git push -u origin bam_sort_stats_samtools
   ```

   Once the PR has been accepted you should delete the branch and checkout master again.

   ```bash
   git checkout master
   git branch -d bam_sort_stats_samtools
   ```

   In case there are commits on the local branch that didn't make it into the PR (usually commits made after the PR), git will warn about this and not delete the branch. If you are sure you want to delete, use the following command

   ```bash
   git branch -D bam_sort_stats_samtools
   ```

### Publishing results

This system uses Nextflow's native [`publishDir`](https://www.nextflow.io/docs/latest/process.html#publishdir) defined directly in a pipeline workflow's `modules.config` (see [here](https://github.com/nf-core/rnaseq/blob/f7702d5b76a1351e2e7796a5ed3f59943a139fbf/conf/modules.config#L100-L106) for a simple example)

### Test data

In order to test that each subworkflow added to `nf-core/modules` is actually working and to be able to track any changes to results files between subworkflow updates we have set-up a number of Github Actions CI tests to run each subworkflow on a minimal test dataset using Docker, Singularity and Conda.

- All test data for the `nf-core/modules` repository MUST be added to the `modules` branch of [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) and organised by filename extension.

- Please adhere to the [test-data specifications](/docs/contributing/component-specifications/test_data.md) when adding new test-data

- In order to keep the size of the test data repository as minimal as possible, pre-existing files from [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) MUST be reused if at all possible.

- Test files MUST be kept as tiny as possible.

- If the appropriate test data doesn't exist in the `modules` branch of [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) please contact us on the [nf-core Slack `#subworkflows` channel](https://nfcore.slack.com/channels/subworkflows) (you can join with [this invite](https://nf-co.re/join/slack)) to discuss possible options.

- It may not be possible to add test data for some subworkflows e.g. if the input data is too large or requires a local database. In these scenarios, it is recommended to use the Nextflow [`stub`](https://www.nextflow.io/docs/latest/process.html#stub) feature to test the subworkflow. Please refer to the [`gtdbtk/classify`](https://github.com/nf-core/modules/blob/79d38a306bdaf07000e0d6f300684d3ed38c8919/modules/gtdbtk/classifywf/main.nf#L66) module and its corresponding [test script](https://github.com/nf-core/modules/blob/79d38a306bdaf07000e0d6f300684d3ed38c8919/tests/modules/gtdbtk/classifywf/main.nf#L20) to understand how to use this feature for your subworkflow development.

If a new test dataset is added to [`tests/config/test_data.config`](https://github.com/nf-core/modules/blob/master/tests/config/test_data.config), check that the config name of the added file(s) follows the scheme of the entire file name with dots replaced with underscores.

For example: the nf-core/test-datasets file `genomics/sarscov2/genome/genome.fasta` labelled as `genome_fasta`, or `genomics/sarscov2/genome/genome.fasta.fai` as `genome_fasta_fai`.

### Writing nf-test tests

We recently decided to use nf-test instead of pytest for testing modules & subworkflows. This is because nf-test is more flexible and allows us to test subworkflows in a more realistic way. You can find more information at [nf-test official docs](https://code.askimed.com/nf-test/) and [in this bytesize talk](https://nf-co.re/events/2022/bytesize_nftest).

#### Philosophy of nf-tests

- Each subworkflow contains a `tests/` folder beside the `main.nf` containing the test files
- Test files come with a [snapshot](https://code.askimed.com/nf-test/docs/assertions/snapshots/) of subworkflows output channels

#### Guidelines for creating nf-test for a subworkflow

- Some subworkflows MAY require additional parameters added to the test command to successfully run. These can be specified with an `ext.args` variable within the process scope of the `nextflow.config` file that exists alongside the test files themselves (and is automatically loaded when the test workflow `main.nf` is executed).

If your subworkflow requires a a `nextflow.config` file to run, create the file to the subworkflow's `tests/` directory and add the additional parameters there.

```bash
touch subworkflows/nf-core/bam_sort_stats_samtools/tests/nextflow.config
```

Then add the path to the `main.nf.test` file.

```groovy title="main.nf.test"
process "SUBWORKFLOW"
config "./nextflow.config"
```

- You can find examples of different nf-tests assertions on [this tutorial](https://nf-co.re/docs/contributing/tutorials/nf-test_assertions).

#### Guidelines for creating nf-test for subworkflow chained with modules

- For subworkflows that involve running a module in advance to generate required test-data, nf-test provides a [setup](https://code.askimed.com/nf-test/docs/testcases/setup/) method.

- Implementing [setup](https://code.askimed.com/nf-test/docs/testcases/setup/) with a subworkflow is very similar as with modules. For this [see docs of nf-test with chained modules](#steps-for-creating-nf-test-for-chained-modules)

### Using a stub test when required test data is too big

If the subworkflow absolutely cannot run using tiny test data, there is a possibility to add [stub-run](https://www.nextflow.io/docs/edge/process.html#stub) to the `test.yml`. In this case it is required to test the subworkflow using larger scale data and document how this is done. In addition, an extra script-block labeled `stub:` must be added, and this block must create dummy versions of all expected output files as well as the `versions.yml`. An example for modules is found in the [ascat module](https://github.com/nf-core/modules/blob/master/tests/modules/nf-core/ascat/main.nf). In the `test.yml` the `-stub-run` argument is written as well as the md5sums for each of the files that are added in the stub-block. This causes the stub-code block to be activated when the unit test is run ([example](https://github.com/nf-core/modules/blob/master/tests/modules/nf-core/ascat/test.yml)):

```bash
nextflow run tests/subworkflows/nf-core/<name_of_subworkflow> -entry test_<name_of_subworkflow> -c tests/config/nextflow.config -stub-run
```

### Migrating from pytest to nf-test

#### Steps for creating nf-test for a subworkflow

- Git checkout a new branch for your subworkflow tests

```bash
git checkout -b <branch>
```

To create the necessary files for nf-test and ensure a smooth transition, we will use the template provided by nf-core/tools.

Here are the steps to follow:

- Use nf-core/tools to create a new subworkflow with the same name as the old one with the option `--migrate-pytest`.
  This command will rename the current subworkflow directory to `<subworkflow>_old` to avoid conflicts with the new subworkflow, create a new subworkflow, and copy the `main.nf` and `meta.yml` files over to preserve the original subworkflow code.

```bash
nf-core subworkflows create <subworkflow> --migrate-pytest
```

- When using the `--migrate-pytest` option you will be asked if you want to delete the old subworkflow directory and see the content of the old pytests in the terminal, or to keep the old subworkflow directory. For the following steps, use the information from the pytest tests to create the new nf-test tests.

- Provide a test name preferably indicating the test-data and file-format used. Example: `test("homo_sapiens - [bam, bai, bed] - fasta - fai")`

:::note
Multiple tests are allowed in a single test file.
:::

- If migrating an existing subworkflow, get the inputs from current pytest files `tests/subworkflow/nf-core/subworkflow/main.nf` and provide as positional inputs `input[0]` in nf-test file

```groovy
input[0] = [
            [id:"ref"],
            file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
           ]
```

- Next, in the `then` block we can write our assertions that are used to verify the test. A test can have multiple assertions but, we recommend enclosing all assertions in a `assertAll()` block as shown below:

```groovy
assertAll(
            { assert workflow.success },
            { assert snapshot(workflow.out).match() }
          )
```

:::note
It's `workflow.` whereas with modules it's `process.`.
:::

- Run the test to create a snapshot of your subworkflow test. This will create a `main.nf.test.snap` file

```bash
nf-core subworkflows test <subworkflow>
```

:::note
The tag in `tags.yml` has to contain both `subworkflows/<subworkflow>` and not just `<subworkflow>` in contrast to modules.
:::

Time for some cleanup!
If you chose to not remove the old module directory with nf-core/tools:

- Remove the corresponding tags from `tests/config/pytest_modules.yml` so that py-tests for the subworkflow will be skipped on github CI.

- Remove the corresponding pytest files in `tests/subworkflow/nf-core`

- Remove the old subworkflow

```bash
rm -r subworkflows/nf-core/<subworkflow>_old
```

- Check if everything is according to the nf-core guidelines with:

```bash
nf-core subworkflowss lint <subworkflows>
```

- create a PR and add the `nf-test` label to it.

#### Steps for creating nf-test for subworkflow chained with modules

- Follow the steps listed above for simple subworkflows for test generation, tags and test-name

- Refer to the section [guidelines for creating nf-test for subworkflow chained with modules](guidelines-for-creating-nf-test-for-subworkflow-chained-with-modules)

:::note
Remove the corresponding tags from `tests/config/pytest_modules.yml` so that py-tests for the module will be skipped on github CI
:::

- create PR and add the `nf-test` label to it.

:::info
The implementation of nf-test in nf-core is still evolving. Things might still change and the information might here might be outdated. Please report any issues you encounter [on the nf-core/website repository](https://github.com/nf-core/website/issues/new?assignees=&labels=bug&projects=&template=bug_report.md) and the `nf-test` channel on nf-core slack.

<!-- NOTE: update when nf-core/tools gets nf-test support -->

:::

### Uploading to `nf-core/modules`

[Fork](https://help.github.com/articles/fork-a-repo/) the `nf-core/modules` repository to your own GitHub account. Within the local clone of your fork add the subworkflow files to the `subworkflows/` directory. Please try and keep PRs as atomic as possible to aid the reviewing process - ideally, one subworkflow addition/update per PR.

Commit and push these changes to your local clone on GitHub, and then [create a pull request](https://help.github.com/articles/creating-a-pull-request-from-a-fork/) on the `nf-core/modules` GitHub repo with the appropriate information.

When you are happy with your pull request, please <span class="x x-first x-last">select </span>the `Ready for Review` label on the GitHub PR tab, and providing that everything adheres to nf-core guidelines we will endeavour to approve your pull request as soon as possible. We also recommend to request reviews from the `nf-core/maintainers-team`<span class="x x-first x-last"> so </span>a core team of volunteers <span class="x x-first x-last">can try</span> to <span class="x x-first x-last">review </span>your <span class="x x-first x-last">PR</span> as fast as possible.

Once you<span class="x x-first x-last"> are </span>familiar with the subworkflow submission process, please consider joining the<span class="x x-first x-last"> reviewing</span> team by asking on the `#subworkflows` Slack channel.

### Talks

:::warning
these may include references to an older syntax, however the general idea remains the same
:::

<div class="ratio ratio-16x9">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/xuNYATGFuw4" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
</div>

<div class="ratio ratio-16x9 pb-3">
     <iframe src="https://widgets.figshare.com/articles/16825369/embed?show_title=1" width="568" height="351" allowfullscreen frameborder="0"></iframe>
</div>

### PR Review Checklist

A PR review is the process of examining a new subworkflow submission or the changes proposed to a subworkflow. The reviewer provides constructive feedback on those changes before they are merged into the nf-core repository.The goal of a PR review is to ensure that the code meets the coding standards of the project, is consistent and of high-quality.

While the team of [maintainers](https://github.com/orgs/nf-core/teams/maintainers/members) is responsible for overseeing the PR review process for subworkflows, these guidelines can assist community members in reviewing PRs and ensure that the review process is consistent and effective. The following is a collection of community suggestions to have into account during the review process.

#### General reviews of submissions to subworkflows:

In general, the main purpose of the review is to ensure

- All subworkflows adhere to the nf-core [module specifications](/docs/contributing/component-specifications/subworkflows.md)
- Ensure all checks pass, including linting, conda, singularity, and docker.

Otherwise, you can cover most of the specifications by checking for the following:

- The subworkflow is suitable for offline running, without automatic database downloads assumed.
- Check that it adheres to nf-core coding standards (e.g. use of meta map).
- Check that the code is readable and the formatting is correct (e.g. indenting, extra spaces).
- Check that there are tests for all outputs, including optional ones.
- Check that the `meta.yml` file has correct documentation links and patterns of files.
- Run the subworkflow help and check that important input (usually optional) has not been missed.
- Check that all outputs are captured by running nf-test or pytest (e.g. on Gitpod).

## Help

For further information or help, don't hesitate to get in touch on [Slack `#subworkflows` channel](https://nfcore.slack.com/channels/subworkflows) (you can join with [this invite](https://nf-co.re/join/slack)).
