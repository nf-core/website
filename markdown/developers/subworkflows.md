---
title: DSL2 Subworkflows
subtitle: Guidelines and reference for DSL2 subworkflows
---

If you decide to upload a subworkflows to `nf-core/modules` then this will ensure that it will become available to all nf-core pipelines, and to everyone within the Nextflow community! See [`subworkflows/`](https://github.com/nf-core/modules/tree/master/subworkflows) for examples.

_TODO_ strip this out of modules and here and move to shared documents
<!-- ## Terminology

The features offered by Nextflow DSL2 can be used in various ways depending on the granularity with which you would like to write pipelines. Please see the listing below for the hierarchy and associated terminology we have decided to use when referring to DSL2 components.

### Module

A `process` that can be used within different pipelines and is as atomic as possible i.e. cannot be split into another module. An example of this would be a module file containing the process definition for a single tool such as `FastQC`. At present, this repository hosts atomic module files that should be added to the [`modules/`](https://github.com/nf-core/modules/tree/master/modules) directory of nf-core/modules along with the required documentation and tests.

### Sub-workflow

A chain of at least two modules that offer a higher-level of functionality within the context of a pipeline. For example, a sub-workflow to run multiple QC tools with FastQ files as input. Sub-workflows should be shipped with the pipeline implementation and if required they should be shared amongst different pipelines directly from there. At present, this repository hosts subworkflows files that should be added to the [`subworkflows/`](https://github.com/nf-core/modules/tree/master/subworkflows) directory of nf-core/modules along with the required documentation and tests.

### Workflow

What DSL1 users would consider an end-to-end pipeline. For example, from one or more inputs to a series of outputs. This can either be implemented using a large monolithic script as with DSL1, or by using a combination of DSL2 individual modules and sub-workflows.

----------------------------------------------------------------

-->
## Writing a new subworkflow

### Before you start

Please check that the subworkflow you wish to add isn't already on [`nf-core/modules`](https://github.com/nf-core/modules/tree/master/subworkflows):

- Use the [`nf-core subworkflow list`](https://github.com/nf-core/tools#list-subworkflow) command
- Check [open pull requests](https://github.com/nf-core/modules/pulls)
- Search [open issues](https://github.com/nf-core/modules/issues)

If the subworkflow doesn't exist on `nf-core/modules`:

- Please create a [new issue](https://github.com/nf-core/modules/issues/new?assignees=&labels=new%20module&template=new_nodule.md&title=new%20module:) before adding it
- Set an appropriate subject for the issue e.g. `new subworkflow: bam_sort_samtools`
- Add yourself to the `Assignees` so we can track who is working on the subworkflow

### New subworkflow workflow

We have implemented a number of commands in the `nf-core/tools` package to make it incredibly easy for you to create and contribute your own subworkflows to nf-core/modules.

1. Install the latest version of [`nf-core/tools`](https://github.com/nf-core/tools#installation) (`>=2.5.2`)
2. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)
3. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Conda`](https://conda.io/miniconda.html)
4. [Fork and clone the nf-core/modules repo locally](#uploading-to-nf-coremodules)
5. Set up git on your computer by adding a new git remote of the main nf-core git repo called `upstream`

   ```bash
   git remote add upstream https://github.com/nf-core/modules.git
   ```

   Make a new branch for your subworkflow and check it out

   ```bash
   git checkout -b bam_sort_samtools
   ```

6. Create a subworkflow using the [nf-core DSL2 subworkflows template](https://github.com/nf-core/tools/blob/master/nf_core/module-template/modules/main.nf):

TODO: check link, update command line output

   ```console
   $ nf-core subworkflows create bam_sort_samtools --author @joebloggs --label process_low --meta


   ```

   All of the files required to add the module to `nf-core/modules` will be created/edited in the appropriate places. There are at most 5 files to modify:

   1. [`./subworkflow/bam_sort_samtools/main.nf`](https://github.com/nf-core/modules/blob/master/subworkflows/bam_sort_samtools/main.nf)

      This is the main script containing the `workflow` definition for the subworkflow. You will see an extensive number of `TODO` statements to help guide you to fill in the appropriate sections and to ensure that you adhere to the guidelines we have set for subworkflow submissions.

   2. [`./subworkflow/bam_sort_samtools/meta.yml`](https://github.com/nf-core/modules/blob/master/subworkflow/bam_sort_samtools/meta.yml)

      This file will be used to store general information about the subworkflow and author details - the majority of which will already be auto-filled. However, you will need to add a brief description of the files defined in the `input` and `output` section of the main script since these will be unique to each subworkflow.

   3. [`./tests/subworkflow/bam_sort_samtools/main.nf`](https://github.com/nf-core/modules/blob/master/tests/subworkflows/bam_sort_samtools/main.nf)

      Every subworkflow MUST have a test workflow. This file will define one or more Nextflow `workflow` definitions that will be used to unit test the output files created by the subworkflow. By default, one `workflow` definition will be added but please feel free to add as many as possible so we can ensure that the module works on different data types / parameters e.g. separate `workflow` for single-end and paired-end data.

      When writing multiple tests, a common practice is to alias subworkflow names to differentiate them between tests. When using an alias, add a suffix to the process name so the CI tests can still find the output in the folder named after the tool, e.g.

      ```groovy
      include { BAM_SORT_SAMTOOLS as BAM_SORT_SAMTOOLS_POST } from '../../.....' // Good: Output folder is still 'fastqc'
      include { BAM_SORT_SAMTOOLS as POST_BAM_SORT_SAMTOOLS    } from '../../.....' // Bad: Generates problems with CI tests - Output folder is 'post'
      ```

      Minimal test data required for your subworkflow may already exist within the [nf-core/modules repository](https://github.com/nf-core/modules/blob/master/tests/config/test_data.config), in which case you may just have to change a couple of paths in this file - see the [Test data](#test-data) section for more info and guidelines for adding new standardised data if required.

   4. [`./tests/subworkflow/bam_sort_samtools/nextflow.config`](https://github.com/nf-core/modules/blob/master/tests/subworkflow/amps/nextflow.config)

      Some modules MAY require additional parameters added to the test command to successfully run. These can be specified with an `ext.args` variable within the process scope of the `nextflow.config` file that exists alongside the test files themselves (and is automatically loaded when the test workflow `main.nf` is executed).

   5. [`./tests/subworkflows/bam_sort_samtools/test.yml`](https://github.com/nf-core/modules/blob/master/tests/modules/fastqc/test.yml)

      This file will contain all of the details required to unit test the main script in the point above using [pytest-workflow](https://pytest-workflow.readthedocs.io/). If possible, any outputs produced by the test workflow(s) MUST be included and listed in this file along with an appropriate check e.g. md5sum. The different test options are listed in the [pytest-workflow docs](https://pytest-workflow.readthedocs.io/en/stable/#test-options).

TODO: is this command existing pipeline

      As highlighted in the next point, we have added a command to make it much easier to test the workflow(s) defined for the module and to automatically create the `test.yml` with the md5sum hashes for all of the outputs generated by the subworkflow.

      `md5sum` checks are the preferable choice of test to determine file changes, however, this may not be possible for all outputs generated by some tools e.g. if they include time stamps or command-related headers. Please do your best to avoid just checking for the file being present e.g. it may still be possible to check that the file contains the appropriate text snippets.

7. Create a yaml file containing information required for module unit testing

TODO: update command prompt

   ```console
   $ nf-core subworkflows create-test-yml


   ```

   > NB: See docs for [running tests manually](#running-tests-manually) if you would like to run the tests manually.

8. Check that the new subworkflow you've added follows the [new subworkflow guidelines](#new-subworkflow-guidelines-and-pr-review-checklist)

9. Lint the subworkflo locally to check that it adheres to nf-core guidelines before submission

TODO: update prompt

   ```console
   $ nf-core subworkflows lint bam_sort_samtools --dir .


   ```

10. Once ready, the code can be pushed and a pull request (PR) created

    On a regular basis you can pull upstream changes into this branch and it is recommended to do so before pushing and creating a pull request - see below. Rather than merging changes directly from upstream the rebase strategy is recommended so that your changes are applied on top of the latest master branch from the nf-core repo. This can be performed as follows

```bash
git pull --rebase upstream master
```

Once you are ready you can push the code and create a PR

```bash
git push -u origin bam_sort_samtools
```

Once the PR has been accepted you should delete the branch and checkout master again.

```bash
git checkout master
git branch -d bam_sort_samtools
```

In case there are commits on the local branch that didn't make it into the PR (usually commits made after the PR), git will warn about this and not delete the branch. If you are sure you want to delete, use the following command

```bash
git branch -D bam_sort_samtools
```

### Test data

In order to test that each subworkflow added to `nf-core/modules` is actually working and to be able to track any changes to results files between subworkflow and upstream module updates we have set-up a number of Github Actions CI tests to run each subworkflow on a minimal test dataset using Docker, Singularity and Conda.

- All test data for `nf-core/modules` MUST be added to the `modules` branch of [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) and organised by filename extension.

- In order to keep the size of the test data repository as minimal as possible, pre-existing files from [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) MUST be reused if at all possible.

- Test files MUST be kept as tiny as possible.

- If the appropriate test data doesn't exist in the `modules` branch of [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) please contact us on the [nf-core Slack `#modules` channel](https://nfcore.slack.com/channels/modules) (you can join with [this invite](https://nf-co.re/join/slack)) to discuss possible options.

- It may not be possible to add test data for some modules e.g. if the input data is too large or requires a local database. In these scenarios, it is recommended to use the Nextflow [`stub`](https://www.nextflow.io/docs/latest/process.html#stub) feature to test the module. Please refer to the [`gtdbtk/classify`](https://github.com/nf-core/modules/blob/79d38a306bdaf07000e0d6f300684d3ed38c8919/modules/gtdbtk/classifywf/main.nf#L66) module and its corresponding [test script](https://github.com/nf-core/modules/blob/79d38a306bdaf07000e0d6f300684d3ed38c8919/tests/modules/gtdbtk/classifywf/main.nf#L20) to understand how to use this feature for your module development.

### Running tests manually

As outlined in the [nf-core subworkflows create](#nf-core-subworkflows-create) section we have made it quite trivial to create an initial yaml file (via the `nf-core subworkflows create-test-yml` command) containing a listing of all of the subworkflow output files and their associated md5sums. However, md5sum checks may not be appropriate for all output files if for example they contain timestamps. This is why it is a good idea to re-run the tests locally with `pytest-workflow` before you create your pull request adding the module. If your files do indeed have timestamps or other issues that prevent you from using the md5sum check, then you can edit the `test.yml` file to instead check that the file contains some specific content or as a last resort, if it exists. The different test options are listed in the [pytest-workflow docs](https://pytest-workflow.readthedocs.io/en/stable/#test-options).

Please follow the steps below to run the tests locally:

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Conda`](https://conda.io/miniconda.html)

3. Install [`pytest-workflow`](https://pytest-workflow.readthedocs.io/en/stable/#installation)

4. Start running your own tests using the appropriate [`tag`](https://github.com/nf-core/modules/blob/20d8250d9f39ddb05dfb437603aaf99b5c0b2b41/tests/subworkflows/bam_sort_samtools/test.yml) defined in the `test.yml`:

TODO: update output

   - Run the test with the helper tool `nf-core subworkflows test` from the modules directory.

     ```console
     $ cd /path/to/git/clone/of/nf-core/modules/
     $ nf-core subworkflows test bam_sort_samtools

     ```

   - See [docs on running pytest-workflow](https://pytest-workflow.readthedocs.io/en/stable/#running-pytest-workflow) for more info.

> 🛈 For docker/singularity, setting the environment variable `TMPDIR=~` is an example of a location the containers can mount (you can change this as you prefer). If you get test failures such as with Nextflow errors that end in `work doesn't exist in container`, check your container can mount your `TMPDIR`.
>
> :warning: if you have a module named `build` this can conflict with some pytest internal behaviour. This results in no tests being run (i.e. recieving a message of `collected 0 items`). In this case rename the `tests/<module>/build` directory to `tests/<module>/build_test`, and update the corresponding `test.yml` accordingly. An example can be seen with the [`bowtie2/build` module tests](https://github.com/nf-core/modules/tree/master/tests/modules/bowtie2/build_test).

### Uploading to `nf-core/modules`

[Fork](https://help.github.com/articles/fork-a-repo/) the `nf-core/modules` repository to your own GitHub account. Within the local clone of your fork add the subworkflow file to the `subworkflows/` directory. Please try and keep PRs as atomic as possible to aid the reviewing process - ideally, one subworkflow addition/update per PR.

Commit and push these changes to your local clone on GitHub, and then [create a pull request](https://help.github.com/articles/creating-a-pull-request-from-a-fork/) on the `nf-core/modules` GitHub repo with the appropriate information.

When you are happy with your pull request, please <span class="x x-first x-last">select </span>the `Ready for Review` label on the GitHub PR tab, and providing that everything adheres to nf-core guidelines we will endeavour to approve your pull request as soon as possible. We also recommend to request reviews from the `nf-core/modules-team`<span class="x x-first x-last"> so </span>a core team of volunteers <span class="x x-first x-last">can try</span> to <span class="x x-first x-last">review </span>your <span class="x x-first x-last">PR</span> as fast as possible.

Once you<span class="x x-first x-last"> are </span>familiar with the module submission process, please consider joining the<span class="x x-first x-last"> reviewing</span> team by asking on the `#subworkflows` slack channel.

## New subworkflow guidelines and PR review checklist

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

### General

1. Subworkflows should combine tools that make up a logical unit in an analysis step. .... how many tools are appropriate?

2.  Each `module` emits a channel containing `versions.yml` collecting the tool(s) versions. They MUST be collected within the workflow and added to the output as `versions` :

```bash

take:
  input

main:

  ch_versions = Channel.empty()

  FASTQC(input)

  ch_versions = ch_versions.mix(FASTQC.out.versions())

emit:
  versions = ch_versions
```

### Naming conventions

1. The directory structure for the subworkflow name must be all lowercase e.g. [`subworkflows/bam_sort_samtools/`](https://github.com/nf-core/modules/tree/master/subworkflows/bam_sort_samtools/).

2. The subworkflow name in the subworkflow file MUST be all uppercase e.g. `process BAM_SORT_SAMTOOLS {`. The name of the file type (i.e. `BAM`), the action (i.e. `SORT`), and the main tool (i.e. `SAMTOOLS`) MUST be all one word separated by an underscore.

3. Output channel names SHOULD just consist of the file-format suffix (e.g. `vcf` or `bam`).

### Input/output options

1. Input channel declarations MUST be defined for all _possible_ input files (i.e. both required and optional files).

   - Directly associated auxiliary files to an input file MAY be defined within the same input channel alongside the main input channel (e.g. [BAM and BAI](https://github.com/nf-core/modules/blob/e937c7950af70930d1f34bb961403d9d2aa81c7d/modules/samtools/flagstat/main.nf#L22)).
   - Other generic auxiliary files used across different input files (e.g. common reference sequences) MAY be defined using a dedicated input channel (e.g. [reference files](https://github.com/nf-core/modules/blob/3cabc95d0ed8a5a4e07b8f9b1d1f7ff9a70f61e1/modules/bwa/mem/main.nf#L21-L23)).

2. Named output declarations MUST be emitted for ALL output channels.

3. Optional inputs are not currently supported by Nextflow. However, passing an empty list (`[]`) instead of a file as a subworkflow parameter can be used to work around this issue.

### Subworkflow parameters

1. A subworkflow file SHOULD only define input and output files as command-line parameters to be executed within the process.

2. All `params` that may determine whether a module is executed in a workflow MUST be handled using the `ext.when` directive of the respective module. `params` SHOULD NOT be passed into the workflow to allow developers full flexibility in parameter naming on pipeline level.
TODO: possibly very ugly since nextflow will print a bunch of warnings, but I am not in favor of passing parameters around.

### Test data config file

If a new test dataset is added to [`tests/config/test_data.config`](https://github.com/nf-core/modules/blob/master/tests/config/test_data.config), check that the config name of the added file(s) follows the scheme of the entire file name with dots replaced with underscores.

For example: the nf-core/test-datasets file `genomics/sarscov2/genome/genome.fasta` labelled as `genome_fasta`, or `genomics/sarscov2/genome/genome.fasta.fai` as `genome_fasta_fai`.

### Using a stub test when required test data is too big

If the subworkflow absolute cannot run using tiny test data, there is a possibility to add [stub-run](https://www.nextflow.io/docs/edge/process.html#stub) to the test.yml. In this case it is required to test the subworkflow using larger scale data and document how this is done. In addition, an extra script-block labeled `stub:` must be added, and this block must create dummy versions of all expected output files as well as the `versions.yml`. An example is found in the [ascat module](https://github.com/nf-core/modules/blob/master/modules/ascat/main.nf). In the `test.yml` the `-stub-run` argument is written as well as the md5sums for each of the files that are added in the stub-block. This causes the stub-code block to be activated when the unit test is run ([example](https://github.com/nf-core/modules/blob/master/tests/modules/ascat/test.yml)):

```console
nextflow run tests/subworkflows/<nameofsubworkflow> -entry test_<nameofsubworkflow> -c tests/config/nextflow.config -stub-run
```
## Help

For further information or help, don't hesitate to get in touch on [Slack `#subworkflows` channel](https://nfcore.slack.com/channels/subworkflows) (you can join with [this invite](https://nf-co.re/join/slack)).
