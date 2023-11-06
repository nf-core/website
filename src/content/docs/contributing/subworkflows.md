---
title: DSL2 Subworkflows
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

What DSL1 users would consider an end-to-end pipeline. For example, from one or more inputs to a series of outputs. This can either be implemented using a large monolithic script as with DSL1, or by using a combination of DSL2 modules and sub-workflows. nf-core pipelines can have multiple workflows, such as processing different data types for the same ultimate purpose (such as in [nf-core/viralrecon](https://github.com/nf-core/viralrecon/tree/master/workflows))

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
3. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Conda`](https://conda.io/miniconda.html)
4. [Fork and clone the nf-core/modules repo locally](#uploading-to-nf-coremodules)
5. Set up git on your computer by adding a new git remote of the main nf-core git repo called `upstream`

   ```bash
   git remote add upstream https://github.com/nf-core/modules.git
   ```

   Make a new branch for your subworkflow and check it out

   ```bash
   git checkout -b bam_sort_stats_samtools
   ```

6. Create a subworkflow using the [nf-core DSL2 subworkflow template](https://github.com/nf-core/tools/blob/master/nf_core/subworkflow-template/subworkflows/main.nf) in the root of the clone of the nf-core/modules repository:

   ```console
   $ nf-core subworkflows create bam_sort_stats_samtools --author @joebloggs

                                             ,--./,-.
             ___     __   __   __   ___     /,-._.--~\
       |\ | |__  __ /  ` /  \ |__) |__         }  {
       | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                             `._,._,'

       nf-core/tools version 2.8 - https://nf-co.re


   INFO     Repository type: modules
   INFO     Press enter to use default values (shown in brackets) or type your own responses. ctrl+click underlined text to open links.
   INFO     Created / edited following files:
              ./subworkflows/nf-core/bam_sort_stats_samtools/main.nf
              ./subworkflows/nf-core/bam_sort_stats_samtools/meta.yml
              ./tests/subworkflows/nf-core/bam_sort_stats_samtools/main.nf
              ./tests/subworkflows/nf-core/bam_sort_stats_samtools/test.yml
              ./tests/subworkflows/nf-core/bam_sort_stats_samtools/nextflow.config
              ./tests/config/pytest_modules.yml
   ```

All of the files required to add the subworkflow to `nf-core/modules` will be created/edited in the appropriate places. There are at most 5 files to modify:

1. [`./subworkflows/nf-core/bam_sort_stats_samtools/main.nf`](https://github.com/nf-core/modules/blob/master/subworkflows/nf-core/bam_sort_stats_samtools/main.nf)

   This is the main script containing the `workflow` definition for the subworkflow. You will see an extensive number of `TODO` statements to help guide you to fill in the appropriate sections and to ensure that you adhere to the guidelines we have set for module submissions.

2. [`./subworkflows/nf-core/bam_sort_stats_samtools/meta.yml`](https://github.com/nf-core/modules/blob/master/subworkflows/nf-core/bam_sort_stats_samtools/meta.yml)

   This file will be used to store general information about the subworkflow and author details. You will need to add a brief description of the files defined in the `input` and `output` section of the main script since these will be unique to each subworkflow.

3. [`./tests/subworkflows/nf-core/bam_sort_stats_samtools/main.nf`](https://github.com/nf-core/modules/blob/master/tests/subworkflows/nf-core/bam_sort_stats_samtools/main.nf)

   Every subworkflow MUST have a test workflow. This file will define one or more Nextflow `workflow` definitions that will be used to unit test the output files created by the subworkflow. By default, one `workflow` definition will be added but please feel free to add as many as possible so we can ensure that the subworkflow works on different data types / parameters e.g. separate `workflow` for single-end and paired-end data.

   When writing multiple tests, a common practice is to alias process names to differentiate them between tests. When using an alias, add a suffix to the process name so the CI tests can still find the output in the folder named after the tool, e.g.

   ```groovy
   include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_SINGLE_END } from '../../../../subworkflows/nf-core/bam_sort_stats_samtools/main' // Good: Output folder is still 'fastqc'
   include { BAM_SORT_STATS_SAMTOOLS as SINGLE_END_BAM_SORT_STATS_SAMTOOLS } from '../../../../subworkflows/nf-core/bam_sort_stats_samtools/main' // Bad: Generates problems with CI tests - Output folder is 'post'
   ```

   Minimal test data required for your subworkflow may already exist within the [nf-core/modules repository](https://github.com/nf-core/modules/blob/master/tests/config/test_data.config), in which case you may just have to change a couple of paths in this file - see the [Test data](#test-data) section for more info and guidelines for adding new standardised data if required.

4. [`./tests/subworkflows/nf-core/bam_sort_stats_samtools/nextflow.config`](https://github.com/nf-core/modules/blob/master//tests/subworkflows/nf-core/bam_sort_stats_samtools/nextflow.config)

   Some subworkflows MAY require additional parameters added to the test command to successfully run. These can be specified with an `ext.args` variable within the process scope of the `nextflow.config` file that exists alongside the test files themselves (and is automatically loaded when the test workflow `main.nf` is executed).

5. [`./tests/subworkflows/nf-core/bam_sort_stats_samtools/test.yml`](https://github.com/nf-core/modules/blob/master/tests/subworkflows/nf-core/bam_sort_stats_samtools/test.yml)

   This file will contain all of the details required to unit test the main script in the point above using [pytest-workflow](https://pytest-workflow.readthedocs.io/). If possible, any outputs produced by the test workflow(s) MUST be included and listed in this file along with an appropriate check e.g. md5sum. The different test options are listed in the [pytest-workflow docs](https://pytest-workflow.readthedocs.io/en/stable/#test-options).

   As highlighted in the next point, we have added a command to make it much easier to test the workflow(s) defined for the subworkflow and to automatically create the `test.yml` with the md5sum hashes for all of the outputs generated by the subworkflow.

   `md5sum` checks are the preferable choice of test to determine file changes, however, this may not be possible for all outputs generated by some tools e.g. if they include time stamps or command-related headers. Please do your best to avoid just checking for the file being present e.g. it may still be possible to check that the file contains the appropriate text snippets.

6. Create a yaml file containing information required for subworkflow unit testing

   ```console
   $ nf-core subworkflows create-test-yml bam_sort_stats_samtools

                                             ,--./,-.
             ___     __   __   __   ___     /,-._.--~\
       |\ | |__  __ /  ` /  \ |__) |__         }  {
       | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                             `._,._,'

       nf-core/tools version 2.8 - https://nf-co.re


   INFO     Press enter to use default values (shown in brackets) or type your own
      responses
   Test YAML output path (- for stdout) (tests/subworkflows/nf-core/bam_sort_stats_samtools/test.yml):
   INFO     Looking for test workflow entry points: 'tests/subworkflows/nf-core/bam_sort_stats_samtools/main.nf'
   INFO     Building test meta for entry point 'test_bam_sort_stats_samtools_single_end'
   Test name (bam_sort_stats_samtools test_bam_sort_stats_samtools_single_end):
   Test command (nextflow run ./tests/subworkflows/nf-core/bam_sort_stats_samtools -entry
      test_bam_sort_stats_samtools_single_end -c ./tests/config/nextflow.config):
   Test tags (comma separated):
   Test output folder with results (leave blank to run test):
   ? Choose software profile Docker
   INFO     Setting env var '$PROFILE' to 'docker'
   INFO     Running 'bam_sort_stats_samtools' test with command:
         nextflow run ./tests/subworkflows/nf-core/bam_sort_stats_samtools -entry
            test_bam_sort_stats_samtools_single_end
            -c ./tests/config/nextflow.config --outdir /var/folders/lt/b3cs9y610fg_13q14dckwcvm0000gn/T/tmping28ow_
            -work-dir /var/folders/lt/b3cs9y610fg_13q14dckwcvm0000gn/T/tmportf0uab
   ```

   :::note
   See docs for [running tests manually](#running-tests-manually) if you would like to run the tests manually.
   :::

7. Run [`prettier`](https://nf-co.re/docs/contributing/code_formating) on all edited and generated files
   prettier -w .

8. Check that the new subworkflow you've added follows the [new subworkflow guidelines](#new-subworkflow-guidelines-and-pr-review-checklist)

<!-- TODO: nf-core: Update these guidelines as we develop them -->

    (COMMAND NOT IMPLEMENTED IN NF-CORE/TOOLS YET!!) Lint the subworkflow locally to check that it adheres to nf-core guidelines before submission

    <!-- TODO: nf-core: Update these guidelines as we develop them -->

    ```console
    $ nf-core subworkflows lint

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.7.dev0 - https://nf-co.re


    INFO     Repository type: modules
    ```

9. Once ready, the code can be pushed and a pull request (PR) created

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

### Test data

In order to test that each subworkflow added to `nf-core/modules` is actually working and to be able to track any changes to results files between subworkflow updates we have set-up a number of Github Actions CI tests to run each subworkflow on a minimal test dataset using Docker, Singularity and Conda.

- All test data for the `nf-core/modules` repository MUST be added to the `modules` branch of [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) and organised by filename extension.

- In order to keep the size of the test data repository as minimal as possible, pre-existing files from [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) MUST be reused if at all possible.

- Test files MUST be kept as tiny as possible.

- If the appropriate test data doesn't exist in the `modules` branch of [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) please contact us on the [nf-core Slack `#subworkflows` channel](https://nfcore.slack.com/channels/subworkflows) (you can join with [this invite](https://nf-co.re/join/slack)) to discuss possible options.

- It may not be possible to add test data for some subworkflows e.g. if the input data is too large or requires a local database. In these scenarios, it is recommended to use the Nextflow [`stub`](https://www.nextflow.io/docs/latest/process.html#stub) feature to test the subworkflow. Please refer to the [`gtdbtk/classify`](https://github.com/nf-core/modules/blob/79d38a306bdaf07000e0d6f300684d3ed38c8919/modules/gtdbtk/classifywf/main.nf#L66) module and its corresponding [test script](https://github.com/nf-core/modules/blob/79d38a306bdaf07000e0d6f300684d3ed38c8919/tests/modules/gtdbtk/classifywf/main.nf#L20) to understand how to use this feature for your subworkflow development.

### Running tests manually

As outlined in the [nf-core subworkflows create](#nf-core-subworkflows-create) section we have made it quite trivial to create an initial yaml file (via the `nf-core subworkflows create-test-yml` command) containing a listing of all of the subworkflow output files and their associated md5sums. However, md5sum checks may not be appropriate for all output files if for example they contain timestamps. This is why it is a good idea to re-run the tests locally with `pytest-workflow` before you create your pull request adding the subworkflow. If your files do indeed have timestamps or other issues that prevent you from using the md5sum check, then you can edit the `test.yml` file to instead check that the file contains some specific content or as a last resort, if it exists. The different test options are listed in the [pytest-workflow docs](https://pytest-workflow.readthedocs.io/en/stable/#test-options).

Please follow the steps below to run the tests locally:

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Conda`](https://conda.io/miniconda.html)

3. Install [`pytest-workflow`](https://pytest-workflow.readthedocs.io/en/stable/#installation)

4. Start running your own tests using the appropriate [`tag`](https://github.com/nf-core/modules/blob/20d8250d9f39ddb05dfb437603aaf99b5c0b2b41/tests/modules/fastqc/test.yml) defined in the `test.yml`:

   - Run the test with the helper tool `nf-core subworkflows test` from the modules directory.

    <!-- TODO: nf-core: Update these guidelines as we develop them -->

   ```console
   $ cd /path/to/git/clone/of/nf-core/modules/
   $ nf-core subworkflows test bam_sort_stats_samtools

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 2.8 - https://nf-co.re


    INFO     Press enter to use default values (shown in brackets) or type your own responses
    ? Choose software profile Docker
    INFO     Setting environment variable '$PROFILE' to 'conda'
    NFO     Running pytest for subworkflow 'bam_sort_stats_samtools'
   ```

   - See [docs on running pytest-workflow](https://pytest-workflow.readthedocs.io/en/stable/#running-pytest-workflow) for more info.

:::note
For docker/singularity, setting the environment variable `TMPDIR=~` is an example of a location the containers can mount (you can change this as you prefer). If you get test failures such as with Nextflow errors that end in `work doesn't exist in container`, check your container can mount your `TMPDIR`.
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

## New subworkflow guidelines and PR review checklist

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

### General

1. Subworkflows should combine tools that make up a logical unit in an analysis step. A subworkflow must contain at least two modules.

2. Each `subworkflow` emits a channel containing all `versions.yml` collecting the tool(s) versions. They MUST be collected within the workflow and added to the output as `versions` :

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

1. The directory structure for the subworkflow name must be all lowercase e.g. [`subworkflows/nf-core/bam_sort_stats_samtools/`](https://github.com/nf-core/modules/tree/master/subworkflows/nf-core/bam_sort_stats_samtools/). The naming convention should be of the format `<file_type>_<operation_1>_<operation_n>_<tool_1>_<tool_n>` e.g. `bam_sort_stats_samtools` where `bam` = `<file_type>`, `sort` = `<operation>` and `samtools` = `<tool>`. Not all operations are required in the name if they are routine (e.g. indexing after creation of a BAM). Operations can be collapsed to a general name if the steps are directly related to each other. For example if in a subworkflow, a binning tool has three required steps (e.g. `<tool> split`, `<tool> calculate`, `<tool> merge`) to perform an operation (contig binning) these can be collapsed into one (e.g. `fasta_binning_concoct`, rather than `fasta_split_calculate_merge_concoct`). If in doubt regarding what to name your subworkflow, please contact us on the [nf-core Slack `#subworkflows` channel](https://nfcore.slack.com/channels/subworkflows) (you can join with [this invite](https://nf-co.re/join/slack)) to discuss possible options.

2. All parameter names MUST follow the `snake_case` convention.

3. All function names MUST follow the `camelCase` convention.

4. Channel names MUST follow `snake_case` convention and be all lower case.

5. Input channel names SHOULD signify the input object type. For example, a single value input channel will be prefixed with `val_`, whereas input channels with multiple elements (e.g. meta map + file) should be prefixed with `ch_`.

6. Output channel names SHOULD only be named based on the major output file of that channel (i.e, an output channel of `[[meta], bam]` should be emitted as `bam`, not `ch_bam`). This is for more intuitive use of these output objects downstream with the `.out` attribute.

### Input/output options

1. Input channel declarations MUST be defined for all _possible_ input files that will be required by the subworkflow (i.e. both required and optional files) within the `take` block.

2. Named file extensions MUST be emitted for ALL output channels e.g. `path "*.txt", emit: txt`.

3. Optional inputs are not currently supported by Nextflow. However, passing an empty list (`[]`) instead of a file as a subworkflow parameter can be used to work around this issue.

### Subworkflow parameters

1. Named `params` defined in the parent workflow MUST NOT be assumed to be passed to the subworkflow to allow developers to call their parameters whatever they want. In general, it may be more suitable to use additional `input` value channels to cater for such scenarios.

### Documentation

1. Each input and output channel SHOULD have a comment describing the output structure of the channel e.g

   ```nextflow
   input:
   ch_reads // channel: [mandatory] meta, reads
   val_sort // boolean: [mandatory] false
   <...>

   emit:
   bam = SAMTOOLS_VIEW.out.bam // channel: [ val(meta), path(bam) ]
   versions = ch_versions      // channel: [ path(versions.yml) ]
   ```

2. Each input and output channel structure SHOULD also be described in the `meta.yml` in the description entry.

   ```text
   description: |
     Structure: [ val(meta), path(tsv)]
     (Sub)contig coverage table
   ```

### Publishing results

This system uses Nextflow's native [`publishDir`](https://www.nextflow.io/docs/latest/process.html#publishdir) defined directly in a pipeline workflow's `modules.config` (see [here](https://github.com/nf-core/rnaseq/blob/f7702d5b76a1351e2e7796a5ed3f59943a139fbf/conf/modules.config#L100-L106) for a simple example)

### Test data config file

If a new test dataset is added to [`tests/config/test_data.config`](https://github.com/nf-core/modules/blob/master/tests/config/test_data.config), check that the config name of the added file(s) follows the scheme of the entire file name with dots replaced with underscores.

For example: the nf-core/test-datasets file `genomics/sarscov2/genome/genome.fasta` labelled as `genome_fasta`, or `genomics/sarscov2/genome/genome.fasta.fai` as `genome_fasta_fai`.

### Migrating from pytest to nf-test

We recently decided to use nf-test instead of pytest for testing modules & subworkflows. This is because nf-test is more flexible and allows us to test subworkflows in a more realistic way. You can find more information at [nf-test official docs](https://code.askimed.com/nf-test/) and [in this bytesize talk](https://nf-co.re/events/2022/bytesize_nftest).

#### Philosophy of nf-tests

- Each subworkflow contains a `tests/` folder beside the `main.nf` containing the test files
- Test files come with a [snapshot](https://code.askimed.com/nf-test/docs/assertions/snapshots/) of subworkflows output channels

#### Steps for creating nf-test for a subworkflow

- Install and update to the latest version of [nf-test](https://code.askimed.com/nf-test/installation/)

```bash
conda install -c bioconda nf-test
```

- Git checkout a new branch for your subworkflow tests

```bash
git checkout -b <branch>
```

- Create a new tests directory within your subworkflow directory

```bash
mkdir subworkflows/nf-core/<subworkflow>/tests
```

- Generate a test file from template for your subworkflow using nf-test

```bash
nf-test generate workflow subworkflows/nf-core/<subworkflow>/main.nf
```

:::note
we use `workflow` instead of `process` as it is a subworkflow not a module
:::

- Move the generated test file to the tests directory

```bash
mv subworkflows/nf-core/<subworkflows>/main.nf.test subworkflows/nf-core/<subworkflows>/tests/
```

- Open the `main.nf.test` file and change the path for the script to a relative path `../main.nf`

```diff title="main.nf.test"
nextflow_workflow {
     name "Test Workflow SUBWORKFLOW"
-    script "subworkflows/nf-core/<subworkflow>/main.nf"
+    script "../main.nf"
     process "SUBWORKFLOW"
```

- Then add tags to identify this module

```groovy title="main.nf.test"
tag "subworkflows"
tag "subworkflows_nfcore"
tag "subworkflows/<subworkflow>"
tag "<subworkflow>"
tag "<tool1>"
tag "<tool2>"
```

:::note
We require to have `tag subworkflows/<subworkflow>` as so it's picked up correctly
:::

:::note
Include the used tools in the tags
:::

- Provide a test name preferably indicating the test-data and file-format used. Example: `test("homo_sapiens - [bam, bai, bed] - fasta - fai")`

:::note
multiple tests are allowed in a single test file
:::

- If migrating an existing subworkflow, get the inputs from current pytest files `tests/subworkflow/nf-core/subworkflow/main.nf` and provide as positional inputs `input[0]` in nf-test file

```groovy
input[2] = [
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
It's `workflow.` and whereas with modules it's `process.`.
:::

- Run the test to create a snapshot of your module test. This will create a `.nf.test.snap` file

```bash
nf-test test --tag "<subworkflow>" --profile docker --update-snapshot
```

- Re-run the test again to verify if snapshots match

```bash
nf-test test --tag "<subworkflow>" --profile docker
```

- Create a new `tags.yml` in the `subworkflows/nf-core/<subworkflow>/tests/` folder and add only the corresponding subworkflow tag from `tests/config/pytest_modules.yml`

```yaml
subworkflows/<subworkflow>:
  - subworkflows/nf-core/<subworkflow>/**
```

:::note
Remove the corresponding tags from `tests/config/pytest_modules.yml` so that py-tests for the module will be skipped on github CI
:::

:::note
The tag has to contain both `subworkflows/<subworkflow>` and not just `<subworkflow>` as supposed with modules
:::

- create PR and add the `nf-test` label to it.

#### Steps for creating nf-test for subworkflow chained with modules

- Follow the steps listed above for simple subworkflows for test generation, tags and test-name

- For subworkflows that involve running a module in advance to generate required test-data, nf-test provides a [setup](https://code.askimed.com/nf-test/docs/testcases/setup/) method.

- Implementing [setup](https://code.askimed.com/nf-test/docs/testcases/setup/) with a subworkflow is very similar as with modules. For this [see docs of nf-test with chained modules](./modules#steps-for-creating-nf-test-for-chained-modules)

:::note
Remove the corresponding tags from `tests/config/pytest_modules.yml` so that py-tests for the module will be skipped on github CI
:::

- create PR and add the `nf-test` label to it.

:::info
The implementation of nf-test in nf-core is still evolving. Things might still change and the information might here might be outdated. Please report any issues you encounter [on the nf-core/website repository](https://github.com/nf-core/website/issues/new?assignees=&labels=bug&projects=&template=bug_report.md) and the `nf-test` channel on nf-core slack. Additionally, nf-core/tools will help you create nf-tests in the future, making some of the steps here obsolete.

<!-- NOTE: update when nf-core/tools gets nf-test support -->

:::

### Using a stub test when required test data is too big

If the subworkflow absolutely cannot run using tiny test data, there is a possibility to add [stub-run](https://www.nextflow.io/docs/edge/process.html#stub) to the `test.yml`. In this case it is required to test the subworkflow using larger scale data and document how this is done. In addition, an extra script-block labeled `stub:` must be added, and this block must create dummy versions of all expected output files as well as the `versions.yml`. An example for modules is found in the [ascat module](https://github.com/nf-core/modules/blob/master/tests/modules/nf-core/ascat/main.nf). In the `test.yml` the `-stub-run` argument is written as well as the md5sums for each of the files that are added in the stub-block. This causes the stub-code block to be activated when the unit test is run ([example](https://github.com/nf-core/modules/blob/master/tests/modules/nf-core/ascat/test.yml)):

```bash
nextflow run tests/subworkflows/nf-core/<name_of_subworkflow> -entry test_<name_of_subworkflow> -c tests/config/nextflow.config -stub-run
```

## What is the `meta` map?

In nf-core DSL2 pipelines, to add sample-specific information and metadata that is carried throughout the pipeline, we use a meta variable. This avoids the need to create separate channels for each new characteristic.
The meta variable can be passed down to processes as a tuple of the channel containing the actual samples, e.g. FastQ files, and the meta variable. The `meta map` is a [groovy map](https://www.tutorialspoint.com/groovy/groovy_maps.htm), which is like a python dictionary.

<!-- TODO: nf-core: Link to DSL2 modules docs section for this instead of duplicating here -->

## Help

For further information or help, don't hesitate to get in touch on [Slack `#subworkflows` channel](https://nfcore.slack.com/channels/subworkflows) (you can join with [this invite](https://nf-co.re/join/slack)).
