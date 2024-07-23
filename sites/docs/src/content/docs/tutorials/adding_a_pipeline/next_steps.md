---
title: Next steps
subtitle: Follow this walkthrough to add a new pipeline to nf-core.
weight: 50
type: 'tutorial'
---

### Adding new pipeline features to existing pipelines

<!-- TODO: Delete / link to another section of the docs -->

nf-core is an open and inclusive community, and welcomes contributions to current pipelines. In many cases, the original developers might not have experience with some new fancy method you would like to add or simply doesn't have the time to implement everything themselves - so they might be really happy to have you actively contributing!

Basic rules for contributions:

- Ask in the pipeline channel whether there is an open issue on the respective pipeline's issue tracker for the feature you're planning to implement
    - If not, create a new issue describing your ideas and wait for someone to comment/discuss
- If there is community consensus, start implementing the feature in your [fork](https://help.github.com/en/articles/fork-a-repo) of the respective pipeline

:::note
Please do not write to multiple channels in the Slack community. Rather, collect all of the information in a single GitHub issue, which makes it also much easier to follow up on your proposal
:::

## Continuous integration testing

<!-- TODO: Delete / link to another section of the docs -->

Automated continuous integration (CI) testing is used to check pipelines are not broken after a change has been pushed to GitHub. The CI testing is initiated via a series of GitHub actions that are defined in the `.github/workflows` directory. The CI testing utilizes  two main test profiles (`conf/test.config` and `conf/test_full.config`). However, a pipeline can may require several profiles for testing different parts of the pipeline and can be created and implemented at the developers discretion. Find more information about how to set-up the test workflow for your pipeline [here](/docs/tutorials/adding_a_pipeline/test_data).

## DSL2 and modules

<!-- TODO: Delete / link to another section of the docs -->

Nextflow DSL2 allows for a more modularized design of pipelines and the reuse of components. nf-core has developed design best practices for DSL2 pipelines which should be used by all nf-core pipelines in order to assure standardization and the reuse of components. The following descriptions describe certain design standards for how a nf-core DSL2 pipeline should be build.

### Modules

Each nf-core pipeline should have a `modules` directory which contains the code for all modules. A module depicts a single process which involves, if possible, only a single tool/software. The `modules` directory is furthermore divided into `local`and `nf-core` subdirectories. The `local` subdirectory contains modules that are specific to the pipeline. The `nf-core` subdirectory contains modules that are installed from the `nf-core/modules` repository. For instance, most pipelines that involve FastQ files will run the FastQC tool for quality control. The module needed for this can be easily reused from the `nf-core/modules` directory using the `nf-core/tools`package.

More information and a comprehensive guide of how to implement modules in a pipeline can be found [here](https://nf-co.re/docs/contributing/modules).

### Sample meta information

In nf-core DSL2 pipelines, every channel that contains sample data in some way should also contain a `meta`variable, which must contain the fields `meta.id`, `meta.single_end` and `meta.strandedness`. The `meta` variable can be passed down to processes as a tuple of the channel containing the actual samples, e.g. FastQ files, and the `meta` variable. This meta information can easily be extracted from a samplesheet which specifies the input files. Use the Nextflow plugin [nf-validation](https://nextflow-io.github.io/nf-validation/) to transform the samplesheet into a channel with the meta information.
