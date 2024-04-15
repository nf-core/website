---
title: Appendices
subtitle: Follow this walkthrough to add a new pipeline to nf-core.
weight: 50
type: 'tutorial'
---

### Adding new pipeline features to existing pipelines

<!-- TODO: Delete / link to another section of the docs -->

We are an open and inclusive community, welcoming any contributions to pipelines already present in nf-core. In many cases, the original developers might either not have experience with some new fancy method or simply doesn't have the time to implement everything themselves - so they might be really happy to see you actively contributing!

Basic rules for such contributions:

- Ask in the [Slack](https://nf-co.re/join/slack) channel for the specific pipeline whether there is an open issue on the respective pipeline's issue tracker for the feature you're planning to
- If not, create a new issue there, describing the purpose and ideas you have and wait for someone to comment/discuss
- If everyone is happy or there is some consensus in the community, start implementing the feature in your [fork](https://help.github.com/en/articles/fork-a-repo) of the respective pipeline
- Please do not write to multiple channels in the Slack community, rather collect all of the information in a single GitHub issue, which makes it also much easier to follow up on your proposal

### Adding new dependencies to an existing pipeline

<!-- TODO: Delete / link to another section of the docs -->

Sometimes, especially when adding new features to a pipeline, the dependencies change as well. In such cases, you might want to have an updated Docker Container available before submitting a pull request, in order to have the GitHub Actions tests run through when testing your updated code. To achieve that, please follow these steps:

- Add _only_ the newly required dependencies to the `environment.yml` in the pipeline code
- If you only add new processes to an already existing pipeline however, you can simply specify the container in the `nextflow.config` file, like so:

```nextflow
process {
    withName:foo {
        container = 'image_name_1'
    }
    withName:bar {
        container = 'image_name_2'
    }
}
charliecloud {
    enabled = true
}
```

An extensive guide on how to handle containers can be found [here](https://www.nextflow.io/docs/latest/container.html)

- List this new dependency as something new in the `CHANGELOG`
- Create a Pull Request including only these two changes against the `dev` branch of the pipeline you're working on

This way, a review process will be very fast and we can merge the changes into the `dev` branch, updating the Docker Image for that pipeline automatically. After ~30 Minutes, the Docker Image for that pipeline is then updated, and you can open your Pull Request containing your actual pipeline code changes.

## Continuous integration testing

<!-- TODO: Delete / link to another section of the docs -->

To assure that nf-core pipelines don't break after some change is made to the code, we use automated continuous integration (CI) testing. This is done via GitHub actions, which are defined in the `.github/workflows` directory. Parameters and file paths are set in the `conf/test.config` and `conf/test_full.config`. Please see also [here](/docs/contributing/adding_pipelines#add-some-test-data) for how to set-up the test workflow for your pipeline.

## DSL2 and modules

<!-- TODO: Delete / link to another section of the docs -->

Nextflow DSL2 allows for a more modularized design of pipelines and the reuse of components. The nf-core team has developed a set of design patterns on how to best implement DSL2 pipelines, which should be used by all nf-core pipelines in order to assure standardization and the reuse of components. The following is meant to help understand certain design choices and how a nf-core DSL2 pipeline should be build.

### Modules

Each pipeline has a `modules` directory which contains all the module code. A module here depicts a single process which involves - if possible - only a single tool/software. The `modules` directory is furthermore divided into `local`and `nf-core` sub-directories, where local contains the `samplesheet_check.nf`. Modules contained in the `local` directory are specific to the pipeline, whereas `nf-core` modules are installed from the `nf-core/modules` repository. For instance, most pipelines that involve FastQ files will run the FastQC tool for quality control. The module needed for this can be easily reused from the `nf-core/modules` directory using the `nf-core/tools`package.

For more information and a comprehensive guide on the guidelines of how to implement modules in pipelines please refer to the [DSL 2 Modules](https://nf-co.re/docs/contributing/modules) page

### Sample meta information

In nf-core DSL2 pipelines, every channel that contains sample data in some way should also contain a `meta`variable, which must contain the fields `meta.id`, `meta.single_end` and `meta.strandedness`. The `meta` variable can be passed down to processes as a tuple of the channel containing the actual samples, e.g. FastQ files, and the `meta` variable. This meta information can easily be extracted from a samplesheet which specifies the input files. For an example process that reads a samplesheet, creates a `meta` variable and returns it along with the filepaths, have a look at the [rnaseq pipeline](https://github.com/nf-core/rnaseq/blob/master/modules/local/samplesheet_check.nf).
