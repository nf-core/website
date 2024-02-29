---
title: Pipeline Review Guidelines
subtitle: Suggestions for reviewing pipeline pull requests
---

The aim is to have standardised best-practice pipelines.
To ensure this standardisation, we maintain a set of guidelines which all sanger-tol pipelines must adhere to.
These are adopted from [nf-core guidelines](https://nf-co.re/docs/contributing/guidelines) and [review checklist](https://nf-co.re/docs/contributing/pipeline_release_review_guidelines).

Pipeline developers are recommended to create **modular and small pull requests (PRs)** to get the most out of the review process.
Think about that _before_ writing the code and opening the pull-request, as breaking down a PR into multiple ones can be tricky.
As a rule of thumb, a PR should not add more than one sub-workflow, a sub-workflow should not contain more than ten steps. A PR can modify multiple sub-workflows, as long as the changes are related.

The role of the reviewer is to check for adherence to the central principles of nf-core and sanger-tol (reproducibility, excellent reporting, documented, keeping to the template etc.,). Here we provide a general set of suggestions when doing pipeline reviews:

The instructions below are subject to interpretation and specific scenarios. If in doubt, please ask for feedback.

## Do: nf-core principles {#requirements}

All sanger-tol pipelines _must_ follow the following guidelines:

- [Identity and branding](#identity-and-branding): Primary development must on the sanger-tol organisation.
- [Workflow size](#workflow-size): Not too big, not too small.
- [Workflow name](#workflow-name): Names should be lower case and without punctuation.
- [Use the template](#use-the-template): All sanger-tol pipelines must be built using the nf-core template and sanger-tol branding.
- [Software license](#software-license): Pipelines must be open source, released with the MIT license.
- [Bundled documentation](#bundled-documentation): Pipeline documentation must be stored in the repository and viewable on the pipeline website.
- [Docker support](#docker-support): Software must be bundled using Docker and versioned.
- [Continuous integration testing](#continuous-integration-testing): Pipelines must pass CI tests.
- [Semantic versioning](#semantic-versioning): Pipelines must use stable release tags.
- [Standardised parameters](#standardised-parameters): Strive to have standardised usage.
- [Single command](#single-command): Pipelines should run in a single command.
- [Keywords](#keywords): Excellent documentation and GitHub repository keywords.
- [Pass lint tests](#pass-lint-tests): The pipeline must not have any failures in the `nf-core lint` tests.
- [Credits and Acknowledgements](#credits-and-acknowledgements): Pipelines must properly acknowledge prior work.
- [Minimum inputs](#minimum-inputs): Pipelines should be able to run with as little input as possible.
- [Use sanger-tol git branches](#git-branches): Use `main`, `dev` and `TEMPLATE`.

## Do: nf-core recommendations {#recommendations}

All sanger-tol pipelines _should_ follow the following guidelines, if possible / appropriate:

- [Use Bioconda](#use-bioconda): Package software using bioconda and biocontainers.
- [File formats](#file-formats): Use community accepted modern file formats such as `CRAM`.
- [DOIs](#dois): Pipelines should have digital object identifiers (DOIs).
- [Publication credit](#publication-credit): Pipeline publications should acknowledge the sanger-tol community and contributing members.

## Do: Local code and modules {#local}

- Do local scripts in `bin/` have author and license embedded?
  - Local scripts must be licensed with the MIT license, like the pipeline code itself.
- Do all local modules have docker/singularity/conda declarations?
  - Are they ideally in bioconda/biocontainers ?
- Do all local modules conda/container tool declarations have versions? (and _not_ `latest`, `dev` etc.)
- Do all local modules report versions (if applicable)?
  - Simple modules with e.g. single `grep` operations not necessary
  - It would be good to add with more complex operations such as `awk`
- Should any local modules be in nf-core/modules?

## Do: Documentation {#documentation}

- Documentation is only on the pipelines website (not pointed to other places, e.g. not readthedocs )
- Is documentation sufficiently described (`usage.md`, `output.md`, `nextflow_schema.json`)?
  - nextflow_schema.json: check if types are correct and that `default` and `enum` are used where applicable
- Are there any typos in the documentation (`usage.md`, `output.md`, `nextflow_schema.json`)
- Is CHANGELOG sufficiently filled in?
  - Check version system is three-point SemVer e.g. 2.1.0
  - Has the date been updated?
- Check citation formatting consistency in `CIATIONS.md`
- Check that all tools are cited
- Check that (all) pipeline author(s) listed themselves in the manifest and other contributors are added in the README

## Do: Code {#code}

- Check no overly non-template components (no readthedocs, entirely custom logo etc.)
- Check for general code readability
- Check for possible code bugs
- Check for consistency in parameters
  - i.e. `snake_case`
  - All boolean parameters evaluate to `false` by default (e.g. bad: `params.run_step = true`, good: `params.skip_step = false` )
- Check manifest includes DOI (if present) etc.

## Don't have to do {#not-needed}

- Review module code from nf-core/modules
- Comment on scientific content (unless you are familiar with the topic)
- Major code optimisation
  - You _can_ suggest small code optimisations
  - Larger ones you can recommend, but should not necessarily be _required_ for release

## If the guidelines don't fit {#exceptions}

We appreciate that the above guidelines are relatively rigid and may not always fit. If that's the case, please discuss at one of the [pipelines meeting](/public_html/events.php).

We hope that the nf-core best practices, tooling and community are helpful for anyone building Nextflow pipelines, even if they are not a good fit for being listed as official nf-core pipelines.

If a pipeline is found to be violating the standards and guidelines, you should try to address the problems with the pipeline maintainers through discussion. Hopefully the pipeline can then be updated so that it adheres to the guidelines.

All members of the sanger-tol community must adhere to the [sanger-tol code of conduct](/code_of_conduct).
The guidelines and actions within the code of conduct take precedence over the development guidelines described in this page.

## Guidelines

### Identity and branding

Please don't call your pipeline `nf-core/<yourpipeline>`, it _must_ be `sanger-tol/<yourpipeline>`.
Please say that your pipeline _"uses"_ nf-core rather than rather than _"is"_ nf-core.
When you generate a pipeline with `nf-core create`, exclude nf-core branding and select custom prefix `sanger-tol`.

**Development must on the sanger-tol organisation.**

All ToL developers have got write access to the sanger-tol repositories so that all development can happen directly there.

Do _not_ fork sanger-tol repositories.

When new pipelines are added to sanger-tol, please transfer ownership to sanger-tol instead of forking it.

If you have already forked your pipeline to sanger-tol, you can [email GitHub support](https://support.github.com/contact?subject=Reroute%20a%20Fork&tags=rr-forks) and request that they reroute the fork. Alternatively, [contact the IT team](https://github.com/sanger-tol/pipelines-website/issues/new?assignees=muffato%2Cmuffato&labels=connect&projects=&template=contact_us.yaml&title=%5BContact+Us%5D%3A+) and we may be able to help.

**Disable GitHub features for forks**

To encourage contributors to focus on the sanger-tol repository, please disable GitHub issues / wiki / projects on your forked repository. You'll find these options under the GitHub repository settings.

### Workflow size

We aim to have a _"not too big, not too small"_ rule. This is deliberately fuzzy, but as a rule of thumb workflows should contain at least three processes and be simple enough to run that a new user can realistically run the pipeline after spending ten minutes reading the docs.

Most pipelines sizes depends on the scope of the work but please consider disk, RAM and other resources. Larger the usage, fewer parallel work jobs can be managed.

### Workflow name

All sanger-tol pipelines should be lower case and without punctuation. This is to maximise compatibility with other platforms such as Docker Hub, which enforce such rules. We prefer that they are descriptive towards the data or analysis type the pipeline will be using or performing, and should be approved. In documentation, please refer to your pipeline as sanger-tol/pipeline.

### Use the template

All sanger-tol pipelines must be built using the nf-core template with a custom prefix `sanger-tol`.

Workflows should be started using the `nf-core create` command which makes a new git repository and the initial commits and branches. This is to ensure that the sync process can work. See the [sync docs](https://nf-co.re/docs/contributing/sync) for details.

Where possible, workflow authors should do their best to follow nf-core conventions for filenames and code locations.

### Software license

All sanger-tol pipelines must be released with an MIT license. The copyrights belong to Genome Research Ltd. as per Wellcome Sanger Institute policy.

Please try not bundle any third party scripts within the workflow, in case they have a different or incompatible license (for example, in the `bin` directory). If you need such a script, even a simple one, please release it on bioconda instead and reference it like any other software.

### Bundled documentation

All documentation must be bundled with the pipeline code in the main repository, within a directory called `docs`.

Documentation must _only_ be hosted on the GitHub repository, which is automatically synchronised to the pipelines website.
Hosting the documentation at a second location (such as custom readthedocs website, or GitHub pages etc) is not allowed.
This is to ensure that users of sanger-tol pipelines can always intuitively find the documentation for all sanger-tol pipelines in the same way.

Documentation must include at least the following files:

- `README.md`
- `docs/usage.md`
- `docs/output.md`

### Docker support

Pipelines must have all software bundled using [Docker](https://www.docker.com/) - that is, it must be possible to run the pipeline with `-profile docker` and have all software requirements satisfied.

Tools should use docker images from [biocontainers](https://biocontainers.pro) where possible, as using Bioconda / Biocontainers gives support for conda + docker + singularity.

All containers must have specific, stable versions pinned.
These should preferably be named after a software release, but it can also be by commit or some other identifier.
Software versions must be static and stable. Labels such as `latest`, `dev`, `master` and so on are not reproducible over time and so not allowed.

### Continuous integration testing

Pipelines must have automated continuous integration testing, running using GitHub Actions. There must be a small dataset that can be tested on GitHub directly, and a larger one that can be tested on the Sanger farm using Nextflow Tower.

There must be a config `profile` called `test` that should be as comprehensive as possible - that is, it should run as much of the pipeline as possible. It should use as tiny test data set as possible (even if the output that it creates is meaningless).

Then, we configure the integration with Nextflow Tower to allow testing the larger dataset (`test_full`) on the Sanger LSF farm. To set up that up, first add the profile `cleanup { cleanup = true }` to your `nextflow.config` (right at the beginning of the `profiles` section). This is to control the amount of space taken on Lustre. Then, copy the two files [`sanger_test.yml`](https://github.com/sanger-tol/insdcdownload/blob/dev/.github/workflows/sanger_test.yml) and [`sanger_test_full.yml`](https://github.com/sanger-tol/insdcdownload/blob/dev/.github/workflows/sanger_test_full.yml) to your `.github/workflows/`. Ask @muffato to enable the Tower integration for your repository.

### Semantic versioning

Pipelines must be released with stable release tags.
Releases must use GitHub releases and [keep a detailed changelog](https://keepachangelog.com/en/1.0.0/) file.

Release version tags must be numerical only (no `v` prefix) and should follow [semantic versioning](https://semver.org/) rules: `[major].[minor].[patch]`

For example, starting with with a release version `1.4.3`, bumping the version to:

- `1.4.4` would be a patch release for minor things such as fixing bugs.
- `1.5.0` would be a minor release, for example adding some new features, but still being backwards compatible.
- `2.0.0` would correspond to the _major_ release where results would no longer be backwards compatible.

### Standardised parameters

Where possible pipelines should use the same command line option names as other pipelines for comparable options.
For example, `--input` and `--fasta`.

In addition to the names of parameters, they should ideally work in a similar way.
For example, `--input` typically takes a `.csv` sample sheet file (but not always, where not appropriate).

> nf-core are planning to build a tool that lists every parameter used by every pipeline, so that you can check for existing parameters with similar names
> You can track the progress of this feature request here: [nf-core/website#1251](https://github.com/nf-core/website/issues/1251)
> Once comple, we will try to implement this for sanger-tol pipelines as well.

### Single command

Every sanger-tol pipeline repository must contain a single pipeline.
That is, there should be a `main.nf` file that is the single way to launch a pipeline.

- It is ok to have multiple 'tracks' within the pipeline, selectable with configuration options.
- It is ok to have workflows that use the output of _another_ nf-core pipeline as input

It should be possible to run all parts of the workflow using `nextflow run sanger-tol/<pipeline>`, without any specific `.nf` filename.

### Keywords

Pipelines should have excellent documentation.

Repositories should have [GitHub Topics](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/classifying-your-repository-with-topics) set on the sanger-tol repository.
These are then shown on the pipelines website and used for categorisation and searching. They are important for workflow visibility and findability.

Topics can be related to the workflow, tools, data or any related aspect.

### Pass lint tests

In order to automate and standardise the nf-core best practices, there is a code linting tool. These tests are run by the [nf-core/tools](https://github.com/nf-core/tools) package. The `nf-core lint` command must be run by continuous integration tests on GitHub Actions and must pass for each pull request and before release.

You can see the list of tests and how to pass them on the [error codes page](https://nf-co.re/tools).

In some exceptional circumstances, it is ok to ignore certain tests using `nf-core.yml`. If that's the case, please discuss first at one of the [pipelines meeting](/public_html/events.php).

### Credits and Acknowledgements

Please acknowledge all developers, reviewers and anyone who has contributed to the finished product (even through verbal discussions and suggestions). Everyone should be credited in alphabetical order by last name. This is independent of authorship in any associated manuscript.

Where previous work from other pipelines / projects is used within a pipeline, the original author(s) must be properly acknowledged. Some examples on how you could do that to make sure they feel valued:

- Send them a message via Slack and let them know that you use their work and had to change something to fit your own purpose. If in doubt, check with them to see how they would like to be acknowledged.
- Check the license of their code and/or graphics components and make sure you obey the rules that this license imposes (e.g. `CC-BY` means you have to attribute the original creator).
- If you use portions of pipeline code, even if its just tiny pieces:
  - Link to the original repository and/or authors.
  - Leave existing credits and acknowledgement sections intact - there may be more than just a single author involved.
- If you find bugs / issues, report and fix them upstream in the main project.

If in doubt about what to do, ask on Slack or discuss at the fortnightly pipeline meetings.

### Minimum inputs

Pipelines can accept as many input files as you like, but it should be possible to run with as few as possible.

For example, pipelines should auto-generate missing reference files, where possible.
So given a reference genome Fasta file a pipeline would build the reference index files. The pipeline should also be able to optionally accept the reference index files in this case, if available.

### Git branches

The latest stable release should be on the main `main` branch.
No additional changes should be pushed to master after each release.

The main development code should be kept in a branch called `dev`.
The sanger-tol `dev` branch should be set as the default branch up until the first release.

For minor bugfixes a `patch` branch may be used and merged directly into `main`, leaving `dev` for continued development work.

The `TEMPLATE` branch should only contain vanilla nf-core template code.
It is used for automated synchronisation of template updates.

### Use Bioconda

All pipeline software should be packaged using [bioconda](https://bioconda.github.io/).
Bioconda packages are automatically available as Docker and Singularity images through [biocontainers](https://biocontainers.pro/).

### File formats

Pipelines should work with best practice modern file formats, as accepted by the community.

Where possible, genomics pipelines should generate `CRAM` alignment files by default, but have a `--bam` option to generate `BAM` outputs if required by the user.

### DOIs

Pipelines must have digital object identifiers (DOIs) for easy referencing in literature.

Typically each release should have a DOI generated by [Zenodo](https://zenodo.org/).
This can be automated through linkage with the GitHub repository.

### Cloud compatible

Pipelines should have explicit support for running in cloud environments. The pipelines created with [nf-core template](https://nf-co.re/tools/#creating-a-new-pipeline) comes with all required code to support this setup.

### Publication credit

**sanger-tol**

Pipeline publications should acknowledge all developers, reviewers and anyone who has contributed to the finished product (even through verbal discussions and suggestions). This can be done by either listing people by name or as a collective.

**nf-core**

**At a minimum**, the nf-core community should be thanked in the acknowledgment section.

> We would like to thank the nf-core community for developing the nf-core infrastructure and resources for Nextflow pipelines. A full list of nf-core community members is available at https://nf-co.re/community.

**Optionally**, the nf-core community can be included as a consortium co-author of the publication ([example](https://doi.org/10.3390/ijms232314512)). This route is a good idea when a pipeline has used extensive existing nf-core pipeline infrastructure (e.g. modules) that were written by other members of the community not directly involved in the pipeline itself.

If any members of the nf-core community have provided significant input to the creation of the pipeline, please consider adding them as coauthors on the paper directly.

If in doubt, contact the core team for guidance.
