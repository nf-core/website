---
title: Release checklist
subtitle: A step-by-step guide for releasing a nf-core pipeline
markdownPlugin: checklist
weight: 10
---

## Before you release

:::warning
If this is your first release, make sure to follow the [first release tutorial](/docs/tutorials/adding_a_pipeline/first_release) for extra review requirements!
:::

- [ ] Check the pipeline follows [nf-core guidelines](/docs/guidelines/pipelines/overview).
- [ ] All continuous-integration tests on the `dev` branch should be passing
  - [ ] Check the output of `nf-core pipelines lint` for warnings and address all that you can
  - [ ] Update any software dependencies that are out of date.
        The linting will warn about available updates via automated API calls to (bio-)conda
- [ ] Check that there are no outstanding issues that need to be addressed, especially bug reports.
- [ ] Finalize the description of the pipeline in the GitHub repository and ensure that you remove any "under development" labels in the description. This wording will be used when a new pipeline gets announced on our socials ([Bluesky](https://bsky.app/profile/nf-co.re) and [Mastodon](https://mstdn.science/@nf_core)).
- [ ] If there is a [release milestone](https://help.github.com/en/github/managing-your-work-on-github/about-milestones), have a look and see if all issues are closed, or can be resolved
  - [ ] It's fine to decide that some things should be postponed until the next release - just move them on to the next milestone
  - [ ] Check that the full-size tests have successfully completed within the `dev` branch

## Steps to release

- [ ] On your fork, bump the version number of the `dev` branch to a release version
  - [ ] For example, `1.0.0dev` becomes `1.0.0`
  - [ ] Use the `nf-core pipelines bump-version` command to make the changes, eg: navigate to the pipeline directory and run `nf-core pipelines bump-version 1.0.0`
  - [ ] Please make sure to use strictly numeric release numbers (must be `1.0.0` NOT `v1.0.1rc`)
  - [ ] Use [Semantic Versioning](https://semver.org/) - see the [nf-core semantic versioning guidelines](/docs/guidelines/pipelines/requirements/semantic_versioning.md) to decide if you need a patch, minor, or major release.
  - [ ] Make sure to update the version in any pipeline diagrams and other figures
- [ ] Run `nf-core pipelines lint --release` and check that there are no test failures for release.
- [ ] Check that `CHANGELOG.md` includes everything that has been added/fixed in this release, update the version number above the changes, and optionally add a human-readable release name (e.g. using a [code name generator](http://www.codenamegenerator.com/))
  - [ ] We recommend you also add the GitHub handle of the main contributors of each CHANGELOG entry (author, and significant reviewers etc.). This will mean each release on GitHub will display each contributors icons for extra visibility and recognition.
- [ ] [Open a Pull Request (PR)](https://help.github.com/en/articles/creating-a-pull-request) with these changes from your fork to the `dev` branch on the nf-core repository.
- [ ] Once merged, open another PR from the nf-core `dev` branch to the nf-core `master`
  - [ ] Make sure that all of the CI tests are passing - this is a special case PR and the tests are different
- [ ] Request PR reviews from at least two people on [#release-review-trading](https://nfcore.slack.com/archives/C08K66XCZSL)
  - [ ] If a first release, one review must come from one core or maintainers team member - feel free to ping the two teams on slack under your post on [#release-review-trading](https://nfcore.slack.com/archives/C08K66XCZSL).
  - [ ] For any other release, one review can come from a member of the same development team, but we strongly recommend the second to be an another external/neutral person from the community
  - [ ] While you're waiting for reviews on your PR, review another community member's release PR (see pinned message in the channel above)
- [ ] Once your PR approved by two reviewers, merge your PR into `master`/`main`
  - [ ] (**First release only**) And finally delete any label of types: "under development", "under construction" or variants of these on the repository
- [ ] (**First release only**) Ask on the [#request-core](https://nfcore.slack.com/archives/C09H6NYHR9T) channel to activate the Zenodo functionality for this repository, which will be used to generate a DOI.
- [ ] Go to GitHub and [create a new release for your pipeline](https://help.github.com/en/articles/creating-releases)

  :::note
  Use _exactly_ the same version as in the code (e.g. `1.0.0`) - **do not prefix with v** (e.g. not `v1.0.0`).
  :::
  - [ ] Optional: Also include your [nice code name](http://www.codenamegenerator.com/) in your pipeline release title (see above with `CHANGELOG.md`).
        For example releases in nf-core/rnaseq follow the pattern:
    - Prefix = Metal
    - Dictionary = Animals
    - Suffix = Don't use a suffix

- [ ] Celebrate! But not too much - you still have a few things left to do...

### Automated events

A number of events are automatically triggered after the pipeline is released:

- [ ] The [nf-core website](https://nf-co.re/pipelines) will be updated automatically with the release information.
- [ ] Our social media accounts ([Bluesky](https://bsky.app/profile/nf-co.re) and [Mastodon](https://mstdn.science/@nf_core)) will send out an automated post about the pipeline release within minutes.
- [ ] A [Zenodo DOI](https://zenodo.org/) is automatically generated that provides a persistent means with which to cite the pipeline.

## After release

The last step is to bump up the pipeline version number in the development branch:

- [ ] Make sure the `dev` branch on your fork is up to date
- [ ] On a new branch, bump to a new `dev` version
  - [ ] For example, `1.0.0` becomes `1.1.0dev`
  - [ ] Use the `nf-core pipelines bump-version` command to make the changes, eg: navigate to the pipeline directory and run `nf-core pipelines bump-version 1.1.0dev`
- [ ] Update the `CHANGELOG.md` to include a new section for this new version
- [ ] If your pipeline is using pipeline-level nf-test tests, update snapshots: `nf-test test tests/ --profile=+<docker/singularity/conda etc.>`
- [ ] [Open a Pull Request (PR)](https://help.github.com/en/articles/creating-a-pull-request) with these changes from your fork to the `dev` branch on the nf-core repository.
- [ ] (**First release only**) After the first release of the pipeline you will need to add the DOI manually into the main `README.md` for the pipeline:
  - [ ] Search for your pipeline on Zenodo and find the DOI that allows you to _"Cite all versions"_ of the pipeline.
  - [ ] Ask a core member to copy the DOI information you added to dev via the PR above to the master branch. The core member will uncomment the Zenodo-related `TODO` statement in the `Citation` section of the main `README.md` and add the DOI, as well as as updating the badge for the Zenodo DOI at the top of the main `README.md` e.g. [nf-core/atacseq](https://github.com/nf-core/atacseq/blob/fa1e3f8993cd20e249b9df09d29c5498eff311d2/README.md).
- [ ] (**First release only**) Ask a core member to change default branch from `dev` to `master`.
- [ ] Check that AWS megatest run succesfully and that the results are made available on the website e.g.[nf-core/atacseq](https://nf-co.re/rnaseq/3.18.0/results).
- [ ] (publication only) If a publication of the pipeline is being prepared, recommended [nf-core guidelines](/docs/guidelines/pipelines/recommendations/publication_credit) are followed.
- [ ] Post on the [#bytesize_suggestion](https://nfcore.slack.com/archives/C081F8J2X8R) slack channel of your release, to begin arranging a 15m introductory bytesize talk about your shiny new pipeline!
