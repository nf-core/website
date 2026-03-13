---
title: Release procedure
subtitle: A step-by-step guide for releasing a nf-core pipeline
markdownPlugin: checklist
---

## Before you release

:::warning
If this is your first release, follow the [first release guide](/docs/contributing/contribute-new-pipelines#make-the-first-release) for extra review requirements.
:::

- [ ] Check on the issues tab of GitHub repository that there are no outstanding bug reports that should be resolved.
- [ ] If you have a [release milestone](https://help.github.com/en/github/managing-your-work-on-github/about-milestones), check all planned functionality and improvements are addressed.
  - [ ] Move any postposed issues to the next milestone
- [ ] Check the pipeline follows [nf-core guidelines](/docs/specifications/pipelines/overview) and [nf-core pipeline release review guidelines](/docs/contributing/reviewing-pull-requests/review_checklists/pipeline-release).
- [ ] Check that all continuous-integration tests on the `dev` branch pass
  - [ ] Check the output of linting for warnings and address all that you can resolve.

    ```bash
    nf-core pipelines lint
    ```

  - [ ] If possible, update any software dependencies that are reported to be out of date.

- [ ] Check that the full-size tests have successfully completed within the `dev` branch.
- [ ] Review the `CHANGELOG.md` to ensure it includes everything that has been added/fixed in this release.
  - Add the GitHub handle of the main contributors of each CHANGELOG entry (author, significant reviewers, etc.).
  - This will mean each release on GitHub will display each contributors icons for extra visibility and recognition.
- [ ] (**First release**) Finalise the description of the pipeline in the GitHub repository.
- [ ] (**First release**) Remove any "under development" labels in the description, README and so on.
  - Ensure wording will fit in a ([Bluesky](https://bsky.app/profile/nf-co.re) and/or [Mastodon](https://mstdn.science/@nf_core)) post for the automated announcements.

## Steps to release

### Prepare release PR

- [ ] Make a branch on the repository or fork to bump the version number of the `dev` branch to a release version.
  - For example, `1.0.0dev` becomes `1.0.0`.
  - Use nf-core/tools make the changes:

    ```bash
    nf-core pipelines bump-version <1.0.0>
    ```

  - Use strictly numeric release numbers that follow [Semantic Versioning](https://semver.org/) (for example, must be `1.0.0` NOT `v1.0.1rc`).
  - Review the [nf-core semantic versioning guidelines](/docs/specifications/pipelines/requirements/semantic_versioning) to decide if you need a patch, minor, or major release (for example `1.0.0` not `1.0`).

- [ ] Update the version in any other files not covered by the nf-core/tools automation (for example, pipeline diagrams, tutorials, and other figures).
- [ ] Run the nf-core release linting and check that there are no test failures for release.

  ```bash
  nf-core pipelines lint --release
  ```

- [ ] Update the version number and date in the `CHANGELOG.md` entry and optionally add a human-readable release name.
  - For the human-readable release name, you can use any scheme you would like, for example nf-core/funcscan uses national dishes, nf-core/airrflow uses magic spells from Harry Potter, and nf-core/eager uses the C14th [Swabian League of Cities](https://en.wikipedia.org/wiki/Swabian_League_of_Cities)
  - If you don't know what schema to use, you can use a [code name generator](http://www.codenamegenerator.com/). For example releases in nf-core/rnaseq use: A Metal (Prefix), and Animals (Dictionary), and no suffix.
- [ ] [Open a Pull Request (PR)](https://help.github.com/en/articles/creating-a-pull-request) with these changes from your fork to the `dev` branch on the nf-core repository.
- [ ] Request a review on the version bump PR and merge once approved.

### Open release PR

Once all checks pass and pipeline version updated you are now ready to get your release review:

- [ ] Open the release PR from the nf-core `dev` branch to the nf-core `main`
  - [ ] Make sure that all of the CI tests pass. This is a special case PR and the tests are different
- [ ] Request PR reviews from at least two people on [#release-review-trading](https://nfcore.slack.com/archives/C08K66XCZSL)
  - (**First release only**) Request one review from one core or maintainers team member - feel free to ping the two teams on slack under your post on [#release-review-trading](https://nfcore.slack.com/archives/C08K66XCZSL).
  - One review can come from a member of the same development team, but we strongly recommend the second to be an another external/neutral person from the community.
- [ ] While you're waiting for reviews on your PR, review another community member's release PR requested in [#release-review-trading](https://nfcore.slack.com/archives/C08K66XCZSL).
- [ ] (**First release only**) Ask a core team member on [#release-review-trading](https://nfcore.slack.com/channels/release-review-trading) to remove the **"Block releases"** GitHub repository ruleset.
- [ ] Once your PR is approved by two reviewers, update `CHANGELOG.md` for any further changes and the release date, and merge your PR into `main`.
- [ ] (**First release only**) Delete any label of types: "under development", "under construction" or variants of these on the GitHub repository itself.
- [ ] (**First release only**) Ask on the [#request-core](https://nfcore.slack.com/archives/C09H6NYHR9T) channel to activate the Zenodo functionality to generate DOIs for this repository .
- [ ] (**First release only**) Ask on the [#request-core](https://nfcore.slack.com/archives/C09H6NYHR9T) channel to be added to the nf-core Seqera Platform AWS Megatest workspace.

### Make release

You can now formalise the release:

- [ ] Go to GitHub and [create a new release for your pipeline](https://help.github.com/en/articles/creating-releases).
  - We recommend that you copy your `CHANGELOG.md` entry for the title and release description.

  :::note
  Use _exactly_ the same version as in the code (for example: `1.0.0`). **Do not prefix with v** (for example: not `v1.0.0`).
  :::

- [ ] 🎉 Celebrate! But not too much - you still have a few things left to do...

## Post-release tasks

### Automated events

A number of events are automatically triggered after the pipeline is released:

- [ ] The [nf-core website](https://nf-co.re/pipelines) will be updated automatically with the release information.
- [ ] The full-size AWS megatests for your pipeline's `test_full` should be executed.
- [ ] Our social media accounts ([Bluesky](https://bsky.app/profile/nf-co.re) and [Mastodon](https://mstdn.science/@nf_core)) will send out an automated post about the pipeline release within minutes.
- [ ] A [Zenodo DOI](https://zenodo.org/) is automatically generated that provides a persistent means with which to cite the pipeline.

### After release

Once all the automated steps are completed you should:

- [ ] Check that AWS megatest run succesfully on Seqera Platform and that the results are made available on the website.
  - For example for [nf-core/atacseq](https://nf-co.re/rnaseq/3.18.0/results).
- [ ] (**First release only**) Ask on the [#request-core](https://nfcore.slack.com/archives/C09H6NYHR9T) channel to add on the `main` branch the new DOI to the `nextflow.config` manifest, the badge to top of the the `README.md`, and update the `Citation` section of the `README.md`.
- [ ] (**First release only**) Ask a core member to change default branch from `dev` to `main`.

You should also do some more manual advertising about your pipeline:

- [ ] (**First release only**) Post on the [#bytesize_suggestion](https://nfcore.slack.com/archives/C081F8J2X8R) Slack channel of your release to arrange a 15 minute introductory bytesize talk about your new pipeline.
- [ ] Announce the new release to the members of your dedicated pipeline channel on Slack, as well as any other of your own social media accounts.
- [ ] (**Publication only**) If a publication of the pipeline is being prepared, review the nf-core guidelines on [publication credit](/docs/guidelines/pipelines/recommendations/publication_credit).

Finally, the last step is to bump up the pipeline version number in the development branch:

- [ ] Make sure the `dev` branch on your fork is up to date with `main`.
- [ ] On a new branch, bump to a new `dev` version as in [Prepare release PR](#prepare-release-pr).
  - [ ] For example, `1.0.0` becomes `1.1.0dev`.

    ```bash
    nf-core pipelines bump-version <1.1.0>dev
    ```

- [ ] Update the `CHANGELOG.md` to include a new section for this new version.
- [ ] Update the version in any other files not covered by the nf-core/tools automation (for example, pipeline diagrams, tutorials, and other figures).
- [ ] Update snapshots to include the new version and any other changes:

  ```bash
  nf-test test tests/ --profile=+<docker/singularity/conda-etc.>
  ```

- [ ] [Open a Pull Request (PR)](https://help.github.com/en/articles/creating-a-pull-request) with these changes from your fork to the `dev` branch on the nf-core repository.
