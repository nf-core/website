---
title: First release
subtitle: Follow this walkthrough to add a new pipeline to nf-core.
weight: 40
type: 'tutorial'
markdownPlugin: 'checklist'
---

# Making the first release

## Reset the default branch

When the code is stable and the pipeline is ready to be released, set the default branch to `master`.

## Bump the version

When developing the pipeline, the version numbers should follow numeric semantic versioning with `dev` at the end,
e.g. `0.0.0dev`. When preparing for a release, use `nf-core bump-version` to increment the version number on the `dev` branch and remove
the `dev` suffix.

The pipeline version number MUST use [semantic versioning](https://semver.org/).

:::tip
Instead of manually changing the version number, use the `nf-core bump-version` command to increment the version numbers. The version number
exists in many places in the codebase, and this tool consistently updates all of them.
:::

## Lint the pipeline

Use the `nf-core lint --release` command to check that the code conforms to nf-core release expectations.

## Core pipeline review

Ok - now the tough bit - does the pipeline stand up to the scrutiny of the nf-core team?!
Not to worry, we're a friendly bunch.

To get the pipeline reviewed for its initial release, do the following:

1. Make a pull-request from the `dev` branch to `master` on the nf-core fork. This is a
   special case and all of the tests should pass.

1. When the tests are passing, request a review from the core team. You can use the
   [#request-review](https://app.slack.com/client/TE6CZUZPH/CQY2U5QU9) slack channel for this.

What happens next depends on the state of the master branch:

- If you have developed in such a way that the master branch is clean, .i.e., it doesn't have
  any commits since its creation, the pull request will represent all changes
  associated with the proposed release, and the core team will use it for review and
  feedback.
- If the master branch already contains changes associated with the release, the core
  team may merge the PR and create a pseudo-PR against the first commit in the
  pipeline. This enables the reviewer to see all of the new code that has been written.

In either case, the code will be reviewed and some changes may be requested. Common things that are flagged at this point include:

- A clear, short but descriptive readme
- Good documentation, especially describing the output files and all parameters
- Pipeline code

Two reviewers are typically required for most code changes, e.g., adding
new major features to an existing pipeline or making an entirely new pipeline release. You
can ping people from the nf-core core or maintainers team to review the pipeline
code by `@`ing them.

## Making the release

The pipeline can be released once any requested changes have been made and the associated PR approved. Add a short changelog entry describing the general
functionality at release.

When you're ready, follow the instructions in the nf-core
[release checklist](/docs/checklists/pipeline_release). We recommend that contributors are explicitly
tagged with their GitHub handles, so each release on GitHub will display their avatars.

The nf-core website and helper tools will automatically detect the new release and will be updated accordingly.

That's it! You're finished! Congratulations!

:::note{title=Publications}
Make sure you [give credit to nf-core](/docs/guidelines/pipelines/recommendations/publication_credit) if you are publishing the pipeline.
:::

# Subsequent releases

After a first release you can continue to work on your fork and make pull-requests
against the `dev` branch on the nf-core repository. Now that there is a stable `master` branch,
there should be reviews of each pull request against `dev` before merging.

When ready to make new releases, make sure that the version number is increased with
`nf-core bump-version` and create a pull-request against `master`. If tests pass, it
can be merged and a new release made.

The `master` branch should always have only the commit from the latest release. This is important
because the commit ID is used to reference whether the pipeline is up to date or not.
