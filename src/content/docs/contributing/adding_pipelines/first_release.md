---
title: Adding a new pipeline
subtitle: Follow this walkthrough to add a new pipeline to nf-core.
---

# Making the first release

## Reset the default branch

When the code is stable and ready for a release, set the `master` branch to be the default
branch again.

## Bump the version

At this point you should bump the version numbers on `dev`.

When developing the pipeline, the version numbers should be numeric with `dev` at the end.
Use the `nf-core bump-version` command to do this - there are quite a few locations in the
code that need updating and this ensures that it happens in the correct places.

When making a release, version numbers should all be numeric. Use `nf-core lint --release`
when ready - this will check that everything looks correct. Pipeline release numbers MUST
use [Semantic Versioning](https://semver.org/).

## Core pipeline review

Ok - now the tough bit - does your workflow stand up to the scrutiny of the nf-core team?!
Not to worry, we're a friendly bunch, just let us know about the new pipeline, when you're
ready, following the process below.

Make a pull-request from the `dev` branch to `master` on the nf-core fork. This is a
special case and the tests should pass, and once they do you can request a review from the
core team.

What happens next depends on the state of your master branch:

- If you have developed in such a way that your master branch is clean, .i.e. doesn't have
  any commits since the initial one, the PR created above will represent all changes
  associated with the proposed release, and the core team will use it for review and
  feedback.
- If your master branch already contains changes associated with the release, the core
  team may merge your PR and create a pseudo-PR against the first commit in the
  pipeline. This gives the PR review interface showing all code that you've written.

In either case we will go through everything and request changes that we think are
necessary until you're good to go. Common things that are flagged at this point are:

- A clear, short but descriptive readme
- Good documentation, especially describing the output files and all parameters
- Pipeline code

We typically tend to have two reviewers for most of the crucial code changes, e.g. adding
new major features to an existing pipeline or making an entirely new pipeline release. You
can also ping people from the nf-core core or maintainers team to review your pipeline
code by `@`ing them.

## Making the release

Once any requested changes have been made and the associated PR approved, you can go ahead
with releasing the pipeline. Put in a basic changelog entry describing the general
functionality at release. When you're ready, follow the instructions in the nf-core
[release checklist](/docs/contributing/release_checklist). We recommend you also explicitly tag contributors with their GitHub handles, so each release on GitHub will display their icons.

The nf-core website and helper tools will automatically detect the new release and be updated accordingly.

That's it, you're finished! Congratulations!

:::note{title=Publications}
If you wish to make a publication based on the pipeline, make sure you [give credit to nf-core](/docs/contributing/guidelines/recommendations/publication_credit).
:::

# Subsequent releases

Once you've made your first release you can continue to work on your fork and make pull-requests
against the `dev` branch on the nf-core repository. Now that we have a stable `master` branch,
there should be reviews of each PR against `dev` before merging.

When ready to make new releases, make sure that the version number is increased and create a
pull-request against `master`. If tests pass, it can be merged and a new release made.

The `master` branch should always have only the commit from the latest release. This is important
because the commit ID is used to reference whether the pipeline is up to date or not.
