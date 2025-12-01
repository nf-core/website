---
title: First release
subtitle: Follow this walkthrough to add a new pipeline to nf-core.
weight: 40
type: "tutorial"
markdownPlugin: "checklist"
---

# Making the first release

## Reset the default branch

When the code is stable and ready for a release, set the default branch to `master`.

## Bump the version

Use `nf-core pipelines bump-version` to increment the version number on the `dev` branch to remove
the `dev` suffix. The first release of a pipeline should normally be version `1.0.0` (we discourage
pre-releases).

When developing the pipeline, the version numbers should be numeric with `dev` at the end,
e.g. `0.0.0dev`. When making a release, version numbers should all be numeric. Pipeline
release numbers MUST use [Semantic Versioning](https://semver.org/).

:::tip
Instead of manually changing the version number, use the `nf-core pipelines bump-version` command to increment the version numbers. The version number
exists in many places in the codebase, and this tool consistently updates all of them.
:::

## Lint the pipeline

Use `nf-core pipelines lint --release`. This will check that your code conforms to nf-core expectations.

## Core pipeline review

Ok - now the tough bit - does your workflow stand up to the scrutiny of the nf-core team?!
Not to worry, we're a friendly bunch.

To get your pipeline reviewed for its initial release, do the following:

1. Make a pull-request from the `dev` branch to `master` on the nf-core fork. This is a
   special case and all of the tests should pass.

1. Once the tests are passing, request two reviews, with one from the core or maintainers team. You can use the
   [#release-review-trading](https://nfcore.slack.com/archives/C08K66XCZSL) slack channel for this.

What happens next depends on the state of your master branch:

- If you have developed in such a way that your master branch is clean, .i.e. doesn't have
  any commits since the initial one, the PR will represent all changes
  associated with the proposed release, and the core team will use it for review and
  feedback.
- If your master branch already contains changes associated with the release, the core
  team may merge your PR and create a pseudo-PR against the first commit in the
  pipeline. This gives the PR review interface showing all code that you've written.

In either case the two reviewers will go through everything and request changes that we think are
necessary until you're good to go. Common things that are flagged at this point are:

- A clear, short but descriptive readme
- Good documentation, especially describing the output files and all parameters
- Pipeline code

We typically tend to have two reviewers for most of the crucial code changes, e.g. adding
new major features to an existing pipeline or making an entirely new pipeline release. You
can also ping people from the nf-core core or maintainers team to review your pipeline
code by `@`ing them.

:::tip
While you're waiting, try reviewing another release pull request from another community member!

Don't feel shy, it's a good way to start and get a better understanding reviews work, and the more eyes the better!

See the documentation pinned on the [#release-review-trading](https://nfcore.slack.com/archives/C08K66XCZSL) for instructions.
:::

## Making the release

Once any requested changes have been made and the associated PR approved, you can go ahead
with releasing the pipeline. Put in a basic changelog entry describing the general
functionality at release. When you're ready, follow the instructions in the nf-core
[release checklist](/docs/checklists/pipeline_release). We recommend that you also explicitly
tag contributors with their GitHub handles, so each release on GitHub will display their icons.

The nf-core website and helper tools will automatically detect the new release and be updated accordingly.

That's it, you're finished! Congratulations!

:::note{title=Publications}
If you wish to make a publication based on the pipeline, make sure you [give credit to nf-core](/docs/guidelines/pipelines/recommendations/publication_credit).
:::

Once you've finished your release - it's time to tell the world about it!
Post on the [#bytesize_suggestion](https://nfcore.slack.com/archives/C081F8J2X8R) slack channel of your release, to begin arranging a 15 minute introductory bytesize talk about your shiny new pipeline!

### Developing nf-core pipelines privately

Generally we would like everyone to develop nf-core pipelines in the open, but we understand that this is not always possible. If you are developing a pipeline privately, you can still release it through nf-core, but you will need to be aware of a few things.

Pipelines developed in private won't be given any special treatment and will be reviewed in the same way as other new pipelines. As such, we recommend that you follow the same steps as above, but with a private repository. Once the pipeline is ready for release, you can make it public and follow the steps above. If you can not make it public until after release, you can request community members to review the pipeline in the nf-core Slack channel `#request-review` before making it public.

Additionally, private pipelines will not be given preference over another pipeline that has been developed in the open. If you are developing a pipeline privately, please be aware that another pipeline may be released before yours. In this case, we may not be able to accept your pipeline into nf-core and you may need to release it independently or merge it with the pipeline that has already been released.

If you are considering developing a pipeline in private that you hope to release it through nf-core, please contact the nf-core core-team on Slack so we can discuss a pre-release strategy.

# Showing your pipeline in the nf-core website

Make sure that your pipeline does not appear in the [list of repositories](https://github.com/nf-core/website/blob/75a92b76129b0e2c3f644bf26482142ba9e41f6f/sites/main-site/src/config/ignored_repos.yaml) ignored form the website build.

If your pipeline appear there, open a PR removing it, and ask a core-team member to approve the PR.

# Subsequent releases

Once you've made your first release you can continue to work on your fork and make pull-requests
against the `dev` branch on the nf-core repository. Now that we have a stable `master` branch,
there should be reviews of each PR against `dev` before merging.

When ready to make new releases, make sure that the version number is increased with
`nf-core pipelines bump-version` and create a pull-request against `master`. If tests pass, it
can be merged and a new release made.

The `master` branch should always have only the commit from the latest release. This is important
because the commit ID is used to reference whether the pipeline is up to date or not.
