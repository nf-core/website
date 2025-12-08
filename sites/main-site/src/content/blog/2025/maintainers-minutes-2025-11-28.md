---
title: "Maintainers Minutes: November 2025"
subtitle: "Keeping you informed of the latest maintainers discussions"
pubDate: 2025-12-08T10:00:00+01:00
headerImage: "/assets/images/blog/maintainers-minutes-2024-07-26/maintainers-wide.png"
headerImageAlt: "Cartoon yellow rubber duck with nf-core logo badge on its body with the nf-core logo."
embedHeaderImage: false
authors:
  - "mashehu"
label:
  - "maintainers"
---

The 'Maintainers Minutes' aims to give further insight into the workings of the [nf-core maintainers team](/governance#maintainers) by providing brief summaries of the monthly team meetings.

## Overview

This month's maintainers meeting covered several topics around testing infrastructure, release processes, and community contributions.

## nf-test improvements

We now have four nf-core core members with write access to the nf-test repository: [Edmund](https://github.com/edmundmiller), [Sateesh](https://github.com/sateeshperi), [Nicolas](https://github.com/nvnieuwk), and [Maxime](https://github.com/maxulysse).
This increased access should help us better align nf-test development with our nf-core's needs and help push nf-test development forward.

We discussed priority improvements we'd like to see in nf-test, including:

- [The ability to exclude tags](https://github.com/askimed/nf-test/pull/283)
- [Topic channels support](https://github.com/askimed/nf-test/issues/258)
- [Compatibility with strict syntax](https://github.com/askimed/nf-test/issues/326)
- [Local storing of downloaded files/containers/conda environments](https://github.com/askimed/nf-test/issues/231)

These enhancements will make testing more efficient across our pipelines and will allow us to remove some boilerplate code.

## Improving the release process

A significant portion of the meeting focused on making our release process more sustainable and less burdensome for maintainers.

### Trunk-based development

We discussed the proposal for [trunk-based development](https://github.com/nf-core/proposals/issues/49).
While this approach could streamline our workflow, it raises concerns about traceability. Since most beginners run pipelines without the `-r TAG` flag, they would automatically pull the most recent version, making it harder to track which version users are actually running.

### Review bottlenecks

Getting code reviews remains challenging, even for small PRs.
Maintainers often need to personally reach out to get reviews, and the review trade can be quite uneven.

Many contributors also don't feel confident judging scientific accuracy.
We emphasized that maintainers and the core team should focus primarily on the Nextflow code quality and nf-core guidelines compliance, not necessarily the underlying science, which the community around the pipeline should be responsible for.

### Automating guideline checks

For the `dev` → `main`/`master` release process, we discussed focusing on overall guideline compliance rather than re-reviewing individual PRs. Some ideas to reduce this burden include:

- Using AI to automate guideline checking
- Improving our linting tools
- Creating automated linting reports that get pushed to repository issues or PRs

The goal is to catch and fix linting issues earlier in the development process, rather than discovering them during the release.

## Infrastructure updates

[@mashehu](https://github.com/mashehu) has implemented new CI for the nf-core/test-datasets repository that now warns when files are being deleted.
This should help prevent accidental removal of test data.

## Topic channels

Exciting news - the first batch of modules with topic channels for version reporting are now merged and ready to use!
See the [blog post](https://nf-co.re/blog/2025/version_topics) and [migration guide](https://nf-co.re/docs/tutorials/migrate_to_topics/update_modules) for more information.

## The end

We will continue to work on these topics in the coming months, particularly around making the release process more sustainable and less dependent on individual maintainers.

As always, we are looking for community involvement, so we encourage anyone to join the discussion in relevant PRs and Slack channels and threads!

\- :heart: from your #maintainers team!
