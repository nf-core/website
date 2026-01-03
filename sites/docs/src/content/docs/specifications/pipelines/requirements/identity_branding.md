---
title: Identity and branding
subtitle: Primary development must on the nf-core organisation.
menu:
  main:
    weight: 5
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

All nf-core pipelines are [community owned](/docs/guidelines/pipelines/requirements/community_owned).
nf-core pipelines MUST be owned by nf-core, and nf-core alone.

:::info{title="Rationale" collapse}
At its inception, the first step was to remove institutional branding from existing workflows and migrate them to the nf-core organisation and nf-core branding.
This was done to remove institutional and organisational ownership of pipelines and code, to encourage collaboration and openness.
:::

Pipelines MAY be forked on GitHub for personal development work, but the nf-core repository MUST be the primary source for all development.

Specifically:

### The nf-core repository must be set as the head repo

Personal or organisation repos MUST show on GitHub as being forked from nf-core, not the other way around.
This clarifies where the primary development location is.
It also means that any pull-requests that are created will automatically select the nf-core repository as the target.

When new pipelines are added to nf-core, the repository SHOULD be moved to nf-core instead of being forked.

If you have already forked your pipeline to nf-core, you can [email GitHub support](https://support.github.com/contact?subject=Reroute%20a%20Fork&tags=rr-forks) and request that they reroute the fork.

### Disable GitHub features for forks

To encourage contributors to focus on the nf-core repository, you SHOULD disable GitHub issues, wiki, and projects on your forked repository.
You'll find these options under the GitHub repository settings.
