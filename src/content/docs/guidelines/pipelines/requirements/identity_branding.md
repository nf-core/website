---
title: Identity and branding
subtitle: Primary development must on the nf-core organisation.
menu:
  main:
    weight: 5
---

The nf-core community was founded to allow different groups to collaborate on pipelines.
At its inception, the first step was to remove institutional branding from existing workflows and migrate them to the nf-core organisation and nf-core branding.
This was done to remove institutional and organisational ownership of pipelines and code, to encourage collaboration and openness.
In this spirit, we require nf-core pipelines to be _owned_ by nf-core, and nf-core alone.

Whilst pipelines can be forked on GitHub for **personal** development work, the nf-core repository should be the primary source for all development.

Specifically:

### The nf-core repository should be set as the _head_ repo

Your personal / organisation repos should show on GitHub as being forked from nf-core, not the other way around.
This is so that it's clear where the primary development location is.
It also means that any pull-requests that are created will automatically select the nf-core repository as the target.

When new pipelines are added to nf-core, please _move_ the repository to nf-core instead of forking it.

If you have already forked your pipeline to nf-core, you can [email GitHub support](https://support.github.com/contact?subject=Reroute%20a%20Fork&tags=rr-forks) and request that they reroute the fork.

### Disable GitHub features for forks

To encourage contributors to focus on the nf-core repository, please disable GitHub issues / wiki / projects on your forked repository.
You'll find these options under the GitHub repository settings.
