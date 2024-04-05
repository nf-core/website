---
title: Move to nf-core
subtitle: Follow this walkthrough to add a new pipeline to nf-core.
weight: 30
---

## Adding your pipeline to the nf-core organisation

Ok, so you're essentially finished. Your pipeline is written, the tests pass and
you're ready to add your workflow to nf-core.

First, go to the settings of your repository. Under the General page, in the 'Danger Zone' you should have an option to Transfer Ownership. Transfer this to the nf-core organisation.

:::tip
You must make sure you are already a part of the nf-core organisation to allow transferring to nf-core. Alternatively you can add a core-team member to your repository, and ask them to do the transfer you.
:::

Once transferred, go to the transferred repository on nf-core, and make a new fork back to your user name or organisation to continue development on a fork.

:::info
Repositories should be _transferred_ instead of _forking_ to nf-core.
If we fork the original repository to nf-core whenever anyone opens a pull-request, it defaults to going to the original user's fork of the repository, not the nf-core repo.
In this case the only way to fix to request manual detachment from GitHub support.
:::

### Branch setup

All nf-core pipelines use branches called `dev` and `master`.
The `master` branch should contain the code from the latest stable release, `dev` should have the latest development code.
We want people to run the latest development code by default up until the first release.
To do this, we set `dev` as the default repository branch.
After an initial release is created, we set the default branch back to `master` so that the default
action is to run the latest stable release code.

Once you have forked the repository, create a new branch called `dev` for the active development.
In the repository settings, set `dev` to be the default branch.

### Repository setup

Remember to configure the repository on the GitHub website with the following:

- A description, the [https://nf-co.re](https://nf-co.re) URL and lots of keywords!
- Issues enabled, disable Wiki and Projects
- A protected `master` branch that requires review and passing tests
- Write permissions for `nf-core/all` and admin permissions for `nf-core/admin`

You can check that all of these settings are done correctly by referring to your pipeline
in the nf-core [Repository health web page](https://nf-co.re/pipeline_health).
This reports the status of various checks and also has the option of fixing errors for you via the GitHub API.

### Differences to your own fork

The main difference when working with the main nf-core fork of your workflow is
that tests for pull-requests against the `master` branch will fail. This is because
the `master` branch should only ever contain code from the last release.
Instead, use the `dev` branch for new work and always make pull-requests against
that. Then the tests should pass.
