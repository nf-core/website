---
title: Move to nf-core
subtitle: Follow this walkthrough to add a new pipeline to nf-core.
weight: 30
type: 'tutorial'
---

## Adding your pipeline to the nf-core organisation

Ok, so you're essentially finished. Your pipeline is written, the tests pass and
you're ready to add your workflow to nf-core.

First, go to the settings of your repository. Under the _General_ page, in the _Danger Zone_ you there is an option to _Transfer Ownership_. Click this option and transfer the ownership of the pipeline to the nf-core organisation.

:::tip
You must make sure you are already a part of the nf-core organisation to enable transferring your pipeline to nf-core. Alternatively, you can add a core-team member to your repository and ask them to do the transfer the pipeline for you.
:::

Once transferred, go to the transferred repository on nf-core and make a new fork back to your user name or organisation to continue development.

:::info
Repositories should be _transferred_ instead of _forking_ to nf-core.
If a repository if forked to nf-core whenever a pull-request is made it will default to the original user's fork of the repository, not the nf-core repository.
A manual detachment must be requested from GitHub support to resolve this issue.
:::

### Branch setup

All nf-core pipelines use branches called `dev` and `master`.
The `master` branch should contain the code from the latest stable release. The `dev` branch should contain the latest development code.
Pipelines should default to the latest development code by default up until the first release.
To do this, the `dev` should be set as the default repository branch.
After an initial release, the default branch can be set back to `master` so that the default
action is to run the latest stable release.

Once you have forked the repository back to your own repository, create a new branch called `dev` for the active development.
In the repository settings, set `dev` to be the default branch.

### Repository setup

Remember to configure the repository on the GitHub website with the following:

- A description, the [https://nf-co.re](https://nf-co.re) URL and lots of keywords!
- Issues enabled
- Wiki and Projects disabled
- A protected `master` branch that requires review and passing tests
- Write permissions for `nf-core/all` and admin permissions for `nf-core/admin`

:::note
You can check that these settings correct by referring to your pipeline
in the nf-core [Repository health web page](https://nf-co.re/pipeline_health).
This reports the status of various checks and also has the option of fixing errors for you via the GitHub API.
:::

### Differences to your own fork

The main difference when working with the main nf-core fork of your workflow is
that tests for pull-requests against the `master` branch will fail. This is because
the `master` branch should only ever contain code from the last release.
Instead, use the `dev` branch for new work and always make pull-requests against
it to allow tests to pass.
