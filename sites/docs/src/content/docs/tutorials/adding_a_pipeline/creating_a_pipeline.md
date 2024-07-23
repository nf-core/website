---
title: Create a new pipeline
subtitle: Follow this walkthrough to add a new pipeline to nf-core.
weight: 10
type: 'tutorial'
---

# Create a pipeline from the template

## Request a new pipeline

Before you get too carried away, if you want to share your pipeline with the community, your first task is to propose the new pipeline in the nf-core Slack [`#new-pipelines`](https://nfcore.slack.com/archives/CE6SDEDAA) channel.

There, you will find a _Workflow_ bookmarked to top of the Slack window called _Pipeline proposal_.
This workflow creates a form for you to describe key information about the pipeline you are planning to write.

The proposal will be discussed and checked for uniqueness (i.e., not too much overlap with any existing pipelines)
and added to our [new pipeline tracking board](https://github.com/orgs/nf-core/projects/35/) on GitHub.

Once accepted, the core team will create a new Slack channel for your pipeline.

At this point, you should also join the [#pipeline-maintainers](https://nfcore.slack.com/channels/pipeline-maintainers) channel to keep up to date with announcements and general discussion about pipeline development.

## Create the pipeline

Start by making a new pipeline locally and working with it on your own GitHub account.
Once you have a version of the pipeline that runs, ask the core team to _move_ the repository to the nf-core GitHub organisation for final development, review, and release.

Pipeline repositories are not created within the nf-core organisation from the start in case development takes longer than expected.
This avoids having a lot of unfinished pipelines listed in nf-core.
See [Adding your pipeline to the nf-core organisation](/docs/tutorials/adding_a_pipeline/move_to_nf-core_org) for further information about this process.

All nf-core pipelines must start from the [nf-core template](/docs/guidelines/pipelines/requirements/use_the_template).
The template is created using the `nf-core create` command and will create a base pipeline with the correct file structure, boiler plate code, and `git` infrastructure for you to help keep your pipeline in sync.
See [the docs](/docs/nf-core-tools/pipelines/create) for further information about this nf-core tools command.

Even if you already have a working Nextflow pipeline that was developed independently, it is usually easier to start fresh and copy your code into the template.
Some exceptions can be made if it is not possible for you to start your pipeline again using the template.
However, you will need to set up [manual synchronisation](/docs/tutorials/sync/overview), and this is not for the faint hearted!
Ask the core team on Slack if you would like some advice.

Once the template for your pipeline is created, make sure to switch branch to the `dev` branch with the `git checkout dev` command.
All development should happen on `dev` (or on other feature branches that get merged into `dev`).

:::note
Pipeline names must be all lower-case and contain no punctuation.</br>
This is to allow consistent names between platforms.
:::

## Push to GitHub

Create an empty repository on GitHub for your new pipeline under your personal account by going to the GitHub website and clicking _+_, then _New Repository_.

Make sure _not_ to initialise your new repository with _any_ file, `README`, or `LICENSE`, as these are already included in nf-core template.

Once created, copy the git URL and add it as a remote to your local git repository.
The `nf-core create` command will have initialised a git repository for you,
so all you need to do is add the remote:

```bash
## Add a remote called 'origin' - this is the default name for a primary remote
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPOSITORY.git
```

The `nf-core create` command also generates three standard nf-core branches (`master`, `dev`, and `TEMPLATE`) with an initial shared commit.
This git structure is required for automatic template synchronisation in the future.

You need to push all new branches to the remote GitHub repository:

```bash
git push --all origin
```

You should now see the vanilla nf-core template and branches in your new GitHub repository.

## Work on your pipeline

You are now ready to start your pipeline development!

It is recommended that you follow good git development practices and work on the `dev` or feature branches, committing, and pushing code regularly.

Remember to run the `nf-core lint` command (see [docs](/docs/nf-core-tools/pipelines/lint))
to ensure that your pipeline passes all of the nf-core compatibility tests.
The automated tests on Github Actions also run the linting command, so if something isn't correct you will receive a
GitHub notification.
