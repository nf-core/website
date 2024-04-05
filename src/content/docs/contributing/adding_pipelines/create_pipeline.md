---
title: Adding a new pipeline
subtitle: Follow this walkthrough to add a new pipeline to nf-core.
---

# Create a pipeline from the template

## Request a new pipeline

Before you get too carried away, the first task is to request the new pipeline in the nf-core Slack.
We have a Slack channel dedicated to this: `#new-pipelines`.

There, you will find a _Workflow_ bookmarked to top of the Slack window called _Pipeline proposal_.
This gives you a form to fill in with key information about the pipeline you want to write.

The proposal will be discussed and checked for uniqueness (not too much overlap with any existing pipelines)
and added to our [new pipeline tracking board](https://github.com/orgs/nf-core/projects/35/) on GitHub.

Once accepted, someone from the core team will create a Slack channel for your pipeline and you can get started on the next steps.

You should also at this point additionally join the [#pipeline-maintainers](https://nfcore.slack.com/channels/pipeline-maintainers) channel for major change announcements as well as general discussion on pipeline development related topics.

## Create the pipeline

You'll start by making a new pipeline locally and working with it on your own GitHub account.
Once you have a version of the pipeline that runs, ask the core team to _move_ the repo to the nf-core GitHub organisation for final development, review and release.

We generally don't create repositories within the nf-core organisation from the start, in case development takes longer than expected.
This way we avoid having a lot of "empty" pipelines listed which are basically just the template.
See [Adding your pipeline to the nf-core organisation](/docs/contributing/adding_pipelines/move-to-nf-core-organisation) below for details on this process.

All nf-core pipelines [_must_ use the nf-core template](https://nf-co.re/docs/contributing/guidelines/requirements/use_the_template).
This is done by using the `nf-core create` command - see [the docs](https://nf-co.re/tools#creating-a-new-pipeline) for detailed instructions.
This tool does lots of things for you: it gives you the correct file structure and boiler plate code
and also sets up the required `git` infrastructure for you to keep your pipeline in sync in the future.

If you already have a working Nextflow pipeline that you're starting from, it's usually easier in the long run to start fresh with the nf-core template and copy over your code in the relevant places.
Some exceptions can be made - ask the core team on Slack if you're unsure.
You'll need to set up [manual synchronisation](sync.md), not for the faint hearted!

Once the template for your pipeline is created, make sure to switch branch to the `dev` branch with `git checkout dev`.
All development should happen on dev (or on other branches that get merged into dev).

:::note
Pipeline names must be all lower-case and contain no punctuation.</br>
This is to allow consistent names between platforms.
:::

## Push to GitHub

Create an empty repository on GitHub for your new pipeline under your personal account.
Do this by going to the GitHub website and clicking + then _New Repository_.

Make sure _not_ to initialise it with _any_ file, `README` or `LICENSE`: you just want an empty repository.
You already have these files generated from the nf-core template.

Once created, copy the git URL and add this as a remote to your local git repository.
The `nf-core create` command will have initialised a git repository for you,
so all you need to do is add the remote:

```bash
## Add a remote called 'origin' - this is the default name for a primary remote
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPOSITORY.git
```

The create command also generated the three standard nf-core branches (`master`, `dev` and `TEMPLATE`),
together with an initial commit which is shared between them.
This git structure is required for automatic template synchronisation in the future.

You need to push these new branches to the remote GitHub repository:

```bash
git push --all origin
```

You should now see the vanilla nf-core template and branches in the github.com web interface.

## Work on your pipeline

Ok, now you're all set with your own personal nf-core pipeline!
You can now start writing code for real.

Follow usual git development practices, working on the `dev` branch and committing + pushing code as normal.

Remember to run the `nf-core lint` command (see [docs](https://nf-co.re/tools#linting-a-workflow))
to make sure that your workflow passes all of the nf-core compatibility tests.
The automated tests on Github Actions also run this, so you should get a
notification from GitHub if something breaks.

When testing the pipeline you can add the `debug` profile (`-profile debug`) to the Nextflow command line,
to enable warnings about process selectors, show additional debug output and disable cleanup.
