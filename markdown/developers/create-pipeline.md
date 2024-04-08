---
title: Pipeline Creation Instructions
subtitle: Instructions for creating a sanger-tol pipeline
---

## nf-core instructions

This page is heavily inspired by the nf-core page [Adding a new pipeline](https://nf-co.re/docs/contributing/adding_pipelines).

## Create the pipeline

All pipelines _must_ use the [nf-core template](https://nf-co.re/docs/contributing/guidelines/requirements/use_the_template).
This is done by using the `nf-core create` command - see [the docs](https://nf-co.re/tools#creating-a-new-pipeline) for detailed instructions.
This tool does lots of things for you: it gives you the correct file structure and boiler plate code
and also sets up the required `git` infrastructure for you to keep your pipeline in sync in the future.

When asked _Do you want to customize which parts of the template are used ?_, answer `Y`.
Then, set the "Pipeline prefix" to "sanger-tol", and when asked to "Skip template areas",
disable "iGenomes config" by:

- pressing the down arrow to choose the option,
- pressing the space-bar to select it for disablement,
- pressing Enter

## Push to GitHub

Create an empty repository on GitHub for your new pipeline under <https://github.com/sanger-tol>.
Do this by:

1. going to the GitHub website,
2. clicking `+` then _New Repository_,
3. selecting "sanger-tol" as the Owner.

Make sure _not_ to initialise it with _any_ file, `README` or `LICENSE`: you just want an empty repository.
You already have these files generated from the nf-core template.

Leave the repository as "Public". We don't want to hide our pipelines, even when they're in progress.

Once created, copy the git URL and add this as a remote to your local git repository.
The `nf-core create` command will have initialised a git repository for you,
so all you need to do is add the remote:

```bash
## Add a remote called 'origin' - this is the default name for a primary remote
git remote add origin https://github.com/sanger-tol/PIPELINE_NAME.git
```

The create command also generated the three standard nf-core branches (`master`, `dev` and `TEMPLATE`),
together with an initial commit which is shared between them.
This git structure is required for automatic template synchronisation in the future.

You first need to rename the `master` branch:

```bash
git branch -m master main
```

Then, you can push these new branches to the remote GitHub repository:

```bash
git push --all origin
```

You should now see the vanilla nf-core template and branches in the github.com web interface.

## GitHub configuration

Head up to your repository on GitHub and do the following.

In the About section on the right, click on the cog wheel and:

1. Set the URL to <https://pipelines.tol.sanger.ac.uk/$PIPELINE_NAME>.
2. Add the topics `pipeline` and `nextflow`. This is required to enable it on the pipelines website.
   - Most pipelines also have `workflow` and `genomics`.
3. Enter a description.

Then, ask [@muffato](https://github.com/muffato) or [@mcshane](https://github.com/mcshane) to add the repository to:

1. The ["nextflow_all"](https://github.com/orgs/sanger-tol/teams/nextflow_all) team with the "write" permission
2. The ["nextflow_admin"](https://github.com/orgs/sanger-tol/teams/nextflow_admin) team with the "admin" permission
3. Remove your individual access to the repository

Finally, ask [@gq1](https://github.com/gq1) to set up the pipeline settings via <https://pipelines.tol.sanger.ac.uk/pipeline_health>.

## Other bits

The repository needs to be integrated with Zenodo before making the first release.
Better to do it now before anyone forgets !
Ask [@muffato](https://github.com/muffato) to enable the Zenodo integration.
