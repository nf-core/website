---
title: Gitpod
subtitle: Testing, Code review and Website dev.
weight: 1
---

# Introduction

Gitpod is an open-source developer platform that can be spin up a little virtual machine quickly from a Git repository, containing all the programs needed to run your pipelines.
It has many purposes, including: testing pipelines, learning Nextflow, reviewing/editing code and live views of edited markdown/HTML for webpage development.
Any changes within Gitpod can be pushed to branches for review.

Gitpod is free (with a GitHub/GitLab/Bitbucket account) for use up to 50hours/month, with 30GB temporary storage and allows up to 4 parallel workspaces (as of March 2022).

> Gitpod has kindly recognised nf-core as an open-source project, so if you are a member of the [nf-core GitHub organisation](https://github.com/nf-core/) (ask on Slack in the `#github-invitations` channel) then you will automatically be moved onto the _'Professional Open Source plan'_ with unlimited hours.

Gitpod runs in your browser (recommended with Chrome, Firefox or Edge). Their website has [extensive documentation](https://www.gitpod.io/docs).

Each nf-core repository (including the nf-core webpage) should have a ready-to-go Gitpod environment waiting to be started. Since this is a new nf-core feature, check for the `.gitpod.yml` file in the root of the repository.
If it's there, you're ready to start using Gitpod.
If the file is absent, then please [synchronise your pipeline with the latest version of the template](https://nf-co.re/docs/contributing/sync).
Alternatively, if this is not possible, you can add it manually.
Make a fork of the repository, start a new branch, and then in your new branch copy the `.gitpod.yml` file from the [nf-core tools repository](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/.gitpod.yml). Then start Gitpod from this branch.

By clicking on the Gitpod button (see below), your repository is cloned into your Gitpod instance with all the software necessary to run tests and display files. A terminal is also available to run test and workflow debugging.

## Running Gitpod on an nf-core repo

We recommend that you first install a [Gitpod browser extension](https://www.gitpod.io/docs/browser-extension).
This adds the green Gitpod button to each Git repository for easy access to the environment.
This button simply appends `https://gitpod.io/#` onto the git repo URL (e.g. `https://gitpod.io/#https://github.com/nf-core/rnaseq`).
![image of gitpod button](@assets/contributing/gitpod/gitpodbutton.png)

Once you click on the green button, you will asked to sign in with Github.
Then you should see something similar to the following:
![Screenshot of the gitpod interface](@assets/contributing/gitpod/nf-core-gitpod.png)

1. **The sidebar** allows you to customise your environment and perform basic tasks (Copy/Paste, Open files, search, git, etc.)
   Click the Explorer button to see which files are in this repository:
   ![PNG](@assets/contributing/gitpod/explorer.png)

2. **The terminal** allows you to run all the programs in the repository.
   The base image comes with `nextflow`, `nf-core` and `docker` already installed.
   The terminal may not appear automatically, in which case navigate to the top icon of the sidebar, and choose _Terminal_ → _New Terminal_.

3. **The main window** allows you to view and edit files.
   Clicking on a file in the Explorer will open it within the main window.
   Once a file is open, Markdown or HTML can be rendered using the preview option.
   ![PNG](@assets/contributing/gitpod/preview.png)

## Guides

See the following guides for Gitpod use in different scenarios:

- [Using git within Gitpod](git_in_gitpod.md)
- [Configuration of a Gitpod environment](config.md)
- [How to develop Markdown for the nf-core website](webdev.md)
- [How to develop nf core pipelines with Gitpod](nf_core_repo.md)
- [Other Gitpod tutorials](other_tutorials.md)
