---
title: Testing with Gitpod
subtitle: Testing, Code review and Website dev.
---

# Introduction

Gitpod is an open-source developer platform that can be spin up a little virtual machine quickly from a Git repository, containing all the programs needed to run your pipelines. It has many purposes, including: testing pipelines, learning Nextflow, reviewing/editing code and live views of edited markdown/HTML for webpage development. In addition, within Gitpod, any change can be pushed to branches for review.

Gitpod is free (with a GitHub/GitLab/Bitbucket account) for use up to 50hours/month, with 30GB temporary storage and allows up to 4 parellel workspaces (as of March 2022). It runs in your browser (recommended with Chrome, Firefox or Edge). You can find their extensive documentation [here](https://gitpod.io/). nf-core members are eligible for professional source plans with unlimited hours.

Each nf-core repository (including the nf-core webpage) should have a ready-to-go Gitpod environment waiting to be started. Since this is a new nf-core feature, check for the `.gitpod.yml` file in the root of the repository. If it's there, you're ready to start using Gitpod. If the file is absent, then make a fork of the repository, start a new branch, and then in your new branch copy the `.gitpod.yml` file from the [nf-core tools repository](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/.gitpod.yml). Then start Gitpod from this branch.

By clicking on the Gitpod button (see below), your repository is cloned into your Gitpod instance with all the software necessary to run tests and display files. A terminal is also available to run test and workflow debugging.

## Running Gitpod on an nf-core repo

First install the browser extension, follow the instructions [here](https://www.gitpod.io/docs/browser-extension). This adds the green Gitpod button to each Git repository for easy access to the environment. This button simply appends `https://gitpod.io/#` onto the git repo URL (e.g. https://gitpod.io/#https://github.com/nf-core/rnaseq).

![PNG](/public_html/assets/markdown_assets/developers/gitpod/gitpodbutton.png)

Once you click on the green button, you will asked to sign in with either Github.

Once you have signed in, you should see something similar to the following:

![PNG](/public_html/assets/markdown_assets/developers/gitpod/nf-core-gitpod.png)

1. **The sidebar** allows you to customise your environment and perform basic tasks (Copy/Paste, Open files, search, git, etc.) Click the Explorer button to see which files are in this repository:
   ![PNG](/public_html/assets/markdown_assets/developers/gitpod/explorer.png)

2. **The terminal** allows you to run all the programs in the repository, for example `nextflow`, `nf-core` and `docker` are installed in the nf-core rnaseq repository. The terminal may not appear automatically, in which case navigate to the top icon of the sidebar, and choose Terminal/New Terminal.

3. **The main window** allows you to view and edit files. Clicking on a file in the explorer will open it within the main window. Once a file is open, Markdown or HTML can be rendered using the preview option.

![PNG](/public_html/assets/markdown_assets/developers/gitpod/preview.png)

## Guides

See the following guides for Gitpod use in different scenarios:

1. [Configuration of a Gitpod environment](gitpod/config.md)
2. [Using git within Gitpod](gitpod/git_in_gitpod.md)
3. [How to test a module with pytest in Gitpod](gitpod/pytest.md)
4. [How to develop Markdown for the nf-core website](gitpod/webdev.md)
5. [How to develop nf core pipelines with Gitpod](gitpod/nf_core_repo.md)
6. [Other Gitpod tutorials](gitpod/other.md)
