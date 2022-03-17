---
title: Testing with Gitpod
subtitle: Testing, Code review and Website dev.
---

# Introduction

Gitpod is an open-source developer platform that can be spun up quickly from a Git repository, containing all the programs needed to run your pipelines. It has many purposes, including: testing pipelines, learning Nextflow, reviewing/editing code and live views of edited markdown/HTML for webpage development. In addition, within Gitpod, any change can be pushed to branches for review. 

Gitpod is free for use up to 50hours/month (currently). You can find their extensive documentation [here](https://gitpod.io/). 

There are several topics covered in this documentation page:

1. Downloading and running Gitpod on an nf-core repo
2. Configuration of Gitpod
3. Making changes to a repository and commiting to git
4. How to do code review (TBC)
5. How to develop the tools package (TBC)
6. How to develop markdown for the website
7. Other Nextflow tutorials using Gitpod


## Downloading and running Gitpod on an nf-core repo

First install the browser extension: https://www.gitpod.io/docs/browser-extension

Now go to an nfcore repository (e.g. https://github.com/nf-core/rnaseq; all repos should have working gitpod instances) and you should find a green Gitpod button toward the upper right of the screen. This link appends `https://gitpod.io/#` onto the git repo URL (e.g. https://gitpod.io/#https://github.com/nf-core/rnaseq), which will open Gitpod and ask you to sign in with either Github,GitLab or Bitbucket. 

Once you have signed in, you should see the following:

![PNG](/assets/markdown_assets/developers/gitpod/nf-core-gitpod.png)


1. **The terminal** allows you to run all the programs in the repository, for example Nextflow, nf-core tools and docker are installed in the nf-core RNAseq repository.

2. **The sidebar** allows you to customise your environment and perform basic tasks (Copy/Paste, Open files, search, etc.)

3. **The main window** allows you to view and edit files.

Now you can run Nextflow with the chosen nf-core repository. 

For example, for nf-core/rnaseq pipeline we can type the following into the terminal:

    ```console
		nextflow run nf-core/rnaseq \
		-profile test,docker \
		--outdir my_result
    ```

This should run the test data through nf-core rnaseq, using docker with your results in the folder: "my_result". This will take some time to complete.


## Configuration of Gitpod

Within each git repository, the main file that controls the gitpod environment is the `.gitpod.yml` file, that contains the instructions on which environment to build and which tools to install. 

Check out the nf core `.gitpod.yml` file [here](https://github.com/nf-core/nf-co.re/blob/master/.gitpod.yml). You will often see five main sections:

1. **github** - allows configuration of github. e.g. Allows gitpod to create prebuilds for branches.
2. **vscode** - allows vscode extensions within your environment.
3. **ports**  - opens a port to serve traffic to a public URL
4. **tasks**  - this tells gitpod to run particular jobs. Usually you will see the following:

`- init:` sections can be used to install packages as a pre-build, so it doesn't have to run each time you open an environment.

`- command:` sections execute the given lines of code on every workspace startup.

5. **image** - a container image to pull into Gitpod. Many nf-core pipelines use the image `nfcore/gitpod:latest`. This allows the gitpod environment to contain working nextflow and nf-core scripts and other essential tools such as docker. 

For more detailed information about these settings, check out the extensive docs at Gitpod [here](https://www.gitpod.io/docs/config-gitpod-file).


## Making changes to a repository and commiting to git.

Gitpod environments are a handy place to try out Nextflow and nf-core tools, test new features and make suggested changes to the underlying code. 

Once you have made changes, you can push these to a new branch for review by the nf core team. First, make sure your edited files are saved. Then, click on the source control button (on the left hand panel). Choose create new branch and give an informative name. Then click the tick button to commit the code changes to your new git branch.


## How to do code review



## How to develop the tools package



## How to develop markdown for the website

To develop code for the nfcore website (https://github.com/nf-core/nf-co.re) click the green Gitpod button in upper right of the screen of the repo (once you have downloaded the browser extension). Or click the following link (https://gitpod.io/#https://github.com/nf-core/nf-co.re), which will open Gitpod and ask you to sign in with either Github,GitLab or Bitbucket. 

Each page of the nf-core webpage is written in markdown. Gitpod convinienty can render the markdown within the window. In the nf-core repo (linked above), go to Explorer (left hand panel), then navigate to `./markdown/testing_with_gitpod.md`. Click on this file to open the raw code for this page. 

Next, at the top right of the markdown window, click the preview button. This should open up in a new window the rendering of the markdown. Now you can edit your raw code and see the changes happening live in the preview. 

## Other Nextflow tutorials using Gitpod

1. Official Seqera (Nextflow) Training  (Paolo di Tommaso and Evan Floden). -> Link Pending

2. Variant Calling Nextflow Tutorial (Sateesh Peri and Michael Cipriano)  https://sateeshperi.github.io/nextflow_varcal/nextflow/. With a gipod environment at this link:https://gitpod.io/#https://github.com/sateeshperi/nextflow_tutorial.git . This one comes with instructions for running a tutorial completely within Gitpod.