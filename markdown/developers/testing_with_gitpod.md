---
title: Testing with Gitpod
subtitle: Testing, Code review and Website dev.
---

# Introduction

Gitpod is an open-source developer platform that can be spun up in one click from a Git repository. It has many purposes, including: running practise pipelines (e.g. training sessions), editing code and live review of html documentation. We can then push changes to git branches for review. Use of Gitpod is free for use up to 50hours/month (currently). 

There are several topics covered in this documentation page:

1. Downloading and running Gitpod with the nf core repo
2. Configuration of Gitpod
3. Testing nextflow and nf core tools in Gitpod
4. How to do code review
5. How to develop the tools package
6. How to develop for the website
7. How to update the Gitpod environment

### Downloading and running Gitpod with the nf core reposistory

The first step is to install the browser extension: https://www.gitpod.io/docs/browser-extension

Now go to the nfcore repository (https://github.com/nf-core/nf-co.re) and you should find a green Gitpod button in upper right of the screen. If you click on this it will link to an appended URL (https://gitpod.io/#https://github.com/nf-core/nf-co.re), which will open Gitpod and ask you to sign in with either Github,GitLab or Bitbucket. 

Once you have signed in, you should see the following:

<Layout_to_be_completed.img>

The terminal can be used ...
The left hand panel allows you to adjust settings ...
The ...

If you have made changes and want to push these changes to a new branch for review by the nf core team. You next need to click on the source control button (left panel), then choose create new branch and give an informative name. Then click the tick button to commit the code changes to git.

### Configuration of Gitpod

Within a git repository, the main file that controls the gitpod environment is the `.gitpod.yml` file, that contains the instructions on which environment to build and which tools to install. 

Check out the nf core `.gitpod.yml` file https://github.com/nf-core/nf-co.re/blob/master/.gitpod.yml[here]. You can see four main sections:

1. github - allows you to set prebuild gitpod environments for the master or branches. It also allows checks, comments and pull requests.
2. vscode - allows vscode extensions within your environment.
3. ports - opens a port to serve traffic to your on an authenticated URL
4. tasks - this tells gitpod to run particular jobs. In this case we pull a docker container so that we have all of the tools needed in this environment. We can also include lines to download nextflow, nfcore and any other tool you need (which we will explore in a later section).

For more detailed information about these settings, check out the extensive docs at Gitpod (https://www.gitpod.io/docs/config-gitpod-file[here])

### Testing nextflow and nf core tools in Gitpod

To configure a Gitpod environment to test nextflow/nfcore scripts or for use in training events. 

### How to do code review

### How to develop the tools package

### How to develop for the website

### How to update the Gitpod environment