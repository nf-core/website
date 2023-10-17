---
title: Basic training to create an nf-core pipeline
subtitle: A guide to create Nextflow pipelines using nf-core tools
---

# Introduction

## Objectives / Overview (which is better?)

## Prerequisites

## Follow the training videos

This training can be followed either based on this documentation alone, or via a training video hosted on youtube. You can find the youtube video in the Youtube playlist below:

(no such video yet)

# Using gitpod

For this tutorial we are going to use Gitpod, which is best for first-timers as this platform contains all the programs and data required.
Gitpod will contain a preconfigured Nextflow development environment and has the following requirements:

- A GitHub account
- Web browser (Google Chrome, Firefox)
- Internet connection

Simply click the link and log in using your GitHub account to start the tutorial:

(link to non-existing gitpod)

For more information about gitpod, including how to make your own gitpod environement, see the gitpod bytesize talk on youtube (link to the bytesize talk)

## Explore your Gitpod interface

You should now see something similar to the following:

(insert Gitpod welcome image)

- **The sidebar** allows you to customize your Gitpod environment and perform basic tasks (copy, paste, open files, search, git, etc.). Click the Explorer button to see which files are in this repository.
- **The terminal** allows you to run all the programs in the repository. For example, both `nextflow` and `docker` are installed and can be executed.
- **The main window** allows you to view and edit files. Clicking on a file in the explorer will open it within the main window. You should also see the nf-training material browser (<https://training.nextflow.io/>).

To test that the environment is working correctly, type the following into the terminal:

```bash
nextflow info
```

This should come up with the Nextflow version and runtime information:

```
Version: 22.10.4 build 5836
Created: 09-12-2022 09:58 UTC
System: Linux 5.15.0-47-generic
Runtime: Groovy 3.0.13 on OpenJDK 64-Bit Server VM 17.0.3-internal+0-adhoc..src
Encoding: UTF-8 (UTF-8)
```

## Reopening a Gitpod session

You can reopen an environment from <https://gitpod.io/workspaces>. Find your previous environment in the list, then select the ellipsis (three dots icon) and select Open.

If you have saved the URL for your previous Gitpod environment, you can simply open it in your browser.

Alternatively, you can start a new workspace by following the Gitpod URL: <https://gitpod.io/#https://github.com/nextflow-io/training>

If you have lost your environment, you can find the main scripts used in this tutorial in the `nf-training` directory.

## Saving files from Gitpod to your local machine

To save any file from the explorer panel, right-click the file and select Download.

# Explore nf-core/tools

The nf-core/tools package is already installed in the gitpod environment. Now you can check out which pipelines, subworkflows and modules are available via tools.

# Create a pipeline from template

To get started with your new pipeline, run the create command:

```
nf-core create
Although you can provide options on the command line, it’s easiest to use the interactive prompts.
```

### Pipeline git repo

The nf-core create command has made a fully fledged pipeline for you. Before getting too carried away looking at all of the files, note that it has also initiated a git repository:

```
git status
```

It’s actually created three branches for you:

```
git branch
```

Each have the same initial commit, with the vanilla template:

```
git log
```

This is important, because this shared git history with unmodified nf-core template in the TEMPLATE branch is how the nf-core automated template synchronisation works (see the docs for more details).

The main thing to remember with this is that:

When creating a new repository on https://github.com or equivalent, don’t initialise it - leave it bare and push everything from your local clone
Develop your code on either the master or dev branches and leave TEMPLATE alone.

## Run the new pipeline

The new pipeline should run with Nextflow, right out of the box. Let’s try:

```
cd ../
nextflow run nf-core-demo/ -profile test,docker --outdir test_results
```

## Template code walk through
