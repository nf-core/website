---
title: Creating pipelines with nf-core
subtitle: Learn how to use nf-core to build your pipeline
menu:
  main:
    weight: 20
---

> Material originally written for the [Nextflow & nf-core training, October 2022](https://nf-co.re/events/2022/training-october-2022).

## Introduction

Welcome to the nf-core section of the Nextflow & nf-core training course!
You can follow the commands and examples covered within this training in the document below.

### Using GitPod

As with the rest of the Nextflow training up until this point, we will be using GitPod for this training.
We will use a different GitPod environment in order to get the very latest releases of the nf-core tools.
To launch GitPod, follow the link:

<p class="text-center">
  <a href="https://www.gitpod.io/#https://github.com/nf-core/tools" class="btn btn-lg btn-success">
    Launch GitPod
  </a>
</p>

#### Working in an empty directory

To work with a clean directory, you can do the following:

- In the terminal, make a new directory to work in:
  ```bash
  cd ~
  mkdir training
  cd training
  ```
- In the menu top left, select File > Open Folder (<kbd>^</kbd><kbd>⇧</kbd><kbd>O</kbd>)
- Enter `/home/gitpod/training`
- GitPod will probably reload the browser tab
- The file explorer on the left should now have an expandable box with the title `TRAINING`

> 💡 ♻️ As you create new files you should see the explorer populate with them. If you don't see a file you expect, try hovering over the `TRAINING` section title and clicking the <kbd>⟳</kbd> icon that appears to refresh the file list.

### Local installation

If you prefer, you can install the nf-core/tools package locally.
It is a Python package, available via the [Python Package Index](https://pypi.org/project/nf-core/) (`pip`) and [bioconda](https://bioconda.github.io/recipes/nf-core/README.html) (`conda` / `mamba`).
There is a docker container available, however for the purposes of this tutorial it not recommended to use this.

To install from PyPI:

```bash
pip install nf-core
```

To get the very latest development version of the code:

```bash
pip install --upgrade --force-reinstall git+https://github.com/nf-core/tools.git@dev
```

If using conda, first set up Bioconda as described in the [bioconda docs](https://bioconda.github.io/#usage) (especially setting the channel order) and then install nf-core:

```bash
conda install nf-core
```

### First run of nf-core

Whether running using GitPod or locally, you can confirm that nf-core is correctly installed by running `nf-core --help`.
You should get something that looks like the following output:

<!-- RICH-CODEX
command: nf-core --help
img_paths:
  - public_html/assets/markdown_assets/developers/creating_with_nf_core/nfcore_help.svg
-->

![nf-core --help](/assets/markdown_assets/developers/creating_with_nf_core/nfcore_help.svg)

The first set of subcommands are typically useful for running pipelines, the second are for developing pipelines.

You can try out some commands, for example listing available nf-core pipelines:

<!-- RICH-CODEX
command: nf-core list
head: 19
img_paths:
  - public_html/assets/markdown_assets/developers/creating_with_nf_core/nfcore_list.svg
-->

![nf-core list](/assets/markdown_assets/developers/creating_with_nf_core/nfcore_list.svg)

In this tutorial we will focus on creating a pipeline, but please do look at the functionality that nf-core/tools provides to you as a user - especially `nf-core launch` and `nf-core download`.

## Creating a pipeline

To get started with your new pipeline, run the create command:

```bash
nf-core create
```

Although you can provide options on the command line, it's easiest to use the interactive prompts.
Follow the instructions and you should see a new pipeline appear in your file explorer.

Let's move into the new pipeline directory in the terminal:

```bash
# NB: The path will vary according to what you called your pipeline!
cd nf-core-demo/
```

## Nextflow Schema

## nf-core modules
