---
title: Creating a new pipeline
subtitle: Scaffolding a new pipeline using the nf-core template
shortTitle: create
weight: 60
---

The `create` subcommand makes a new pipeline using the nf-core base template.
With a given pipeline name, description and author, it makes a starter pipeline which follows nf-core best practices.

After creating the files, the command initialises the folder as a git repository and makes an initial commit.
It also allows the creation of a GitHub repository where this commit is pushed.
This first "vanilla" commit which is identical to the output from the templating tool is important, as it allows us to keep your pipeline in sync with the base template in the future.
See the [nf-core syncing docs](https://nf-co.re/developers/sync) for more information.

<!-- RICH-CODEX
working_dir: tmp
-->

![` nf-core pipelines create -n nextbigthing -d "This pipeline analyses data from the next big omics technique" -a "Big Steve"`](/images/tools/nf-core-create.svg)

Please see the [nf-core documentation](/docs/tutorials/adding_a_pipeline/overview) for a full walkthrough of how to create a new nf-core workflow.

> [!TIP]
> As the command documentation says, remember to come and discuss your idea for a pipeline as early as possible!
> See the [documentation](/docs/tutorials/adding_a_pipeline/move_to_nf-core_org) for instructions.

Note that if the required arguments for `nf-core pipelines create` are not given, it will start a user interface to create the pipeline interactively. If you prefer, you can supply them as command line arguments. See `nf-core pipelines create --help` for more information.

## Customization

The `nf-core pipelines create` command comes with a number of options that allow you to customize the creation of a pipeline if you intend to not publish it as an
nf-core pipeline. This can be done in two ways: by using the graphical interface, or by supplying a `template.yml` file using the `--template-yaml <file>` option.
Both options allow you to specify a custom pipeline prefix to use instead of the common `nf-core`, as well as selecting parts of the template to be excluded during pipeline creation.
The interface app will guide you through the pipeline creation process.

### Creation of a template.yml

The `template.yml` file can include different options to provide to the pipeline creation command.

- `name`: (str) The name of your pipeline [required if not provided through command line options]
- `description`: (str) The description of your pipeline [required if not provided through command line options]
- `author`: (str) The author(s) of the pipeline [required if not provided through command line options]
- `org`: (str) The name of your organisation. Use `nf-core` if you want to submit your pipeline to nf-core [default: nf-core]
- `version`: (str) The initial version of the pipeline [default: 1.0.0]
- `force`: (bool) If the pipeline exists, overwrite it [default: False]
- `outdir`: (str) Output directory where the pipeline will be created [default: .]
- `is_nfcore`: (bool) If the pipeline is an nf-core or custom pipeline [default: True]
- `skip_features`: (list) List of template features which can be skipped

Skippable template features for both nf-core and custom pipelines:

- `igenomes`: removes pipeline options related to iGenomes. Including the `conf/igenomes.config` file and all references to it.

Skipable template features for custom pipelines:

- `github`: removed all files required for GitHub hosting of the pipeline. Specifically, the `.github` folder and `.gitignore` file.
- `ci`: removes the GitHub continuous integration tests from the pipeline. Specifically, the `.github/workflows/` folder.
- `github_badges`: removes GitHub badges from the `README.md` file.
- `nf_core_configs`: excludes `nf_core/configs` repository options, which make multiple config profiles for various institutional clusters available.

:::tip
If you have an old pipeline created with the `nf-core pipelines create` command and you would like to skip some features, you can add this options to the `.nf-core.yml` file in the root directory of your pipeline.
:::

An example of a `template.yml` file is shown below.

```yaml
name: coolpipe
description: A cool pipeline
author: me
org: myorg
skip_features:
  - github
  - ci
  - github_badges
  - igenomes
  - nf_core_configs
```

This will create a pipeline called `coolpipe` in the directory `myorg-coolpipe` (`<prefix>-<name>`) with `me` as the author. It will exclude all possible parts of the template.
