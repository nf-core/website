---
title: "Chapter 4: Generating boilerplate files"
subtitle: "How to generate all the file skeletons you need for writing your module"
shortTitle: "Chapter 4: Boilerplate files"
---

## Introduction

In this chapter, we will explain how to generate an nf-core module template and what template files you will need, to help you start writing your nf-core module.

## Preparation

To create an nf-core module, you first need to fork and clone to your working environment the [nf-core/modules GitHub repository](https://github/nf-core/modules).

You should then make a new git branch corresponding to your new module.

```bash
git switch master ## to make sure we have the latest state of the repository
git switch -c <MY_NEW_MODULE_BRANCH>
```

## Using nf-core tools to create the boilerplate files

From within your new branch, you can then use the nf-core tools to generate boilerplate template files for your specific module.
To do this, run the following command in the root of the repository:

```bash
nf-core modules create <toolname>/<subcommand>
```

The convention in nf-core modules is that for:

- 1️⃣ Single-command tools (i.e., a tool that has no subcommand) you can specify just the tool name.
  - Example: a software that is executed with `fastp -i <input> -o <output>` that has no subcommands would be created with `nf-core modules create fastp`
- 2️⃣ Tool with subcommands (even if you do not plan to add these) should include both the tool name, and then the name of the subcommand after a `/`.
  - Example: a software that is executed with `samtools view -i <input> -o <output>` would be created with `nf-core modules create samtools/view`
- 3️⃣ In cases where you have a teriary level of subcommands, you can append this to the subcommand.
  - Example: a tool that executed with `samtools view flagstats -i <input> -o <output>` would be created with `nf-core modules create samtools/viewflagstat`

Another important standard of an nf-core module is how to format names.
Within the command, all parts of the module name should be lowercase, and just alphabetic characters or numbers with no punctuation or special characters.

An example of a tool (`drep`) with the subcommand `compare` would be:

```bash
nf-core modules create drep/compare
```

nf-core/tools will try to automate as much set up of the boilerplate template files as possible.
For example, it will try to search [Bioconda](https://bioconda.github.io/) and [biocontainers](https://biocontainers.pro/) for the latest versions of your bioinformatics tool, to pre-specify this information in the container definitions of your module already.

You will then be prompted to provide a few bits of information that will be added to the boilerplate template files for you, typically:

- ❓ Your Github user name
- ❓ A process resource label
  - These are standardised tags corresponding to certain memory, CPU, and wall time resource defaults
  - You should try and match your process resource label to the expected resource usage of your tool as much as you can.
  - You can see the nf-core default resources per label [here](https://github.com/nf-core/tools/blob/52e810986e382972ffad0aab28e94f828ffd509b/nf_core/pipeline-template/conf/base.config#L29-L54), however each pipeline may customise the values of each resource
- ❓ Whether the module should use a [meta map](https://nf-co.re/docs/contributing/components/meta_map)
  - In most cases, we recommend you say 'yes' for this as this allows you carry additional meta information about the sample of file that is being processed
  - This meta information can then be used for more automated additional processing decisions by the pipeline

:::tip
You can also provide some of this information by dedicated command-line flags, allowing you to skip interacting with your console.
:::

## Output of `nf-core modules create`

Once the command has completed, you should see the following files and directories.

```tree {8-13}
modules/nf-core/drep/compare/
├── environment.yml
├── main.nf
├── meta.yml
└── tests
    └── main.nf.test
```

If you later create a second module for a second subcommand (`dereplicate`), the directory structure would look like this.

```tree {8-13}
modules/nf-core/drep/
├── compare
│   ├── environment.yml
│   ├── main.nf
│   ├── meta.yml
│   └── tests
│       └── main.nf.test
└── dereplicate
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
        └── main.nf.test
```

## Summary

In this chapter we've gone through the nf-core tools command you can use to generate boilerplate template files, given an overview of what files are generated, and provided some naming conventions.

In the next chapter we will go step-by-step through each file and provide guidance what to change and when to complete your nf-core module.
