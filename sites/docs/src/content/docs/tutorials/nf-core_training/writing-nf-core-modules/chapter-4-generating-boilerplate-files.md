---
title: "Chapter 4: Generating boilerplate files"
subtitle: "How to generate all the file skeletons you need for writing your module"
---

In this chapter, we will explain how to generate and what template files you will need to help you start writing your nf-core module.

To create an nf-core module, you first need to fork and clone to your working environment the [nf-core/modules GitHub repository](https://github/nf-core/modules).
You should then make a new branch corresponding to your new module.

You can then use the nf-core tools to generate boilerplate template files for your specific module, by the following command in the root of the repository:

```bash
nf-core modules create <toolname>/<subcommand>
```

The convention in nf-core modules is that for single-command tools (i.e., a tool that has no subcommand) you can specify just the tool name.
If the tool has subcommands (even if you do not plan to add these), you should include both the tool name, and then the name of the subcommand after a `/`.
In cases where you have a teriary level of subcommands, you can append this to the subcommand.

All parts of the name should be lowercase, and just alphabetic characters or numbers with no punctuation or special characters.

An example of a tool (`drep`) with the subcommand `compare` would be:

```bash
nf-core modules create drep/compare
```

nf-core/tools will try to automate as much set up of the boilerplate template files as possible. For example, it will search [Bioconda](https://bioconda.github.io/) and [biocontainers](https://biocontainers.pro/) for the latest versions of your bioinformatics tool, to pre-specify this information in the container definitions of your module already.

You will then be prompted to provide a few bits of information that will be added to the boilerplate template files for you, typically:

- Your Github user name
- A process resource label
  - These are standardised tags corresponding to certain memory, CPU, and wall time resource defaults
- Whether the module should use a meta map

:::tip
You can also provide some of this information by dedicated command-line flags, allowing you to skip interacting with your console.
:::

Once the command has completed, you should see the following files and directories.

```tree
modules/nf-core/drep/compare/
├── environment.yml
├── main.nf
├── meta.yml
└── tests
    └── main.nf.test
```

If you later create a second module for a second subcommand (`dereplicate`), the directory structure would look like this.

```tree
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

In this chapter we've explained the nf-core tools command you can use to generate boilerplate template files, what files are used, and some naming conventions.

In the next chapter we will go step-by-step through each file and provide guidance what to change and when to complete your nf-core module.
