---
title: "Chapter 3: What is a nf-core module?"
subtitle: "Understanding the different file components of an nf-core/module"
shortTitle: "Chapter 3: What is a module?"
---

An nf-core module is an opinionated open-source Nextflow 'wrapper' around command-line tool or script.
A collection of these wrappers are stored in a [central community repository](https://nf-co.re/modules) of more than 1000 modules that can be used by anyone in their own pipelines.

This chapter describes what a nf-core module is, and what are special about them compared to a standard Nextflow process.

nf-core modules (mostly) follow the rule of one tool command or a tool plus a single subcommand per wrapper.

Each 'wrapper' consists of a set of files where the structure and contents of each which follow strict specifications.
These specifications are designed based on a consensus of the wider nf-core community.
This specifications span multiple areas, such as:

- Naming and formatting conventions
- What software environment to use
- Use standardised tags
- What input and output files are allowed
- How to allow pipeline developers to use specific tool parameters
- Use of ['meta maps'](https://nf-co.re/docs/contributing/components/meta_map)
- Use of [stubs](https://www.nextflow.io/docs/latest/process.html#stub)

and many more.

Furthermore, the nf-core specifications also expands the number files from simply just a single Nextflow `.nf` script file, due to the nf-core initiative's focus on standardisation, reproducibility, and documentation.

For example, a 'complete' nf-core/module directory can consist of the following 6 files.

```tree
├── environment.yml
├── main.nf
├── meta.yml
└── tests/
    ├── main.nf.test
    ├── main.nf.test.snap
    └── nextflow.config
```

These can be split up into three main categories:

- Module execution
- Module documentation
- Module testing

The relationship between these files is shown in the following diagram:

![Diagram of three main categories of files in an nf-core module. The central 'execution' box includes two files icons, one called main.nf and an arrow pointing into it from a second file icon called 'environment.yml'. In the left 'Documentation' box an arrow points from the previous main.nf to a file in this section called 'meta.yaml'. On the right 'Testing' box, there are three file icons, an icon called 'main.nf.test' with an arrow coming from the execution section's main.nf, an arrow pointing into the main.nf.test icon from another file icon called 'nextflow.config', and an arrow pointing to a final file icon from the 'main.nf.test' to the final icon 'main.nf.test.snap'.](/images/tutorials/training/nf-core-module-file-relationship.png)

## Summary

In this chapter we've briefly introduced what distinguishes an nf-core module versus a typical Nextflow process, and given an overview of what it looks like.

In the next chapter we will explain how to use nf-core tooling to generate partially completed boilerplate template files you can use to speed up the development of your nf-core module.
