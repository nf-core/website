---
title: "Training: writing nf-core modules"
subtitle: Workshop material for writing an nf-core module
shortTitle: "Chapter 1: Introduction"
---

:::tip{title = "Last updated"}
Material originally written for the **A practical introduction to nf-core components for Nextflow developers** training at the nf-core Hackathon March 2025, Berlin local site _(see [event](https://nf-co.re/events/2025/hackathon-march-2025/robert-koch-institute))_

Duration: **4hr**

Original author: James A. Fellows Yates (@jfy133), with corrections and improvements from Niklas Schandry (@nschan) and Alexandru Mizeranschi (@amizeranschi)
:::

## Introduction

nf-core modules are a large step towards 'plug and play' Nextflow pipeline development.

An nf-core module is an 'atomic', standardised, reproducible, and already tested [Nextflow DSL2 module](https://www.nextflow.io/docs/latest/module.html).

You can install an nf-core module from a central community repository into a Nextflow pipeline with a single command.
You can then quickly integrate the Nextflow process into your own workflow due to the standardised internal structure of the module, matching upstream outputs and downstream inputs.

So why write or use nf-core modules?

- Efficiency : reduce the number of Nextflow modules you have to make yourself
- Consistency: each module follows the same rules and structure, making it easy to connect with others
- Documentation: every module comes with documentation making it easy to understand to use
- Future-proofing: nf-core modules required to cover all input and output files from the beginning, no more having to re-write your own module when you realise you need extra functionality!
- Automation: consistency and documentation will in the future help automate workflow development

By using these pre-made community-contributed modules, you can vastly e your workflow development time by reducing the need to 're-invent' the wheel for many common bioinformatics tasks.

### Scope

The scope of this training covers:

- Writing nf-core modules
- Using nf-core modules in official nf-core pipelines and custom Nextflow pipelines

### Learning Objectives

By the end of this training, you will:

- Understand the different file components of an nf-core module
- Understand the standards specification behind an nf-core module
- Know where to find documentation about developing an nf-core module
- Know how to write an nf-core module
- Know how to write nf-test tests for an nf-core module
- Understand the workflow development for contributing nf-core modules to the community repository
- Know how to use an nf-core module in any Nextflow pipeline

In the next chapter, we will describe what you will need to have ready to make an nf-core module.

## Chapters

This training has the following chapters

- 📖 Chapter 1: [Introduction](/docs/tutorials/nf-core_training/writing-nf-core-modules/chapter-1-introduction)
- 📖 Chapter 2: [Getting started](/docs/tutorials/nf-core_training/writing-nf-core-modules/chapter-2-getting-started)
- 📖 Chapter 3: [What is a nf-core module?](/docs/tutorials/nf-core_training/writing-nf-core-modules/chapter-3-what-is-a-nf-core-module)
- 📖 Chapter 4: [Generating boilerplate files](/docs/tutorials/nf-core_training/writing-nf-core-modules/chapter-4-generating-boilerplate-files)
- 📖 Chapter 5: [Writing an nf-core module](/docs/tutorials/nf-core_training/writing-nf-core-modules/chapter-5-writing-an-nf-core-module)
- 📖 Chapter 6: [Testing an nf-core module](/docs/tutorials/nf-core_training/writing-nf-core-modules/chapter-6-testing-an-nf-core-module)
- 📖 Chapter 7: [Development workflow](/docs/tutorials/nf-core_training/writing-nf-core-modules/chapter-7-development-workflow)
- 📖 Chapter 8: [Using nf-core modules in pipelines](/docs/tutorials/nf-core_training/writing-nf-core-modules/chapter-8-using-in-pipelines)

:::tip
We recommend reading chapters 2-7 once through _before_ developing the module.
When you get to chapter 7, you can then go back to chapter 2 to start developing the module itself.
:::
