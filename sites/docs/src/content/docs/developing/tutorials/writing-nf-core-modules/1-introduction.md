---
title: "Chapter 1: Introduction"
subtitle: "Learn to write nf-core modules"
shortTitle: "Chapter 1: Introduction"
---

:::tip{title="Last updated"}
Material originally written for the **A practical introduction to nf-core components for Nextflow developers** training at the nf-core Hackathon March 2025, Berlin local site _(see [event](https://nf-co.re/events/2025/hackathon-march-2025/robert-koch-institute))_.

Duration: **4hr**

Original author: James A. Fellows Yates (@jfy133), with corrections and improvements from Niklas Schandry (@nschan) and Alexandru Mizeranschi (@amizeranschi).
:::

An nf-core module is an atomic, standardised, reproducible, and tested [Nextflow DSL2 module](https://www.nextflow.io/docs/latest/module.html). You can install an nf-core module into any Nextflow pipeline with a single command, then integrate it into your workflow thanks to its consistent internal structure.

## Why use nf-core modules

- **Efficiency**: reuse community modules instead of writing your own.
- **Consistency**: every module follows the same structure, so modules connect together predictably.
- **Documentation**: each module ships with metadata and descriptions.
- **Future-proofing**: modules cover all inputs and outputs from the start, so you avoid rewrites when requirements change.
- **Automation**: standardisation enables future automated workflow development.

## Scope

This training covers:

- Writing nf-core modules.
- Using nf-core modules in official nf-core pipelines and custom Nextflow pipelines.

## Learning objectives

By the end of this training, you will be able to:

- Identify the file components of an nf-core module.
- Apply the nf-core module specifications.
- Locate documentation for developing an nf-core module.
- Write an nf-core module.
- Write nf-test unit tests for an nf-core module.
- Follow the workflow for contributing modules to the community repository.
- Use an nf-core module in any Nextflow pipeline.

## Chapters

- Chapter 1: [Introduction](/docs/developing/tutorials/writing-nf-core-modules/1-introduction)
- Chapter 2: [Getting started](/docs/developing/tutorials/writing-nf-core-modules/2-getting-started)
- Chapter 3: [What is an nf-core module?](/docs/developing/tutorials/writing-nf-core-modules/3-what-is-a-module)
- Chapter 4: [Generating boilerplate files](/docs/developing/tutorials/writing-nf-core-modules/4-boilerplate)
- Chapter 5: [Writing an nf-core module](/docs/developing/tutorials/writing-nf-core-modules/5-writing-modules)
- Chapter 6: [Testing an nf-core module](/docs/developing/tutorials/writing-nf-core-modules/6-testing)
- Chapter 7: [Development workflow](/docs/developing/tutorials/writing-nf-core-modules/7-development)
- Chapter 8: [Using nf-core modules in pipelines](/docs/developing/tutorials/writing-nf-core-modules/8-using)

:::tip
Read chapters 2–7 through once before starting development. When you reach chapter 7, return to chapter 2 to begin writing the module itself.
:::
