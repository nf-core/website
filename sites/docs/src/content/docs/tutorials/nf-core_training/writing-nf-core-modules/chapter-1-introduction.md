---
title: "Training: writing nf-core modules"
subtitle: Workshop material for writing an nf-core module
---

## Introduction

nf-core modules are a large step towards 'plug and play' Nextflow pipeline development.

An nf-core module is an 'atomic', standardised, reproducible, and already tested [Nextflow DSL2 module](https://www.nextflow.io/docs/latest/module.html) process.
You can install an nf-core module from a central community repository into a Nextflow pipeline with a single command.
Then you can quickly integrate the Nextflow process into your own workflow due to the standardised internal structure of the module; matching upstream process outputs and downstream inputs.

By using these pre-made community contributed modules, you can vastly reduce your workflow development time by reducing the need to 're-invent' the wheel for many common bioinformatics tasks.

### Scope

The scope of this training covers:

- Writing nf-core modules
- Using nf-core modules in official nf-core pipelines and custom Nextflow pipelines

### Learning Objectives

By the end of this training, you will:

- Understand the different file components of an nf-core/module
- Understand the standard of nf-core modules
- Know where to find documentation of developing nf-core modules
- Know how to write an nf-core module
- Know how to write nf-test tests for your nf-core module
- Understand the workflow development for contributing nf-core modules to the community repository
- Know how to use an nf-core pipeline in any Nextflow pipeline

In the next chapter, we will describe what you will need to have ready to make an nf-core module.

## Chapters

This training has the following chapters

- Chapter 1: [Introduction](/docs/tutorials/nf-core_training/writing-nf-core-modules/chapter-1-introduction)
- Chapter 2: [Getting started](/docs/tutorials/nf-core_training/writing-nf-core-modules/chapter-2-getting-started)
- Chapter 3: [What is a nf-core module?](/docs/tutorials/nf-core_training/writing-nf-core-modules/chapter-3-what-is-a-nf-core-module)
- Chapter 4: [Generating boilerplate files](/docs/tutorials/nf-core_training/writing-nf-core-modules/chapter-4-generating-boilerplate-files)
- Chapter 5: [Writing an nf-core module](/docs/tutorials/nf-core_training/writing-nf-core-modules/chapter-5-writing-an-nf-core-module)
- Chapter 6: [Testing an nf-core module](/docs/tutorials/nf-core_training/writing-nf-core-modules/chapter-6-testing-an-nf-core-module)
- Chapter 7: [Development workflow](/docs/tutorials/nf-core_training/writing-nf-core-modules/chapter-7-development-workflow)
- Chapter 8: [Using nf-core modules in pipelines](/docs/tutorials/nf-core_training/writing-nf-core-modules/chapter-8-using-in-pipelines)

We recommend reading chapters 2-7 once through _before_ developing the module.
When you get to chapter 7, you can then go back to chapter 2 to start developing the module itself.
