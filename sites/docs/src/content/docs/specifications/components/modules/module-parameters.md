---
title: Module parameters
subtitle: Guidelines for using parameters in nf-core modules
markdownPlugin: addNumbersToHeadings
shortTitle: Parameters
weight: 5
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Module input and outputs

The module file SHOULD only define input and output files as command-line parameters to be executed within the process.

## Use of parameters within modules

Initialise and use `params` only within the module in the local context of the module.
Do not assume that named `params` defined in the parent workflow will be passed to the module to allow developers to call their parameters whatever they want.
Use additional `input` value channels for such scenarios.

## Specification of multiple-threads or cores

If the tool supports multi-threading, provide the appropriate parameter using the Nextflow `task` variable.
For example, `--threads $task.cpus`.

## Evaluation of parameter within a module

Define within the process any parameters that need to be evaluated in the context of a particular sample.
For example, single-end/paired-end data.
