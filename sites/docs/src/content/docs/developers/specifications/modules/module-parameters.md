---
title: Module parameters
subtitle: Guidelines for using parameters in nf-core modules
markdownPlugin: addNumbersToHeadings
shortTitle: Parameters
weight: 5
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Module input and outputs

Your module file SHOULD only define input and output files as command-line parameters to be executed within the process.

## Use of parameters within modules

You MUST only initialise and use `params` within the module in the local context of the module.
You MUST NOT assume that named `params` defined in the parent workflow will be passed to the module to allow developers to call their parameters whatever they want.
It may be more suitable to use additional `input` value channels to cater for such scenarios.

## Specification of multiple-threads or cores

If the tool supports multi-threading then you MUST provide the appropriate parameter using the Nextflow `task` variable. For example, `--threads $task.cpus`.

## Evaluation of parameter within a module

You MUST also define within the process any parameters that need to be evaluated in the context of a particular sample. For example, single-end/paired-end data.
