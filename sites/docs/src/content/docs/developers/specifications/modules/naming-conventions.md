---
title: Naming conventions
subtitle: Naming standards for nf-core module files, processes, and channels
markdownPlugin: addNumbersToHeadings
shortTitle: Naming conventions
weight: 2
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Name format of module files

You must use all lowercase for the directory structure of your module name, without punctuation. For example, [`modules/nf-core/bwa/mem/`](https://github.com/nf-core/modules/tree/master/modules/nf-core/bwa/mem/).

You MUST make the name of the software (i.e. `bwa`) and tool (i.e. `mem`) all one word.

The nf-core/tools command will validate your suggested name.

## Name format of module processes

You MUST make the process name in the module file all uppercase. For example, `process BWA_MEM {`. You MUST make the name of the software (i.e. `BWA`) and tool (i.e. `MEM`) all one word separated by an underscore.

## Name format of module parameters

Your parameter names MUST follow the `snake_case` convention.

## Name format of module functions

Your function names MUST follow the `camelCase` convention.

## Name format of module channels

Your channel names MUST follow `snake_case` convention and be all lower case.

## Command file output naming

Your output file (and/or directory) names SHOULD consist of only `${prefix}` and the file-format suffix (for example, `${prefix}.fq.gz` or `${prefix}.bam`).

- This provides re-usability, giving developers flexibility to name their output files when using the module.
- As a result of using this syntax, if the module could _potentially_ have the same named inputs and outputs, you should add a line in the `script` section like below (another example [here](https://github.com/nf-core/modules/blob/e20e57f90b6787ac9a010a980cf6ea98bd990046/modules/lima/main.nf#L37)) which will raise an error asking the developer to change the `ext.prefix` variable to rename the output files so they don't clash.

  ```groovy
  script:
  if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
  ```

- If the input and output files are likely to have the same name, you may set an appropriate default prefix, for example:

  ```nextflow
  def prefix = task.ext.prefix ?: "${meta.id}_sorted"
  ```
