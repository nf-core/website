---
title: Input/output options
subtitle: Guidelines for defining module inputs and outputs
markdownPlugin: addNumbersToHeadings
shortTitle: Input/output Options
weight: 3
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Required `path` channel inputs

Input channel `path` declarations MUST be defined for all _possible_ input files (i.e., both required and optional files).

- Directly associated auxiliary files to an input file MAY be defined within the same input channel alongside the main input channel (e.g., [BAM and BAI](https://github.com/nf-core/modules/blob/e937c7950af70930d1f34bb961403d9d2aa81c7d/modules/samtools/flagstat/main.nf#L22)).
- Other generic auxiliary files used across different input files MAY be defined using a dedicated input channel (e.g., [reference files](https://github.com/nf-core/modules/blob/3cabc95d0ed8a5a4e07b8f9b1d1f7ff9a70f61e1/modules/bwa/mem/main.nf#L21-L23)).

## Required `val` channel inputs

Input channel `val` declarations SHOULD be defined for all mandatory non-file inputs that are essential for the functioning of the tool (e.g., parameters, flags etc).

- Mandatory non-file inputs are options that the tool MUST have to be able to be run.
- These non-file inputs are typically booleans or strings, and must be documented as such in the corresponding entry in the `meta.yaml`.
- Options, flags, parameters that are _not_ required by the tool to function should NOT be included - rather these can be passed via `ext.args`.

:::info{title="Rationale" collapse}
In 2023, it was decided by a [vote](https://nfcore.slack.com/archives/C043UU89KKQ/p1677581560661679) amongst interested parties to allow non-file mandatory input channels.

It is important to have documented (using the existing display on the website) the bare minimum information required for a module to run.
It also allows module code to consume parameter values without parsing them out of the `ext.args` string and reduces possible risks of entire breakage of modules with future [expected config changes](https://github.com/nextflow-io/nextflow/issues/2723) at a Nextflow level.

Downsides to this approach are readability (now multiple places must be checked on how to modify a module execution - modules.conf `ext.args`, the module invocation in pipeline code etc.), and reduced user freedom.
However it was felt that it was more important for stability in and 'installation' and 'execution' of modules was preferred (e.g., for tools that require position arguments etc.)
:::

:::info{title="Inputs particular cases" collapse}
When one and only one of multiple argument are required:

- If they all are string argument : use 1 argument that will be equal to the string.

  For example, parameter model of [glimpse2 chunk](https://nf-co.re/modules/glimpse2_chunk).

- If some are files put them all in one channel and test if only one is present.

  For example, grouping output parameters of [glimpse2 concordance](https://nf-co.re/modules/glimpse2_concordance):

  ```
  if (((file1 ? 1:0) + (val1 ? 1:0) + (val2 ? 1:0)) != 1) error "One and only one argument required"{:bash}
  ```

  :::

## Output channel emissions

Named file extensions MUST be emitted for ALL output channels. For example, `path "*.txt", emit: txt`.

## Optional inputs

Optional inputs are not currently supported by Nextflow.
However, passing an empty list (`[]`) instead of a file as a module parameter can be used to work around this issue.

For example, having a module (`MY_MODULE`) that can take a `cram` channel and an optional `fasta` channel as input, can be used in the following ways:

```groovy
MY_MODULE(cram, [])     // fasta is optional, the module will run without the fasta present
MY_MODULE(cram, fasta)  // execution of the module will need an element in the fasta channel
```

## Optional outputs

Optional outputs SHOULD be marked as optional:

```groovy
tuple val(meta), path('*.tab'), emit: tab,  optional: true
```

## One output channel per output file type

Each output file type SHOULD be emitted in its own channel (and no more than one), along with the `meta` map if provided ( the exception is the versions.yml ).

In some cases the file format can be different between files of the same type or for the same function (e.g., indices: `.bai` and `.crai`). These different file formats SHOULD be part of the same output channel since they are they serve the same purpose and are mutually exclusive.

```groovy
tuple val(meta), path("*.{bai,crai}"), emit: index
```

:::info{title="Rationale" collapse}
This approach allows specific output types to be identified and accessed within their designated channel.

So when the output definition of module called `SAMTOOLS_INDEX` looks like this:

```groovy
tuple val(meta), path("*.{bai,crai}"), emit: index
```

The output files can be accessed like this:

```groovy
SAMTOOLS_INDEX.out.index
```

Regardless whether they are a `bai` or `crai` as downstream SAMTOOLS modules should accept either without an issue.

:::
