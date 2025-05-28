---
title: "Using subworkflows with modules from multiple organisations"
subtitle: Guidance on how to use nf-core/tools with subworkflows that use modules from more than one organisation.
---

In order to use cross-organisational subworkflows, you can specify individual remotes to be used for each module within the `meta.yml`
file in that subworkflow.

For example, let's say you have a subworkflow in your own organisation's `modules` repository that uses FastQC and fastp.
But, instead of using the fastp module from your own organisation, you want to use nf-core's fastp module. This is the case for the `fastq_trim_fastp_fastqc` subworkflow in the [nf-core-test](https://github.com/nf-core-test/modules/tree/main/subworkflows/nf-core-test/fastq_trim_fastp_fastqc) organization.

In order to have this subworkflow, with a FastQC module from your organisation and fastp from nf-core, you'd define the
components section of this subworkflow's `meta.yml` file as such:

```yaml title="meta.yml" {11-12}
name: "fastq_trim_fastp_fastqc"
description: Read QC, fastp trimming and read qc
keywords:
  - qc
  - quality_control
  - adapters
  - trimming
  - fastq
components:
  - fastqc
  - fastp:
      git_remote: https://github.com/nf-core/modules.git
```

By using the structure above, `nf-core/tools` knows that it should install FastQC from the same repository that the subworkflow is located
in, but fastp should come from the `nf-core/modules` remote.

:::note
Be aware that to install subworkflows from organisations other than nf-core you must
specify `--git-remote` in the `nf-core subworkflows` command, e.g. `nf-core subworkflows install fastq_trim_fastp_fastqc --git-remote $ORG_NAME{:bash}`
:::

If you're using cross-organisational subworkflows in your repository, be aware that you must also specify a different JSON schema
to lint the `meta.yml` files for your subworkflows.
To do this, just change the `components` section type to allow both objects and strings.

```json title="yaml-schema.json" {5}
"components": {
  "type": "array",
  "description": "Modules and subworkflows used in the subworkflow",
  "items": {
    "type": ["object", "string"]
  },
  "minItems": 0
},
```

:::tip
To get the full content for the schema, see the corresponding
file in the [nf-core-test](https://github.com/nf-core-test/modules/blob/main/subworkflows/yaml-schema.json) modules repository.
:::
