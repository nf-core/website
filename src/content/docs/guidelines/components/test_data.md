---
title: Test-data Specifications
subtitle: Specifications for writing nf-core test dataset files
---

## New test-data specifications

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

### 1 General

#### 1.1 Replication of test data

The new test data file within a branch (modules or pipelines) SHOULD NOT replicate existing test-data unless absolutely necessary

- If you need to make a new file that can be generated from an upstream file
- For example, if you need a particular bioinformatic index file for a tool, index an existing FASTA file on the test-datasets branch

#### 1.2 Size of test data

Test data SHOULD be as small as possible

- It cannot exceed the GitHub file maximum
- Data should be sub-sampled as aggressively as possible

#### 1.3 License of test data

Test data MUST be publicly available and have licenses to allow public reuse

#### 1.4 Documentation of test data

Test data files SHOULD be described on the given branch's README file, describing source, how generated, licenses etc.

### 2 Modules

- In order to keep the size of the test data repository as minimal as possible, pre-existing files from [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) MUST be reused if at all possible.

- If the appropriate test data doesn't exist in the `modules` branch of [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) please contact us on the [nf-core Slack `#modules` channel](https://nfcore.slack.com/channels/modules) (you can join with [this invite](https://nf-co.re/join/slack)) to discuss possible options.

- It may not be possible to add test data for some modules e.g. if the input data is too large or requires a local database. In these scenarios, it is recommended to use the Nextflow [`stub`](https://www.nextflow.io/docs/latest/process.html#stub) feature to test the module. Please refer to the [`gtdbtk/classify`](https://github.com/nf-core/modules/blob/79d38a306bdaf07000e0d6f300684d3ed38c8919/modules/gtdbtk/classifywf/main.nf#L66) module and its corresponding [test script](https://github.com/nf-core/modules/blob/79d38a306bdaf07000e0d6f300684d3ed38c8919/tests/modules/gtdbtk/classifywf/main.nf#L20) to understand how to use this feature for your module development.

#### 2.1 Module test data organisation

Files SHOULD be generally organised based on [existing structure](/docs/contributing/test_data_guidelines.md#field-specific-guidance), typically (for bioinformatics pipelines) by discipline, organism, platform or format

#### 2.2 Relatedness of module test data

Downstream or related test-data files SHOULD be named based on the upstream file name

- For example, if you used `genome.fasta` as the upstream file, your output file should be called `genome.<new_extension>`.

#### 2.3 Module test data documentation

Test data files MUST have an entry in the nf-core/test-datasets repo README
