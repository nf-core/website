---
title: Modules
subtitle: Test data guidelines specific to nf-core modules
shortTitle: Modules
weight: 2
---

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Reuse of existing test data

Pre-existing files from [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data) MUST be reused if at all possible to keep the size of the test data repository minimal.

If appropriate test data does not exist in the `modules` branch of [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/modules/data), contact the nf-core community on the [nf-core Slack `#modules` channel](https://nfcore.slack.com/channels/modules) to discuss possible options.

## Test data alternatives for large datasets

Adding test data may not be possible for some modules if the input data is too large or requires a local database.

In these scenarios, use the Nextflow [`stub`](https://www.nextflow.io/docs/latest/process.html#stub) feature to test the module.

Refer to the [`gtdbtk/classify`](https://github.com/nf-core/modules/blob/79d38a306bdaf07000e0d6f300684d3ed38c8919/modules/gtdbtk/classifywf/main.nf#L66) module and its corresponding [test script](https://github.com/nf-core/modules/blob/79d38a306bdaf07000e0d6f300684d3ed38c8919/tests/modules/gtdbtk/classifywf/main.nf#L20) for an example of how to use this feature for module development.

## Module test data organisation

Files SHOULD be organised based on the [existing structure](/docs/tutorials/tests_and_test_data/test_data#field-specific-guidance).

For bioinformatics pipelines, files are typically organised by discipline, organism, platform, or format.

## Relatedness of module test data

Downstream or related test data files SHOULD be named based on the upstream file name.

For example, if `genome.fasta` is used as the upstream file, the output file should be named `genome.<new_extension>`.

## Module test data documentation

Test data files MUST have an entry in the nf-core/test-datasets repository README.
