---
title: General
subtitle: General guidelines for nf-core test data
shortTitle: General
weight: 2
---

The keywords "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

## Location of test data

Test data MUST be stored in the [nf-core/test-datasets](https://github.com/nf-core/test-datasets) repository.
Full-size test data for pipeline `test_full` profiles MAY instead be loaded from stable links to the original public source.

Test data files MUST NOT be committed to pipeline repositories.

## Review of new pipeline test data

New test data for a pipeline MUST be added to the pipeline's test-datasets branch through a pull request.

The pull request MUST be reviewed by a member of the nf-core core or maintainers team who is not part of the pipeline's development team.

New pipelines SHOULD add and review their test data before pipeline development starts.

:::info{title="Rationale" collapse}
A review of the test data at the start of development catches problems such as unclear licences or identifiable human data before the pipeline and its tests depend on the data.
A reviewer from outside the development team gives an independent check, which matters most when a small team develops a pipeline.
:::

## Replication of test data

New test data files within a branch (modules or pipelines) SHOULD NOT replicate existing test data unless absolutely necessary.

Files that can be generated from upstream files SHOULD be derived from existing test data on the test-datasets branch.
For example, if a particular bioinformatic index file is needed for a tool, index an existing FASTA file from the test-datasets branch rather than uploading a new index file.

:::info{title="Rationale" collapse}
CI tests for nf-core modules, subworkflows, or pipelines are not required to produce meaningful output.

The main goal of nf-core CI tests is to ensure that a given tool executes without errors.

A test may produce nonsense output or find nothing, as long as the tool does not crash or produce an error.

You SHOULD reuse existing test data as much as possible to reduce the size of the test dataset repository.

You SHOULD upload new test data only if there is no other option within the existing test-data archive.
:::

## Size of test data

Test data SHOULD be as small as possible.

Test data files MUST NOT exceed the GitHub file size maximum.

Data SHOULD be sub-sampled as aggressively as possible while still allowing the tool to execute successfully.

## License of test data

Test data MUST be publicly available and have licenses that allow public reuse.

## Personal and sensitive data

Test data MUST NOT contain personal data, including human sequence data from individual donors, unless the donors consented to open public sharing of that data (for example, [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) samples).

Simulated data and public reference genomes MAY be used.

## Documentation of test data

Test data files SHOULD be described in the given branch's README file.

The README SHOULD include the source of the data, how it was generated, and license information.
