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

#### 2.1 Module test data organisation

Files SHOULD be generally organised based on [existing structure](/docs/contributing/test_data_guidelines.md#field-specific-guidance), typically (for bioinformatics pipelines) by discipline, organism, platform or format

#### 2.2 Relatedness of module test data

Downstream or related test-data files SHOULD be named based on the upstream file name

- For example, if you used `genome.fasta` as the upstream file, your output file should be called `genome.<new_extension>`.

#### 2.3 Module test data documentation

Test data files MUST have an entry in the nf-core/test-datasets repo README
