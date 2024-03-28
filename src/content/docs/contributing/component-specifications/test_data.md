---
title: Test-data Specifications
subtitle: Specifications for writing nf-core test dataset files
---

## New test-data specifications

The key words "MUST", "MUST NOT", "SHOULD", etc. are to be interpreted as described in [RFC 2119](https://tools.ietf.org/html/rfc2119).

### General

1. The new test data file within a branch (modules or pipelines) SHOULD NOT replicate existing test-data unless absolutely necessary
   - If you need to make a new file that can be generated from an upstream file
   - For example, if you need a particular bioinformatic index file for a tool, index an existing FASTA file on the test-datasets branch
2. Test data SHOULD be as small as possible
   - It cannot exceed the GitHub file maximum
   - Data should be sub-sampled as aggressively as possible
3. Test data MUST be publically available and have licenses to allow public reuse
4. Test data files SHOULD be described on the given branch's README file, describing source, how generated, licenses etc.

### Modules

1. Files SHOULD be generally organised based on [existing structure](/docs/contributing/test_data_guidelines.md#field-specific-guidance), typically (for bioinformatics pipelines) by discipline, organism, platform or format
2. Test-data files SHOULD be named based on the upstream file name
   - For example, if you used `genome.fasta` as the upstream file, your output file should be called `genome.<new_extension>`.
3. Test data files MUST have an entry in the nf-core/test-datasets repo README
