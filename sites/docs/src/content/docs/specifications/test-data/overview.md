---
title: Test Data
subtitle: Specifications for managing test data in nf-core
weight: 1
parentWeight: 40
---

Test data is essential for automated testing and continuous integration of nf-core components and pipelines. These specifications define standards for organizing, documenting, and managing test datasets to ensure reproducible testing across the nf-core ecosystem.

All test data used by nf-core components and pipelines is centrally stored in the [nf-core/test-datasets](https://github.com/nf-core/test-datasets) repository. This centralized approach ensures test data remains accessible, properly licensed, and efficiently reused across multiple projects.

While these specifications are mandatory for test data contributed to nf-core repositories, they represent proven best practices for managing test data in any Nextflow workflow development project.

## Test data specifications

The following specifications define standards for managing test data in nf-core:

- **[General](/developers/specifications/test-data/general):** Guidelines for test data replication, file size limits, licensing requirements, and documentation standards to ensure test data is accessible and properly managed.
- **[Modules](/developers/specifications/test-data/modules):** Module-specific guidelines for reusing existing test data, handling large datasets with stub runs, proper organization structure in test-datasets repository, and naming conventions for test files.
