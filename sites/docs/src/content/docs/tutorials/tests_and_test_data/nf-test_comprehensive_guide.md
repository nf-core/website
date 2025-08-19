---
title: "nf-test: Comprehensive Testing Guide"
subtitle: Complete guide to testing nf-core components, subworkflows, and pipelines with nf-test
shortTitle: nf-test Comprehensive Guide
weight: 5
---

## Overview

This comprehensive guide covers testing nf-core components using nf-test.
From writing your first test to implementing advanced patterns, this guide provides everything needed to create robust, maintainable tests for modules, subworkflows, and pipelines.

## Guide structure

### Getting started

1. [Installation](./components/01_installation.md) - Setting up nf-test in your development environment
2. [Project setup](./components/02_project_setup.md) - Configuring your nf-core pipeline repository for testing with nf-test

### Component Testing

3. [Testing modules](./components/03_testing_modules.md) - Testing individual nf-core modules
4. [Testing subworkflows](./components/04_testing_subworkflows.md) - Testing nf-core subworkflows
5. [Testing pipelines](./components/05_testing_pipelines.md) - End-to-end pipeline testing

- [nf-test assertions](./components/06_assertions.md) - Comprehensive assertion patterns and verification techniques

### Data management and integration

7. [Test data management](./components/07_test_data_management.md) - Organizing and managing test datasets
8. [CI/CD integration](./components/08_cicd_integration.md) - Integrating nf-test with continuous integration

### Commands and reference

9. [nf-test commands and integration](./components/09_commands_integration.md) - Essential commands and nf-core integration

### Troubleshooting & Best Practices

10. [FAQ and Debugging](./components/10_faq_debugging.md) - Best practices, common issues, and solutions for nf-test

## nf-core testing principles

nf-core testing follows these key principles:

1. Comprehensive coverage: Every module, subworkflow, and pipeline should have tests
2. Efficient test data: Tests use small, centralized datasets from nf-core/test-datasets
3. Snapshot testing: Output verification through snapshots ensures consistency
4. CI/CD ready: Tests run reliably in automated environments
5. Standard patterns: Consistent testing patterns across all nf-core components

## Getting help

For issues not covered in this guide:

- Browse examples: nf-core/methylseq, nf-core/sarek tests for current best practices
- nf-test documentation: [Official nf-test docs](https://code.askimed.com/nf-test/)
- Community support: [nf-core Slack](https://nf-co.re/join) `#nf-test` channel
- nf-core modules: [GitHub repository](https://github.com/nf-core/modules) for examples

Ready to get started? Begin with [Installation](./components/01_installation.md) to set up nf-test in your development environment.
