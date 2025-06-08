---
title: "nf-test: Comprehensive Testing Guide"
subtitle: Complete guide to testing nf-core components, subworkflows, and pipelines with nf-test
shortTitle: nf-test Comprehensive Guide
weight: 5
---

## Overview

This comprehensive guide covers testing nf-core components using nf-test. From writing your first test to implementing advanced patterns, this guide provides everything needed to create robust, maintainable tests for modules, subworkflows, and pipelines.

## Guide Structure

### Getting Started
- **[1. Installation](./components/01_installation.md)** - Setting up nf-test in your development environment
- **[2. nf-test Commands & Integration](./components/02_commands_integration.md)** - Essential commands and nf-core integration
- **[3. Project Setup](./components/03_project_setup.md)** - Configuring your project for testing

### Component Testing
- **[4. Testing Modules](./components/04_testing_modules.md)** - Testing individual nf-core modules
- **[5. Testing Subworkflows](./components/05_testing_subworkflows.md)** - Testing subworkflow components
- **[6. Testing Pipelines](./components/06_testing_pipelines.md)** - End-to-end pipeline testing
- **[7. nf-test Assertions](./components/07_assertions.md)** - Comprehensive assertion patterns and verification techniques

### Configuration & Management
- **[8. Configuration Management](./components/07_configuration_management.md)** - Managing test configurations and parameters
- **[9. Snapshot Management](./components/08_snapshot_management.md)** - Working with snapshots and test outputs

### Advanced Topics
- **[10. Advanced Testing Patterns](./components/09_advanced_testing_patterns.md)** - Complex testing scenarios and patterns
- **[11. Test Data Management](./components/10_test_data_management.md)** - Organizing and managing test datasets
- **[12. External Commands and Tools](./components/11_external_commands_tools.md)** - Testing with external dependencies
- **[13. Custom Utility Classes](./components/12_custom_utility_classes.md)** - Creating reusable test utilities

### Quality & Deployment
- **[14. Error Testing & Edge Cases](./components/13_error_testing_edge_cases.md)** - Testing failure scenarios and edge cases
- **[15. CI/CD Integration](./components/14_cicd_integration.md)** - Integrating tests with continuous integration
- **[16. Best Practices](./components/15_best_practices.md)** - Testing best practices and conventions
- **[17. Troubleshooting](./components/16_troubleshooting.md)** - Common issues and solutions

## nf-core Testing Principles

nf-core testing follows these key principles:

1. **Comprehensive Coverage**: Every module, subworkflow, and pipeline should have tests
2. **Realistic Data**: Tests use representative data from standardized repositories
3. **Snapshot Testing**: Output verification through snapshots ensures consistency
4. **CI/CD Ready**: Tests run reliably in automated environments
5. **Standard Patterns**: Consistent testing patterns across all nf-core components

Current best practices are exemplified in nf-core/methylseq, demonstrating essential pipeline testing patterns, proper plugin usage, and comprehensive `.nftignore` configurations.

## Getting Help

For issues not covered in this guide:

- **Browse examples**: nf-core/methylseq tests for current best practices
- **nf-test documentation**: [Official nf-test docs](https://code.askimed.com/nf-test/)
- **Community support**: [nf-core Slack](https://nf-co.re/join) `#help` channel
- **nf-core modules**: [GitHub repository](https://github.com/nf-core/modules) for examples

Ready to get started? Begin with [Installation](./components/01_installation.md) to set up nf-test in your development environment. 