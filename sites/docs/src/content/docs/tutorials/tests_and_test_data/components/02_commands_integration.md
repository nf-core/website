---
title: "2. nf-test Commands & Integration"
subtitle: Essential commands and nf-core integration
weight: 20
---

## Prerequisites

Before running these commands, ensure you have:
- nf-test installed (see [Installation Guide](./01_installation.md))
- An nf-core pipeline or nf-core/modules repo set up
- Access to test data (local or remote)

## Understanding Testing Contexts

There are **two main contexts** for running nf-test commands in the nf-core ecosystem:

### 1. nf-core/modules Repository Context
When working in the **nf-core/modules repository**, use `nf-core` tools commands:
- Testing individual modules and subworkflows
- Contributing to the shared nf-core modules library
- Use `nf-core modules test` and `nf-core subworkflows test`

### 2. Individual Pipeline Context
When developing **individual nf-core pipelines**, use `nf-test` commands directly:
- Testing complete pipeline workflows
- Pipeline-specific test configurations
- Use `nf-test test` commands

---

## nf-core/modules Repository Commands

> **Note**: These commands only work within the nf-core/modules repository

### Module Testing

```bash
# Test a specific module
nf-core modules test bedtools/bamtobed --profile docker

# Update module test snapshots
nf-core modules test bedtools/bamtobed --profile docker --update

# Test with verbose output
nf-core modules test bedtools/bamtobed --profile docker --verbose
```

### Subworkflow Testing

```bash
# Test a subworkflow
nf-core subworkflows test vcf_impute_glimpse --profile docker

# Update subworkflow snapshots
nf-core subworkflows test vcf_impute_glimpse --profile docker --update
```

### Module Management

```bash
# Create module with test template
nf-core modules create <tool>

# Lint module tests
nf-core modules lint <tool>

# Search for available modules
nf-core modules list

# Get information about a module
nf-core modules info <tool>
```

---

## Individual Pipeline Commands

> **Note**: Use these commands when working on individual nf-core pipelines

### Basic Testing

```bash
# List all available tests
nf-test list

# Run all tests with profiles
nf-test test --profile test,docker

# Run specific test file
nf-test test tests/default.nf-test --profile test,docker
```

### Test Management

```bash
# Update snapshots
nf-test test --profile test,docker --update-snapshot

# Run tests with verbose output
nf-test test --profile test,docker --verbose

# Run specific test by tag
nf-test test --tag pipeline --profile test,docker
```

### Test Selection

```bash
# Run tests matching a pattern
nf-test test --filter "*fastqc*" --profile docker

# Run only changed tests
nf-test test --changed --profile docker

# Dry run to see which tests would execute
nf-test test --dry-run --profile docker
```

---

## Next Steps

Continue to [Project Setup](./03_project_setup.md) to configure your testing environment.
