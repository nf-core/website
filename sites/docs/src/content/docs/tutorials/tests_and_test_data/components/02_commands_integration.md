---
title: "2. nf-test Commands & Integration"
subtitle: Essential commands and nf-core integration
weight: 20
---

## Essential nf-test Commands

### Running Tests

```bash
# Run all tests
nf-test test

# Run specific test file
nf-test test modules/nf-core/fastqc/tests/main.nf.test

# Run tests with specific tag
nf-test test --tag fastqc

# Run tests in verbose mode
nf-test test --verbose
```

### Generating Tests

```bash
# Generate test template
nf-test generate process path/to/main.nf

# Generate test for workflow
nf-test generate workflow path/to/main.nf
```

### Managing Snapshots

```bash
# Update snapshots
nf-test test --update-snapshot

# Clean test directory
nf-test clean
```

## nf-core Integration

### Using nf-core tools with nf-test

```bash
# Create module with test template
nf-core modules create <tool>

# Test specific module
nf-core modules test <tool>

# Lint module tests
nf-core modules lint <tool>

# Search for available modules
nf-core modules list

# Get information about a module
nf-core modules info <tool>
```

### Test Data Integration

nf-core uses centralized test data:

```groovy
// Reference test data
file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/genome.fasta', checkIfExists: true)
```

## Prerequisites

Before running these commands, ensure you have:
- nf-test installed (see [Installation Guide](./01_installation.md))
- An nf-core pipeline or module project set up
- Access to test data (local or remote)

## Command Options

### Common Flags

- `--profile` - Specify Nextflow profile
- `--tag` - Run tests with specific tags
- `--verbose` - Detailed output
- `--silent` - Minimal output
- `--update-snapshot` - Update test snapshots

### Example Usage

```bash
# Run module tests with conda profile
nf-test test --tag modules --profile conda

# Update snapshots for specific module
nf-test test modules/nf-core/fastqc/tests/main.nf.test --update-snapshot
```

## Next Steps

Continue to [Project Setup](./03_project_setup.md) to configure your testing environment. 