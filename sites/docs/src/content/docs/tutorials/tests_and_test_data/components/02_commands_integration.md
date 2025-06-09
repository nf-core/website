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

## nf-core Integration

### Using nf-core tools with nf-test. nf-core tools only provides test commands for the `nf-core/modules` repository context.

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

#### Testing modules

```bash
nf-core modules test [OPTIONS] <tool> or <tool/subtool>
```

Example help output:
```
Usage: nf-core modules test [OPTIONS] <tool> or <tool/subtool>

Run nf-test for a module.

╭─ Options ─────────────────────────────────────────────────────────────────╮
│ --verbose         -v    Print verbose output to the console. Sets         │
│                         `--debug` inside the nf-test command.             │
│ --dir             -d    <nf-core/modules directory>                       │
│ --no-prompts      -p    Use defaults without prompting                    │
│ --update          -u    Update existing snapshots                         │
│ --once            -o    Run tests only once. Don't check snapshot         │
│                         stability                                         │
│ --profile               [docker|singularity|conda] Run tests with a       │
│                         specific profile                                  │
│ --migrate-pytest        Migrate a module with pytest tests to nf-test     │
│ --help            -h    Show this message and exit.                       │
╰───────────────────────────────────────────────────────────────────────────╯
```

#### Testing subworkflows

```bash
nf-core subworkflows test [OPTIONS] subworkflow name
```

Example help output:
```
Usage: nf-core subworkflows test [OPTIONS] subworkflow name

Run nf-test for a subworkflow.

╭─ Options ────────────────────────────────────────────────────────────────╮
│ --dir             -d    <nf-core/modules directory>                      │
│ --no-prompts      -p    Use defaults without prompting                   │
│ --update          -u    Update existing snapshots                        │
│ --once            -o    Run tests only once. Don't check snapshot        │
│                         stability                                        │
│ --profile               [docker|singularity|conda] Run tests with a      │
│                         specific profile                                 │
│ --migrate-pytest        Migrate a subworkflow with pytest tests to       │
│                         nf-test                                          │
│ --help            -h    Show this message and exit.                      │
╰───────────────────────────────────────────────────────────────────────────╯
```


When working within the **nf-core/modules repository**, you can use nf-core tools commands:

```bash
# Test a specific module (only works in nf-core/modules repo)
nf-core modules test bedtools/bamtobed --profile docker

# Test with docker profile
nf-core modules test bedtools/bamtobed --profile docker

# Test a subworkflow (only works in nf-core/modules repo)
nf-core subworkflows test vcf_impute_glimpse --profile docker
```

### In individual pipeline repositories

When developing **individual nf-core pipelines**, use nf-test commands directly:

```bash
# List all available tests
nf-test list

# Run specific tests with profiles
nf-test test tests/default.nf-test --profile test,docker

# Run all tests
nf-test test --profile docker,test

# Update snapshots
nf-test test tests/default.nf-test --profile docker,test --update-snapshot
```


## Test Data Integration


nf-core tools 3.3+ includes a new command for discovering and managing test datasets from the [nf-core/test-datasets repository](https://github.com/nf-core/test-datasets):

```bash
# List all available test datasets
nf-core test-datasets list

# List datasets for a specific branch
nf-core test-datasets list --branch mag

# Search for specific datasets
nf-core test-datasets search --branch mag minigut_reads

# Generate download URLs for datasets
nf-core test-datasets list --branch mag --generate-dl-url

# Generate nf-test compatible paths
nf-core test-datasets list --branch mag --generate-nf-path
```

This command provides streamlined functionality for:
- **Dataset discovery**: Easily search and explore available test datasets
- **Dataset integration**: Get download URLs or nf-test compatible paths for use in your pipeline tests

> **Note:** This feature was introduced in [nf-core/tools 3.3](https://nf-co.re/blog/2025/tools-3_3#new-nf-core-test-datasets-command). Make sure you have the latest version of nf-core tools installed.


## Essential nf-test Commands

### Running Tests

```bash
# List all available tests
nf-test list

# Run all tests in a repository
nf-test test --profile docker
```

```bash
# Run specific test file in nf-core/modules repository
nf-test test modules/nf-core/fastqc/tests/main.nf.test --profile docker

# Run nf-core/modules tests with specific tag
nf-test test --tag bedtools/bamtobed --profile docker

# Run nf-core/modules tests in verbose mode
nf-test test --tag bedtools/bamtobed --profile docker --verbose
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
nf-test test bedtools/bamtobed --profile docker --update-snapshot

# Clean the snapshot for any invalid content
nf-test test --tag bedtools/bamtobed --clean-snapshot --profile docker
```

## Command Options

### Core Options

| Option | Description |
|--------|-------------|
| `--profile` | Specify Nextflow configuration profile(s) to use |
| `--tag` | Run tests with specific tag(s) |
| `--verbose` | Enable verbose output for debugging |
| `--update-snapshot` | Update test snapshots with new expected values |
| `--clean-snapshot` | Clean and regenerate snapshots |

### Test Selection

| Option | Description |
|--------|-------------|
| `--filter` | Run tests matching a specific pattern |
| `--dry-run` | Show which tests would run without executing them |
| `--changed` | Run only tests affected by recent changes |

### Execution Control

| Option | Description |
|--------|-------------|
| `--silent` | Suppress all output except errors |
| `--ci` | Enable CI mode with optimized settings |
| `--shard` | Split tests across multiple shards for parallel execution |
| `--index` | Specify shard index when using `--shard` |

### Configuration

| Option | Description |
|--------|-------------|
| `--config` | Specify custom nf-test configuration file |
| `--global` | Run tests in global mode |
| `--follow-dependencies` | Automatically run dependent tests |

### Output Control

| Option | Description |
|--------|-------------|
| `--tap` | Output results in TAP (Test Anything Protocol) format |
| `--junitxml` | Generate JUnit XML test report |


## Next Steps

Continue to [Project Setup](./03_project_setup.md) to configure your testing environment.