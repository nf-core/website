---
title: "9. CI/CD Integration"
subtitle: Integrating nf-test with continuous integration
weight: 90
---

## nf-core CI/CD Setup

This section provides production-ready examples of CI/CD integration with nf-test, featuring advanced sharding, multiple profiles, and GPU testing.

## Main nf-test Workflow

### Complete GitHub Actions Workflow

## GitHub Actions Workflow: Run nf-tests

| **Component**             | **Details**                    | **Values/Configuration**                                                                |
| ------------------------- | ------------------------------ | --------------------------------------------------------------------------------------- |
| **Workflow Name**         | Run nf-tests                   | -                                                                                       |
| **Triggers**              | Events that start the workflow | `push` (dev branch), `pull_request`, `release` (published), `workflow_dispatch`         |
| **Path Exclusions**       | Files ignored for triggers     | docs/\*_, meta.yml, _.md, _.png, _.svg                                                  |
| **Concurrency**           | Prevents multiple runs         | Group: `${{ github.workflow }}-${{ github.event.pull_request.number \|\| github.ref }}` |
| **Environment Variables** | Global settings                | `GITHUB_TOKEN`, `NFT_VER: 0.9.2`, `NXF_VER: 24.10.2`, `NXF_ANSI_LOG: false`             |

### Jobs Breakdown

| **Job Name**     | **Purpose**             | **Runner**                                        | **Key Outputs/Actions**              |
| ---------------- | ----------------------- | ------------------------------------------------- | ------------------------------------ |
| **get-shards**   | Determine test sharding | `runs-on=${{ github.run_id }}-nf-test-get-shards` | `shard`, `total_shards` outputs      |
| **nf-test**      | Execute tests           | `runs-on=${{ github.run_id }}` (4cpu-linux-x64)   | Run tests across matrix combinations |
| **confirm-pass** | Validate results        | `runs-on=${{ github.run_id }}-confirm-pass`       | Check if all tests passed            |

### Strategy Matrix (nf-test job)

| **Matrix Dimension** | **Options**                 | **Purpose**                             |
| -------------------- | --------------------------- | --------------------------------------- |
| **profile**          | conda, docker, singularity  | Test different containerization methods |
| **shard**            | Dynamic from get-shards job | Parallel test execution                 |
| **NXF_VER**          | 24.10.2                     | Nextflow version                        |
| **filters**          | pipeline                    | Test scope                              |

### Key Steps per Job

| **Job**          | **Step**          | **Action**                                                  | **Purpose**                 |
| ---------------- | ----------------- | ----------------------------------------------------------- | --------------------------- |
| **get-shards**   | Check out code    | `actions/checkout@0ad4b8fadaa221de15dcec353f45205ec38ea70b` | Get pipeline code           |
|                  | Run nf-test-shard | `./.github/actions/nf-test-shard`                           | Calculate test shards       |
| **nf-test**      | Check out code    | `actions/checkout@0ad4b8fadaa221de15dcec353f45205ec38ea70b` | Get pipeline code           |
|                  | Run nf-test       | `./.github/actions/nf-test`                                 | Execute tests with sharding |
| **confirm-pass** | Validate results  | Custom shell commands                                       | Check overall test status   |

### Conditional Logic

| **Condition**        | **Applies To**   | **Logic**                                                                   |
| -------------------- | ---------------- | --------------------------------------------------------------------------- |
| **Job Execution**    | nf-test job      | Only runs if `total_shards > 0` AND (not push OR push to nf-core/methylseq) |
| **Failure Handling** | confirm-pass job | Exits with error if any test failed or was cancelled                        |

## Key CI/CD Features

### 1. **Smart Triggers**

```yaml
on:
  push:
    paths-ignore:
      - "docs/**"
      - "**/meta.yml"
      - "**/*.md"
      - "**/*.png"
      - "**/*.svg"
    branches:
      - dev
  pull_request:
    paths-ignore:
      - "docs/**"
      - "**/meta.yml"
      - "**/*.md"
      - "**/*.png"
      - "**/*.svg"
```

- **Path-based filtering**: Only runs when relevant files change
- **Branch protection**: Only runs on specific branches
- **Documentation exclusion**: Skips runs for documentation changes

### 2. **Resource Optimization**

```yaml
# Cancel if a newer run is started
concurrency:
  group: ${{ github.workflow }}-${{ github.event.pull_request.number || github.ref }}
  cancel-in-progress: true
```

- **Concurrency control**: Cancels redundant runs
- **Resource limits**: Uses appropriate runner sizes
- **Intelligent sharding**: Only creates necessary shards

### 3. **Multi-Profile Testing**

```yaml
strategy:
  fail-fast: false
  matrix:
    profile: [conda, docker, singularity]
    shard: ${{ fromJson(needs.get-shards.outputs.shard) }}
    NXF_VER:
      - "24.10.2"
    filters: [pipeline]
```

- **Container engines**: Tests Docker, Singularity, and Conda
- **Version testing**: Tests specific Nextflow versions
- **Fail-fast disabled**: Allows all profiles to complete

### 4. **Advanced nf-test Options**

```bash
NFT_WORKDIR=~ \
nf-test test \
  --ci \
  --shard ${{ inputs.shard }}/${{ inputs.total_shards }} \
  --changed-since HEAD^ \
  --profile=${{ inputs.profile }} \
  --filter ${{ inputs.filters }} \
  --tap=test.tap \
  --verbose \
  --tag ${{ inputs.tags }}
```

Key options:

- **`--ci`**: CI-optimized mode
- **`--shard`**: Parallel test execution
- **`--changed-since HEAD^`**: Only test changed files
- **`--tap=test.tap`**: TAP output for reporting
- **`--tag`**: Tag-based test filtering

## Environment Configuration

### Essential Environment Variables

```yaml
env:
  GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}
  NFT_VER: "0.9.2" # nf-test version
  NXF_ANSI_LOG: false # Disable ANSI logs for CI
  NXF_SINGULARITY_CACHEDIR: ${{ github.workspace }}/.singularity
  NXF_SINGULARITY_LIBRARYDIR: ${{ github.workspace }}/.singularity
  NXF_VER: "24.10.2" # Nextflow version
```

### Version Management

methylseq uses Renovate for automated version updates:

```yaml
# renovate: datasource=github-releases depName=askimed/nf-test versioning=semver
NFT_VER: "0.9.2"
# renovate: datasource=github-releases depName=nextflow-io/nextflow versioning=semver
NXF_VER: "24.10.2"
```

## Recommended CI/CD Patterns

### 1. **Progressive Testing**

- Use `--changed-since HEAD^` to test only modified components
- Implement sharding for large test suites
- Separate CPU and GPU tests

### 2. **Resource Management**

- Set appropriate `NFT_WORKDIR` for CI environments
- Use concurrency controls to prevent resource conflicts
- Clean up test artifacts after completion

### 3. **Comprehensive Coverage**

- Test multiple container engines (Docker, Singularity, Conda)
- Test different Nextflow versions
- Include both unit and integration tests

### 4. **Monitoring and Reporting**

- Generate detailed test summaries
- Use TAP output for structured reporting
- Include debugging information for failed tests

### 5. **Version Management**

- Pin specific versions of tools for reproducibility
- Use automated dependency updates (Renovate)
- Test against multiple Nextflow versions

## Next Steps

Continue to [FAQ & Debugging](./10_faq_debugging.md) to learn comprehensive testing strategies, troubleshooting, and best practices.
