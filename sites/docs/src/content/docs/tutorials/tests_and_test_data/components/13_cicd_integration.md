---
title: "15. CI/CD Integration"
subtitle: Integrating nf-test with continuous integration
weight: 150
---

## Real nf-core CI/CD Setup

This section provides production-ready examples of CI/CD integration with nf-test, featuring advanced sharding, multiple profiles, and GPU testing.

## Main nf-test Workflow

### Complete GitHub Actions Workflow

Here's a standard nf-core pipeline nf-test workflow:

```yaml
name: Run nf-tests
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
  release:
    types: [published]
  workflow_dispatch:
    inputs:
      runners:
        description: "Runners to test on"
        type: choice
        options:
          - "ubuntu-latest"
          - "self-hosted"
        default: "ubuntu-latest"

# Cancel if a newer run is started
concurrency:
  group: ${{ github.workflow }}-${{ github.event.pull_request.number || github.ref }}
  cancel-in-progress: true

env:
  GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}
  # renovate: datasource=github-releases depName=askimed/nf-test versioning=semver
  NFT_VER: "0.9.2"
  NXF_ANSI_LOG: false
  NXF_SINGULARITY_CACHEDIR: ${{ github.workspace }}/.singularity
  NXF_SINGULARITY_LIBRARYDIR: ${{ github.workspace }}/.singularity
  # renovate: datasource=github-releases depName=nextflow-io/nextflow versioning=semver
  NXF_VER: "24.10.2"

jobs:
  get-shards:
    runs-on:
      - runs-on=${{ github.run_id }}-nf-test-get-shards
      - runner=2cpu-linux-x64
    name: "Get Shards"
    outputs:
      shard: ${{ steps.set-shards.outputs.shard }}
      total_shards: ${{ steps.set-shards.outputs.total_shards }}
    steps:
      - name: Check out pipeline code
        uses: actions/checkout@0ad4b8fadaa221de15dcec353f45205ec38ea70b # v4
        with:
          fetch-depth: 0

      - name: run nf-test-shard action
        id: set-shards
        uses: ./.github/actions/nf-test-shard
        env:
          NFT_VER: ${{ env.NFT_VER }}
        with:
          tags: "cpu"
          max_shards: 15

  nf-test:
    needs: [get-shards]
    runs-on:
      - runs-on=${{ github.run_id }}
      - runner=4cpu-linux-x64
      - image=ubuntu22-full-x64
    name: "Test | ${{ matrix.profile }} | ${{ matrix.shard }} | ${{ matrix.NXF_VER }} | ${{ matrix.filters }}"
    env:
      NXF_VER: ${{ matrix.NXF_VER }}

    # Only run on push if this is the nf-core dev branch (merged PRs)
    # and number of shards is greater than 0
    if: |
      needs.get-shards.outputs.total_shards > 0 &&
      (
        github.event_name != 'push' ||
        ( github.event_name == 'push' && github.repository == 'nf-core/methylseq' )
      )

    strategy:
      fail-fast: false
      matrix:
        profile: [conda, docker, singularity]
        shard: ${{ fromJson(needs.get-shards.outputs.shard) }}
        NXF_VER:
          - "24.10.2"
        filters: [pipeline]

    steps:
      - name: Check out pipeline code
        uses: actions/checkout@0ad4b8fadaa221de15dcec353f45205ec38ea70b # v4
        with:
          fetch-depth: 0

      - name: Run nf-test Action
        uses: ./.github/actions/nf-test
        with:
          profile: ${{ matrix.profile }}
          shard: ${{ matrix.shard }}
          total_shards: ${{ needs.get-shards.outputs.total_shards }}
          filters: ${{ matrix.filters }}
          tags: "cpu"

  confirm-pass:
    runs-on:
      - runs-on=${{ github.run_id }}-confirm-pass
      - runner=2cpu-linux-x64
    needs: [nf-test]
    if: always()
    steps:
      - name: One or more tests failed
        if: ${{ contains(needs.*.result, 'failure') }}
        run: exit 1

      - name: One or more tests cancelled
        if: ${{ contains(needs.*.result, 'cancelled') }}
        run: exit 1

      - name: All tests ok
        if: ${{ contains(needs.*.result, 'success') }}
        run: exit 0
```

## Advanced CI/CD Features

### 1. Test Sharding

The methylseq workflow uses advanced sharding to parallelize tests:

#### Dynamic Shard Generation (`.github/actions/nf-test-shard/action.yml`)

```yaml
name: "Get number of shards"
description: "Get the number of nf-test shards for the current CI job"
inputs:
  tags:
    description: "Tags to pass as argument for nf-test --tag parameter"
    required: true
  max_shards:
    description: "Maximum number of shards allowed"
    required: true
outputs:
  shard:
    description: "Array of shard numbers"
    value: ${{ steps.shards.outputs.shard }}
  total_shards:
    description: "Total number of shards"
    value: ${{ steps.shards.outputs.total_shards }}

runs:
  using: "composite"
  steps:
    - name: Install nf-test
      uses: nf-core/setup-nf-test@v1
      with:
        version: ${{ env.NFT_VER }}
        install-pdiff: true

    - name: Get number of shards
      id: shards
      shell: bash
      run: |
        # Run nf-test to get the number of related tests
        nftest_output=$(nf-test test --dry-run --changed-since HEAD^ --filter pipeline --tag ${{ inputs.tags }})
        echo "nf-test dry-run output: $nftest_output"

        # Default values for shard and total_shards
        shard="[]"
        total_shards=0

        # Check if there are related tests
        if echo "$nftest_output" | grep -q 'Nothing to do'; then
          echo "No related tests found."
        else
          # Extract the number of related tests
          number_of_shards=$(echo "$nftest_output" | sed -n 's|.*Executed \([0-9]*\) tests.*|\1|p')
          if [[ -n "$number_of_shards" && "$number_of_shards" -gt 0 ]]; then
            shards_to_run=$(( $number_of_shards < ${{ inputs.max_shards }} ? $number_of_shards : ${{ inputs.max_shards }} ))
            shard=$(seq 1 "$shards_to_run" | jq -R . | jq -c -s .)
            total_shards="$shards_to_run"
          else
            echo "Unexpected output format. Falling back to default values."
          fi
        fi

        # Write to GitHub Actions outputs
        echo "shard=$shard" >> $GITHUB_OUTPUT
        echo "total_shards=$total_shards" >> $GITHUB_OUTPUT

        # Debugging output
        echo "Final shard array: $shard"
        echo "Total number of shards: $total_shards"
```

#### Key Sharding Features:

- **Dynamic discovery**: Uses `nf-test --dry-run` to discover tests
- **Intelligent sharding**: Only creates shards for changed tests
- **Resource optimization**: Limits maximum shards to prevent resource waste
- **Tag-based filtering**: Supports CPU vs GPU test separation

### 2. Comprehensive nf-test Action

The main test action (`.github/actions/nf-test/action.yml`) handles all test execution:

```yaml
name: "nf-test Action"
description: "Runs nf-test with common setup steps"
inputs:
  profile:
    description: "Profile to use"
    required: true
  shard:
    description: "Shard number for this CI job"
    required: true
  total_shards:
    description: "Total number of test shards(NOT the total number of matrix jobs)"
    required: true
  filters:
    description: "Filter test cases by specified types (e.g., module, pipeline, workflow or function)"
    required: true
  tags:
    description: "Tags to pass as argument for nf-test --tag parameter"
    required: true

runs:
  using: "composite"
  steps:
    - uses: actions/setup-java@8df1039502a15bceb9433410b1a100fbe190c53b # v4
      with:
        distribution: "temurin"
        java-version: "17"

    - name: Setup Nextflow
      uses: nf-core/setup-nextflow@v2
      with:
        version: "${{ env.NXF_VER }}"

    - name: Set up Python
      uses: actions/setup-python@0b93645e9fea7318ecaed2b359559ac225c90a2b # v5
      with:
        python-version: "3.11"

    - name: Install nf-test
      uses: nf-core/setup-nf-test@v1
      with:
        version: "${{ env.NFT_VER }}"
        install-pdiff: true

    - name: Setup apptainer
      if: contains(inputs.profile, 'singularity')
      uses: eWaterCycle/setup-apptainer@main

    - name: Set up Singularity
      if: contains(inputs.profile, 'singularity')
      shell: bash
      run: |
        mkdir -p $NXF_SINGULARITY_CACHEDIR
        mkdir -p $NXF_SINGULARITY_LIBRARYDIR

    - name: Conda setup
      if: ${{inputs.profile == 'conda'}}
      uses: conda-incubator/setup-miniconda@d2e6a045a86077fb6cad6f5adf368e9076ddaa8d # v3
      with:
        auto-update-conda: true
        conda-solver: libmamba
        conda-remove-defaults: true

    - name: Run nf-test
      shell: bash
      run: |
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

    - name: Generate test summary
      if: always()
      shell: bash
      run: |
        # Add header if it doesn't exist (using a token file to track this)
        if [ ! -f ".summary_header" ]; then
          echo "# 🚀 nf-test results" >> $GITHUB_STEP_SUMMARY
          echo "" >> $GITHUB_STEP_SUMMARY
          echo "| Status | Test Name | Profile | Shard |" >> $GITHUB_STEP_SUMMARY
          echo "|:------:|-----------|---------|-------|" >> $GITHUB_STEP_SUMMARY
          touch .summary_header
        fi

        if [ -f test.tap ]; then
          while IFS= read -r line; do
            if [[ $line =~ ^ok ]]; then
              test_name="${line#ok }"
              # Remove the test number from the beginning
              test_name="${test_name#* }"
              echo "| ✅ | ${test_name} | ${{ inputs.profile }} | ${{ inputs.shard }}/${{ inputs.total_shards }} |" >> $GITHUB_STEP_SUMMARY
            elif [[ $line =~ ^not\ ok ]]; then
              test_name="${line#not ok }"
              # Remove the test number from the beginning
              test_name="${test_name#* }"
              echo "| ❌ | ${test_name} | ${{ inputs.profile }} | ${{ inputs.shard }}/${{ inputs.total_shards }} |" >> $GITHUB_STEP_SUMMARY
            fi
          done < test.tap
        else
          echo "| ⚠️ | No test results found | ${{ inputs.profile }} | ${{ inputs.shard }}/${{ inputs.total_shards }} |" >> $GITHUB_STEP_SUMMARY
        fi

    - name: Clean up
      if: always()
      shell: bash
      run: |
        sudo rm -rf /home/ubuntu/tests/
```

## GPU Testing Workflow

methylseq also includes specialized GPU testing:

### GPU nf-test Workflow (`.github/workflows/nf-test-gpu.yml`)

```yaml
name: Run GPU nf-tests
on:
  push:
    branches:
      - dev
  pull_request:
  release:
    types: [published]
  workflow_dispatch:
    inputs:
      runners:
        description: "Runners to test on"
        type: string
        default: "gpu"

# Cancel if a newer run is started
concurrency:
  group: ${{ github.workflow }}-${{ github.event.pull_request.number || github.ref }}
  cancel-in-progress: true

env:
  GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}
  NFT_VER: "0.9.2"
  NXF_ANSI_LOG: false
  NXF_SINGULARITY_CACHEDIR: ${{ github.workspace }}/.singularity
  NXF_SINGULARITY_LIBRARYDIR: ${{ github.workspace }}/.singularity
  NXF_VER: "24.10.2"

jobs:
  get-shards:
    runs-on:
      - runs-on=${{ github.run_id }}-nf-test-get-shards-gpu
      - runner=2cpu-linux-x64
    name: "Get Shards"
    outputs:
      shard: ${{ steps.set-shards.outputs.shard }}
      total_shards: ${{ steps.set-shards.outputs.total_shards }}
    steps:
      - name: Check out pipeline code
        uses: actions/checkout@0ad4b8fadaa221de15dcec353f45205ec38ea70b # v4
        with:
          fetch-depth: 0

      - name: run nf-test-shard action
        id: set-shards
        uses: ./.github/actions/nf-test-shard
        env:
          NFT_VER: ${{ env.NFT_VER }}
        with:
          tags: "gpu"
          max_shards: 2

  nf-test-gpu:
    needs: [get-shards]
    runs-on: "runs-on=${{ github.run_id }}/family=g4dn.xlarge/image=ubuntu24-gpu-x64"
    name: "GPU Test | ${{ matrix.profile }} | ${{ matrix.shard }} | ${{ matrix.NXF_VER }} | ${{ matrix.filters }}"
    env:
      NXF_VER: ${{ matrix.NXF_VER }}

    if: |
      needs.get-shards.outputs.total_shards > 0 &&
      (
        github.event_name != 'push' ||
        ( github.event_name == 'push' && github.repository == 'nf-core/methylseq' )
      )
    strategy:
      fail-fast: false
      matrix:
        profile: [docker, singularity]
        shard: ${{ fromJson(needs.get-shards.outputs.shard) }}
        NXF_VER:
          - "24.10.2"
        filters: [pipeline]

    steps:
      - name: Check out pipeline code
        uses: actions/checkout@0ad4b8fadaa221de15dcec353f45205ec38ea70b # v4
        with:
          fetch-depth: 0

      - name: Test CUDA
        run: |
          nvidia-smi -L

      - name: Run nf-test Action
        uses: ./.github/actions/nf-test
        with:
          profile: ${{ matrix.profile }},gpu
          shard: ${{ matrix.shard }}
          total_shards: ${{ needs.get-shards.outputs.total_shards }}
          filters: ${{ matrix.filters }}
          tags: "gpu"
```

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

## Implementing nf-core CI/CD in Your Project

### 1. Copy Custom Actions

Create `.github/actions/` directory with:

```
.github/actions/
├── nf-test/
│   └── action.yml           # Main test execution action
├── nf-test-shard/
│   └── action.yml           # Dynamic sharding action
```

### 2. Set Up Workflows

#### Basic nf-test Workflow

```yaml
# .github/workflows/nf-test.yml
name: Run nf-tests
on:
  push:
    branches: [dev, main]
  pull_request:

concurrency:
  group: ${{ github.workflow }}-${{ github.event.pull_request.number || github.ref }}
  cancel-in-progress: true

env:
  NFT_VER: "0.9.2"
  NXF_VER: "24.10.2"

jobs:
  nf-test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        profile: [docker, singularity]
    steps:
      - uses: actions/checkout@v4
      - uses: ./.github/actions/nf-test
        with:
          profile: ${{ matrix.profile }}
          shard: "1"
          total_shards: "1"
          filters: "pipeline"
          tags: "pipeline"
```

### 3. Configure Test Reporting

Add TAP output processing for detailed test reports:

```yaml
- name: Generate test summary
  if: always()
  shell: bash
  run: |
    if [ -f test.tap ]; then
      echo "# Test Results" >> $GITHUB_STEP_SUMMARY
      echo "| Status | Test Name |" >> $GITHUB_STEP_SUMMARY
      echo "|:------:|-----------|" >> $GITHUB_STEP_SUMMARY
      
      while IFS= read -r line; do
        if [[ $line =~ ^ok ]]; then
          test_name="${line#ok }"
          test_name="${test_name#* }"
          echo "| ✅ | ${test_name} |" >> $GITHUB_STEP_SUMMARY
        elif [[ $line =~ ^not\ ok ]]; then
          test_name="${line#not ok }"
          test_name="${test_name#* }"
          echo "| ❌ | ${test_name} |" >> $GITHUB_STEP_SUMMARY
        fi
      done < test.tap
    fi
```

## Best Practices from methylseq

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

Continue to [Best Practices](./14_best_practices.md) to learn comprehensive testing strategies and conventions.
