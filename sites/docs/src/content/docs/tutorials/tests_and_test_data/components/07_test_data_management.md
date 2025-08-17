---
title: "7. Test Data Management"
subtitle: Organizing and managing test datasets
weight: 70
---

## Test Data Strategy

Understanding how test data flows through the nf-core ecosystem is crucial for effective testing. The following diagram illustrates the complete test data management strategy:

```mermaid
flowchart TD
    A[nf-core Test Data Strategy] --> B[Centralized Repositories]

    B --> B1[nf-core/test-datasets]
    B1 --> B2[modules branch<br/>Module test data]
    B1 --> B3[pipeline branches<br/>Pipeline-specific data]

    B2 --> C[Module Test Data Access]
    B3 --> D[Pipeline Test Data Access]

    C --> C1[modules_testdata_base_path]
    C1 --> C2[genomics/sarscov2/illumina/fastq]
    C1 --> C3[genomics/homo_sapiens/genome]
    C1 --> C4[proteomics/database]

    D --> D1[pipelines_testdata_base_path]
    D1 --> D2[Pipeline-specific datasets]
    D2 --> D3[Sample sheets]
    D2 --> D4[Configuration files]
    D2 --> D5[Reference data]

    A --> E[Test Data Discovery]
    E --> E1[nf-core test-datasets command]
    E1 --> E2[list branches]
    E1 --> E3[search datasets]
    E1 --> E4[generate URLs]
    E1 --> E5[generate nf-paths]

    E2 --> F[Integration Methods]
    E3 --> F
    E4 --> F
    E5 --> F

    F --> F1[Direct URL References]
    F --> F2[Parameter-based Paths]
    F --> F3[Local Asset Files]

    F1 --> G[Implementation in Tests]
    F2 --> G
    F3 --> G

    G --> G1[checkIfExists: true]
    G --> G2[file function usage]
    G --> G3[Channel creation]

    G1 --> H[Best Practices]
    G2 --> H
    G3 --> H

    H --> H1[Use consistent base paths]
    H --> H2[Verify file existence]
    H --> H3[Organize by data type]
    H --> H4[Document data requirements]

    H1 --> I[Efficient Testing]
    H2 --> I
    H3 --> I
    H4 --> I
```

### nf-core Test Data Repository

nf-core maintains centralized test datasets:

- **Modules test data**: `https://github.com/nf-core/test-datasets/tree/modules`
- **Pipeline test data**: `https://github.com/nf-core/test-datasets/tree/<pipeline-name>`

```groovy
// Standard nf-core test data reference
file(params.modules_testdata_base_path + 'genomics/sarscov2/illumina/fastq/test_1.fastq.gz', checkIfExists: true)
```

### New nf-core test-datasets Command

> **New in nf-core/tools 3.3**: The `nf-core test-datasets` command provides streamlined functionality for dataset discovery and integration.

#### Available Commands

```bash
# List available subcommands
nf-core test-datasets --help

# List remote branches with test data
nf-core test-datasets list-branches

# List files on a specific branch
nf-core test-datasets list --branch <branch-name>

# Search for specific datasets
nf-core test-datasets search <query> --branch <branch-name>
```

#### Dataset Discovery Examples

```bash
# List all datasets for a specific pipeline
nf-core test-datasets list --branch mag

# Search for specific test data
nf-core test-datasets search minigut_reads --branch mag

# Get download URLs for integration
nf-core test-datasets list --branch mag --generate-dl-url

# Get Nextflow-compatible paths
nf-core test-datasets list --branch mag --generate-nf-path
```

#### Command Options

The `list` and `search` commands support these options:

- `--branch`, `-b`: Specify the branch in test-datasets repository
- `--generate-nf-path`, `-p`: Auto-generate file paths for Nextflow code
- `--generate-dl-url`, `-u`: Auto-generate GitHub download URLs

## Next Steps

Continue to [CI/CD Integration](./09_cicd_integration.md) to learn about integrating tests with continuous integration.
