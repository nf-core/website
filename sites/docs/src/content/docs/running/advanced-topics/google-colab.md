---
title: Run nf-core pipelines on Google Colab
subtitle: A guide for using cloud resources with Google Colab
shortTitle: Google Colab
parentWeight: 100
---

This guide enables you to run computationally intensive nf-core pipelines using Google Colab's cloud resources, addressing limitations in local computing environments. Google Colab provides free access to cloud computing resources, making it a useful option for researchers who lack access to high-performance computing infrastructure. You can run pipelines directly in the browser or integrate with VS Code for enhanced development capabilities.

:::warning
Google Colab sessions have limitations including timeouts (typically 12 hours), restricted storage (approximately 100 GB), and no root access. For production workloads, consider dedicated cloud computing resources.
:::

:::note{title="Prerequisites"}
You will need the following to get started:

- A Google account to access [Google Colab](https://colab.research.google.com/)
- Basic familiarity with [Nextflow](/get_started/environment_setup/nextflow/) and nf-core pipelines
- Understanding of [conda environments](/get_started/environment_setup/software-dependencies/) (Docker and Singularity require root access unavailable in Colab)
  :::

## Install Java

Nextflow requires Java 17 or later. Install Java in your Colab environment:

```bash
!apt update
!apt install openjdk-17-jdk
!export JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
!source ~/.bashrc
```

## Install Nextflow

Download and install Nextflow:

```bash
!wget -qO- https://get.nextflow.io | bash
!mv nextflow /usr/bin/nextflow
!chmod +x /usr/bin/nextflow
!nextflow -v
```

The final command verifies the installation by displaying the Nextflow version.

## Configure conda

Docker and Singularity require root access, which is unavailable in Google Colab. Use conda as the software dependency manager instead.

Install condacolab:

```bash
!pip install -q condacolab
import condacolab
condacolab.install()
```

Configure conda channels for bioinformatics packages:

```bash
!conda config --add channels bioconda
!conda config --add channels conda-forge
!conda config --set channel_priority strict
```

## Run a pipeline

Test your setup by running the nf-core demo pipeline:

```bash
!nextflow pull nf-core/demo
!nextflow run nf-core/demo -profile conda,test --outdir demo-results
```

This downloads the demo pipeline and executes it with the test profile, which uses small datasets to verify the installation.

For other pipelines, replace `nf-core/demo` with the pipeline name and adjust parameters according to the pipeline documentation:

```bash
!nextflow run nf-core/<pipeline_name> -profile conda,test --outdir <output_directory>
```

## Integrate with VS Code

You can connect your Colab environment to VS Code for enhanced development capabilities using the vscode-colab extension.

Install and configure vscode-colab:

```python
!pip install vscode-colab
import vscode_colab
vscode_colab.login()  # Authenticate with GitHub
vscode_colab.connect(name="<tunnel_name>", git_user_name="<github_username>")
```

Connect from VS Code desktop:

1. Install the **Remote Tunnels** extension in VS Code
2. Sign in with the same GitHub account used in the Colab authentication
3. Run the **Remote Tunnels: Connect to Tunnel...** command
4. Select your tunnel from the list

## Manage data and outputs

Google Colab provides temporary storage that persists only for the duration of your session. For data persistence, integrate external storage solutions.

### Mount Google Drive

Mount Google Drive to store smaller outputs:

```python
from google.colab import drive
drive.mount('/content/drive')
```

After mounting, specify output directories within your Drive mount point:

```bash
!nextflow run nf-core/<pipeline_name> -profile conda,test --outdir /content/drive/MyDrive/<output_directory>
```

:::warning
Google Drive has storage limitations and may not accommodate large genomic datasets. For larger outputs, use institutional storage solutions or cloud storage services.
:::

### External storage options

For larger datasets, consider:

- Cloud storage services (AWS S3, Google Cloud Storage, Azure Blob Storage)
- Institutional storage systems accessible via network protocols
- Regular synchronisation to local storage during pipeline execution

## Limitations and best practices

Google Colab imposes several constraints that affect pipeline execution. Understanding these limitations helps you work around them effectively.

### Platform constraints

- **Session timeouts**: Sessions typically terminate after 12 hours of inactivity or 24 hours of continuous use
- **No root access**: System-level modifications and certain container technologies are unavailable
- **Storage constraints**: Available storage is limited to approximately 100 GB
- **Computing resources**: Resources are shared and may be less performant than dedicated cloud virtual machines
- **Network restrictions**: Some network operations may be limited or slower than dedicated infrastructure

### Working within constraints

Address these limitations through the following practices:

- **Manage session interruptions**: Use the `-resume` flag to restart pipelines after session timeouts. Commit code changes to version control regularly to prevent loss when sessions terminate
- **Optimize storage usage**: Start with small test datasets to verify configurations. Monitor storage throughout execution and configure intermediate file cleanup to manage consumption
- **Protect your outputs**: Save critical results to Google Drive or external storage throughout execution rather than waiting until completion. Back up results frequently
- **Maintain reproducibility**: Document your Colab setup (package versions, configurations) in version control for consistent results across sessions

## Troubleshooting

### Session timeout during execution

If your session times out mid-execution, you can resume the pipeline:

```bash
!nextflow run nf-core/<pipeline_name> -profile conda,test --outdir <output_directory> -resume
```

Nextflow's resume functionality relies on the work directory remaining intact. If the work directory was cleared, you'll need to restart the pipeline.

### Storage space exhausted

Monitor available storage:

```bash
!df -h
```

If storage is insufficient, consider:

- Using a smaller test dataset
- Configuring automatic cleanup of intermediate files
- Streaming outputs to external storage as they're generated

### Conda environment issues

If conda package installation fails, verify channel configuration:

```bash
!conda config --show channels
```

The output should show `bioconda` and `conda-forge` with strict priority enabled.
