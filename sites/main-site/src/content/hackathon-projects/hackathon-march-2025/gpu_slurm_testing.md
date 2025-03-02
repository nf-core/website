---
title: GPU development / testing
category: components
slack: "https://nfcore.slack.com/archives/C07PPQR415W"
intro_video: ""
image: "/assets/images/events/2025/hackathon-march/gpu_slurm_tests_nvidia_h100.jpg"
image_alt: Picture of NVIDIA H100 Graphics card
leaders:
  Florian Wuennemann:
    name: Florian Wuennemann
    slack: "https://nfcore.slack.com/archives/DTZKT23D1"
---

Thanks to the generous sponsorship from [Nebius](https://nebius.com/), we have access to 16 powerful H100 GPUs for the duration of this hackathon. Nebius is providing these resources through their Slurm-based workload manager called [Soperator](https://nebius.com/services/soperator), enabling participants to test and develop GPU-accelerated tools with Nextflow using a familiar Slurm interface.
We welcome people from all other projects to use these resources if they need access to GPUs.

To enable participants access to the Slurm cluster with the GPU resources, we need to add your public key to the authorized_keys list on the Slurm cluster. If you want access, please contact [Florian](https://nfcore.slack.com/archives/DTZKT23D1) with your public key prior to the hackathon start.


You can create a dedicated ssh keypair with ssh-keygen like this (store in a safe location):

`ssh-keygen -t rsa -C "your_email@example.com"`

Once your public key has been added to the authorized_keys list on the cluster, you can connect to it with:

`ssh -i /path/to/your/private_key username@hostname`

:::warning
**Important information**: Access to the cluster and all data on the cluster will be wiped once the hackathon ends, so make sure to not put any critical data on the Slurm data storage!
:::

## Goal

Use GPUs via Slurm to test and develop GPU accelerated tools in Nextflow.
