---
title: Running GPU-accelerated pipelines
subtitle: Enable and configure GPU acceleration for supported pipeline steps
shortTitle: GPU pipelines
---

Some nf-core pipelines can run specific steps on a GPU instead of a CPU, for example,
tools like NVIDIA Parabricks, GPU-accelerated basecallers, or machine-learning-based classifiers.
GPU support is opt-in and pipeline-specific: there is no single flag that turns on GPU acceleration
across all of nf-core. Always check the pipeline's own `docs/usage.md` for the exact parameter names
before you run.

Currently, GPU acceleration typically assumes you have Nvidia/CUDA capable cards (e.g. AMD GPUs are unfortunately currently not supported).

This page assumes you're already comfortable with `-profile` and custom configuration files — see [Configuration options](./configuration-options) first if not.

## How pipelines expose GPU support

Two patterns currently exist across nf-core pipelines, and you need to know which one your pipeline uses.

### Pattern 1: a standalone parameter, no profile required

Some pipelines gate GPU acceleration entirely behind a boolean parameter, with no separate profile needed.
For example, [nf-core/rnaseq](https://nf-co.re/rnaseq) accelerates rRNA removal and STAR alignment this way:

```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --remove_ribo_rna --ribo_removal_tool ribodetector --use_gpu_ribodetector \
    -profile docker
```

Here, `--use_gpu_ribodetector` is enough on its own — `-profile docker` (or `singularity`) is the same profile you would use for a CPU-only run.
The GPU container flags (`--gpus all` for Docker, `--nv` for Singularity/Apptainer) are applied automatically, and only to the GPU-enabled step.
Other steps in the same run are unaffected and continue to run normally on CPU-only nodes, which matters on a shared or mixed cluster.


### Pattern 2: a feature parameter combined with `-profile gpu`

Other GPU integrations instead require you to add a `gpu` profile alongside a feature parameter that selects the GPU-based tool.
For example, [nf-core/sarek](https://nf-co.re/sarek):

```bash
nextflow run nf-core/sarek \
    --input samplesheet.csv \
    --outdir results \
    --aligner parabricks \
    -profile docker,gpu
```

or [nf-core/methylseq](https://nf-co.re/methylseq):

```bash
nextflow run nf-core/methylseq \
    --aligner bwameth \
    -profile gpu \
    --input samplesheet.csv --genome GRCh38
```

:::caution{title="`-profile gpu` applies to the whole run"}
In this pattern, the `gpu` profile might set the container engine's GPU flags globally (`docker.runOptions`, `singularity.runOptions`, etc.), so every process in the run gets them — not just the GPU-accelerated one - check the pipeline configuration.

Check the pipeline's usage docs for whether you also need to set the `accelerator` directive yourself in a custom config — some pipelines only set it automatically on cloud/Kubernetes executors (`awsbatch`, `google-batch`, `k8s`, `hq`), leaving it up to you elsewhere.
:::

## Running on a shared cluster (HPC)

If your institution has a profile in [nf-core/configs](https://github.com/nf-core/configs), GPU-enabled processes are usually already routed to the right queue or partition, with the correct scheduler flags (SLURM `--gres=gpu:N`, LSF `-gpu num=N`, PBS `-l select=...:ngpus=N`) added for you — you don't need to work these out yourself.

Combine your institution's profile with whichever GPU mechanism your pipeline needs:

```bash
# Pattern 1 pipeline (institution profile only)
nextflow run nf-core/rnaseq --use_gpu_ribodetector -profile <institution>

# Pattern 2 pipeline (institution profile plus the pipeline's own gpu profile)
nextflow run nf-core/sarek --aligner parabricks -profile gpu,<institution>
```

Check [nf-core/configs](https://github.com/nf-core/configs) to see if your cluster already has a profile, and read its `docs/<institution>.md` page — most include a copy-paste command for GPU runs on that system.

:::note{title="Why profile order usually isn't a problem here"}
`-profile a,b,c` order decides which profile's setting wins when two profiles assign the same config key — but that assignment happens once, when Nextflow loads your config. Most GPU-related settings (`clusterOptions`, `containerOptions`, `accelerator` inside a `withLabel:` block) are written as `{ }` closures, so what actually gets decided at load time is _which closure_ is in charge — that closure still runs fresh for every task, using that task's real values, not something baked in early. By convention, a pipeline's own `gpu` profile only sets `accelerator.request` (the GPU count); an institutional profile is responsible for `accelerator.type` (which GPU model) — so the two aren't fighting over the same setting.

The one place order still genuinely matters is a container engine's global `runOptions` (`docker.runOptions`, `singularity.runOptions`, etc.). Unlike `containerOptions` on a process, `runOptions` has no per-process scoping and can't be written as a closure on `task.accelerator` — it applies to every container the engine runs, full stop. That's also why the "applies to the whole run" caution above exists: older Pattern 2 pipelines like sarek set the GPU flag via `runOptions` rather than a process-scoped `containerOptions`. If your institutional profile also sets `runOptions` (for example, to add a bind mount), whichever profile is listed last fully replaces the other's value rather than combining with it — so you could silently lose either the GPU flag or the bind mount. This is specifically a `runOptions` problem: `containerOptions`, as used in the [institutional profile examples](../../developing/institutional-profiles/configuration#gpu-resource-requests), doesn't have it, since it's scoped per-process and only applies where you attach it. Check your institution's own `docs/<institution>.md` for the profile order it recommends for GPU runs.
:::

## Overriding GPU settings for a single run

To request more GPUs or target specific hardware for one run, add a custom config with [`-c`](./configuration-options#custom-configuration-files) that targets the process by name:

```groovy title="custom-gpu.config"
process {
    withName: 'RIBODETECTOR' {
        accelerator = [request: 2, type: 'a100']
    }
}
```



Only set `clusterOptions` yourself if there's no institutional profile doing this translation, or its `clusterOptions` hardcodes a flag instead of reading `task.accelerator`:

```groovy title="custom-gpu.config"
process {
    withName: 'RIBODETECTOR' {
        accelerator = [request: 2, type: 'a100']
        clusterOptions = '--gpus=a100:2'
    }
}
```

See [Customize process resources](./nextflow-for-your-system#customize-process-resources) for how to find the exact process name to target, including cases where a tool is used more than once in the same pipeline.

## See also

- [Configuration options](./configuration-options)
- [Institutional profiles](../../developing/institutional-profiles/overview)
- [GPU-capable modules](../../developing/components/gpu-modules) — for pipeline and module developers
