---
title: Configure a pipeline run
subtitle: A hands-on tutorial for configuring nf-core/demo with environment variables, parameters, and config files
shortTitle: Configure a pipeline run
---

nf-core pipelines are shaped by multiple levels of configuration: environment variables that tune Nextflow itself, parameters that tune the pipeline, and config files that tune execution.
Using [`nf-core/demo`](https://nf-co.re/demo), this tutorial explains each layer as worked examples.

![nf-core/demo subway map](../../../../assets/images/get_started/nf-core-demo-subway.png)

By the end you will have:

- Supplied pipeline parameters on the command line and from a `params.yaml` file.
- Activated profiles to switch between container engines and test scenarios.
- Layered a custom `custom.config` file with process-specific overrides.
- Used `ext.args` (the convention for `conf/modules.config`) to extend a tool's command line without editing the module.
- Used `NXF_*` environment variables to control Nextflow's runtime.
- Understood which configuration source overrides which.

:::note{title="Prerequisites"}
You will need the following to get started:

- An internet connection
- [Nextflow version 25.10.4 or later](../../get_started/environment_setup/nextflow.md)
- A container engine, such as Docker, Singularity, or Conda. See [Software dependencies](../../get_started/environment_setup/software-dependencies.md)

:::

## Baseline run

Before changing anything, run the pipeline with no custom configuration so you have something to compare against.
This baseline run relies entirely on `nf-core/demo`'s built-in defaults.
No environment variables, no custom parameters, and no extra config files are set.
In each of the following steps you will replace one of these defaults and observe the effect.

Start with a plain test invocation as your reference point, replacing `docker` with your preferred software dependency manager (for example, `singularity` or `conda`).

```bash
nextflow run nf-core/demo -r 1.2.0 -profile test,docker --outdir results
```

In the background, this run uses three defaults that you will override in the next steps:

- The pipeline's internal [`nextflow.config`](https://github.com/nf-core/demo/blob/master/nextflow.config) and the [`test`](https://github.com/nf-core/demo/blob/master/conf/test.config) and `docker` profiles
- Default pipeline parameters (path to input samplesheet, the customisable MultiQC title and so on). See the [nf-core/demo parameters page](https://nf-co.re/demo/parameters) and [`nextflow_schema.json`](https://github.com/nf-core/demo/blob/master/nextflow_schema.json)
- Your default Nextflow runtime settings (e.g., work directory and Nextflow version)

:::tip
Keep the `results/` directory from this run.
You can compare it to the output of later steps to confirm your configuration changes took effect.
:::

## Configure with parameters

Pipeline parameters are the knobs and switches the pipeline offers to users to control how the different steps of the pipeline execute.
For example, inputs, outputs, skip flags, reference data choices, and tool toggles.
They're documented per pipeline.
For `nf-core/demo`, see the [parameters reference](https://nf-co.re/demo/parameters) and [`nextflow_schema.json`](https://github.com/nf-core/demo/blob/master/nextflow_schema.json).

You can supply them two ways:

- Directly on the command line
- From a `params.yaml` or `params.json` file.

The two forms are interchangeable, and you can mix them.

:::important
CLI flags override values in a parameter file.
:::

1. Override a single parameter on the command line:

   ```bash
   nextflow run nf-core/demo -r 1.2.0 -profile test,docker --outdir results_customcommandline
   ```

   Output now goes to `results_customcommandline/` instead of `results/`.

2. For longer parameter sets, create a `params.yaml` file:

   ```yaml title="params.yaml"
   outdir: results_customparamsfile
   multiqc_title: "nf-core/demo parameter file configured run"
   ```

3. Apply it with `-params-file`:

   ```bash
   nextflow run nf-core/demo -r 1.2.0 -profile test,docker -params-file params.yaml
   ```

   Output now goes to `results_customparamsfile/` instead of `results/`.
   Open `results_customparamsfile/multiqc/nf-coredemo-parameter-file-configured-run_multiqc_report.html` in a web browser to confirm the report title was set.

:::tip
Both JSON and YAML formats are supported for parameter files.
Parameter files make runs easier to reproduce and share.
You can publish them alongside results so others can re-run the exact same configuration.
See [Pipeline parameters](https://docs.seqera.io/nextflow/cli#pipeline-parameters) for more information.
:::

:::warning
nf-core pipeline parameters must be supplied on the command line or via `-params-file`.
Use config files for executor and computational resource settings, and parameter files for pipeline parameters.
:::

## Configure with config files

Config files are the most flexible configuration layer.
They control both **where** the pipeline runs and **how** it runs.
For example, which executor submits jobs (e.g., `local`, `slurm`, `awsbatch`), how much CPU and memory each process should request, what Nextflow does when a process fails, which container registry to pull from, where intermediate work lives, how to apply different settings to specific steps, and much more.
They're the right place for any setting that depends on your infrastructure rather than on what the pipeline does scientifically.

Nextflow assembles its final configuration from several config files.
Nextflow picks up some config files automatically, some you switch on by name with `-profile`, and others you point it at explicitly with `-c`.

The sections that follow cover all three approaches and then show how to target specific processes inside any config file.

### Auto-loaded `nextflow.config` files

Nextflow looks for config files in three places without being told.

1. **The pipeline's own `nextflow.config` in the project directory.** Included with every nf-core pipeline's source code. For example, [`nf-core/demo`'s `nextflow.config`](https://github.com/nf-core/demo/blob/master/nextflow.config), which loads [`conf/base.config`](https://github.com/nf-core/demo/blob/master/conf/base.config) and declares every shipped profile. **Don't edit this file**. See the warning in [Configuration options](../configuration/configuration-options.md).
2. **`$HOME/.nextflow/config`.** Your personal config, applied to every Nextflow run you launch. Good for settings that should follow you across all pipelines, such as your Singularity cache location or a default executor for your workstation.
3. **`nextflow.config` in the launch directory.** Applied automatically whenever you run Nextflow from that directory. Good for per-project or per-working-directory settings you don't want to retype each time.

Add a `nextflow.config` in your current directory:

```groovy title="nextflow.config"
process {
  cpus = 1
  memory = 2.GB
}
```

Re-run the baseline command. No `-c` flag is needed, Nextflow picks the file up automatically:

```bash
nextflow run nf-core/demo -r 1.2.0 -profile test,docker --outdir results_customnextflowconfig
```

Compare between the 'Tasks' section `results/pipeline_info/results/execution_report_<datetimestamp>.html` and `results_customnextflowconfig/pipeline_info/results/execution_report_<datetimestamp>.html`.
Observe that the process `NFCORE_DEMO:DEMO:COWPY` has changed it's memory from `4.GB` to `2.GB`, as specified in the new `nextflow.config`.

:::tip
Use `$HOME/.nextflow/config` for personal defaults that follow you across machines, and a launch-directory `nextflow.config` for project-specific settings.
Reserve `-c` (covered below) for one-off overrides and shared project configs that should live alongside your pipeline command.
:::

### Activate bundled settings with profiles

A **profile** is a named bundle of configuration that lives inside a config file under the `profiles` scope.
They allow you to switch between different sets of configurations stored in a single config file.
You activate one or more with the `-profile` flag.
You've been using profiles since [Baseline run](#baseline-run). `-profile test,docker` activates two: `test` (a small public dataset) and `docker` (use Docker to manage software).

Every nf-core pipeline comes with a standard set of profiles:

- **Software profiles**: `docker`, `singularity`, `apptainer`, `podman`, `conda`, `charliecloud`, `shifter` (one per container engine or environment manager).
- **Test profiles**: [`test`](https://github.com/nf-core/demo/blob/master/conf/test.config) (a small dataset for quick verification) and [`test_full`](https://github.com/nf-core/demo/blob/master/conf/test_full.config) (the full-size dataset used in CI).

  :::note
  The [`test_full`](https://github.com/nf-core/demo/blob/master/conf/test_full.config) profile in nf-core/demo is also a very small test dataset that can be used for testing.
  In other pipelines these can be much larger, but produce realistic output.
  :::

- **Institutional profiles**: contributed to [nf-core/configs](https://github.com/nf-core/configs) and loaded automatically by every nf-core pipeline. Activate one with `-profile <institution>` if your cluster has one - see [Use shared institutional configs](#use-shared-institutional-configs).

Combine profiles with commas, and remember that order matters - later profiles override earlier ones where they overlap:

```bash
nextflow run nf-core/demo -r 1.2.0 -profile test,docker --outdir results
```

In this example, options in the `docker` profile will override any of the options with the same name as in `test`.

You can also define your own profile for settings you'd like to reuse.
Add a `mymachine` profile to your existing launch-directory `nextflow.config`:

```groovy title="nextflow.config"
profiles {
  mymachine {
    process {
      cpus = 2
    }
  }
}
```

Activate it alongside the existing profiles:

```bash
nextflow run nf-core/demo -r 1.2.0 -profile test,docker,mymachine --outdir results_customprofile
```

Compare between the 'Tasks' section `results_customnextflowconfig/pipeline_info/results/execution_report_<datetimestamp>.html` and `results_customprofile/pipeline_info/results/execution_report_<datetimestamp>.html`.
Observe that the process `NFCORE_DEMO:DEMO:COWPY` has changed the CPUs from `1` to `2`, as specified in the new `nextflow.config`.

:::note
The `mymachine` profile is illustrative — it requests more CPUs and memory than the `test` profile's small dataset needs.
On a workstation without those resources available, Nextflow will fail when it tries to run the processes.
:::

### Use shared institutional configs

If you work on a shared HPC cluster or cloud platform, there's a good chance someone has already written a profile for it.
The [nf-core/configs](https://github.com/nf-core/configs) repository collects over 150 cluster-specific configs contributed by the community, covering the executor, queue, resource limits, container engine, scratch paths, module systems, and any other cluster-specific quirks needed to run nf-core pipelines well in that environment.

Every nf-core pipeline loads these configs automatically — there's nothing to download or copy yourself.
At run time, each pipeline fetches the [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) from the `nf-core/configs` repository and makes every profile in it available alongside the pipeline's own profiles.
A handful of examples from the [configs directory](https://github.com/nf-core/configs/tree/master/conf):

- [`uppmax`](https://github.com/nf-core/configs/blob/master/conf/uppmax.config) — UPPMAX (Sweden)
- [`bih`](https://github.com/nf-core/configs/blob/master/conf/bih.config) — Berlin Institute of Health
- [`aws_tower`](https://github.com/nf-core/configs/blob/master/conf/aws_tower.config) — Seqera Platform on AWS
- [`crick`](https://github.com/nf-core/configs/blob/master/conf/crick.config) — Francis Crick Institute

Browse the full list at [nf-co.re/configs](https://nf-co.re/configs).

You activate one the same way you activate any other profile. If your institution's infrastructure is supported by nf-core/configs, replace `<institution>` with the profile name for your environment to try it out.

Some institutions also contribute **pipeline-specific** overrides — a profile that further tunes resources for a particular pipeline on a particular cluster.
These live under [`conf/pipeline/<pipeline>/<institution>.config`](https://github.com/nf-core/configs/tree/master/conf/pipeline) in the repository and are picked up automatically when you run that pipeline with the matching `-profile`.

:::tip
If your cluster doesn't have a profile yet, the [nf-core/configs contribution guide](https://github.com/nf-core/configs#documentation) describes how to add one.
Contributing a profile to `nf-core/configs` is the recommended way to share cluster configuration with your team — it removes the need for everyone to carry a `custom.config` file around, and benefits every nf-core pipeline at once.
:::

:::note
The shared configs are loaded over the network at run time.
If you're working on an air-gapped system, see [Running pipelines offline](../run-pipelines-offline.md) for how to download them ahead of time and point Nextflow at the local copy via the `custom_config_base` parameter.
:::

### Pass a config explicitly with `-c`

For configuration you don't want loaded by default.
For example, a one-off resource bump, a shared institutional config, or an experimental override, pass it on the command line with `-c`.

1. Create a small file named `custom.config`:

   ```groovy title="custom.config"
   process {
     memory = 3.GB
   }
   ```

2. Apply it with `-c`:

   ```bash
   nextflow run nf-core/demo -r 1.2.0 -profile test,docker -c custom.config --outdir results_customconfig
   ```

Compare between the 'Tasks' section `results_customprofile/pipeline_info/results/execution_report_<datetimestamp>.html` and `results_customcontig/pipeline_info/results/execution_report_<datetimestamp>.html`
Observe that the process `NFCORE_DEMO:DEMO:COWPY` has changed the memory from `2` to `3`, as specified in the new `custom.config`.
Furthermore, the CPUs have gone back from 2 to 1 as we did not specify the `mymachine` profile, thus falling back to the value in `nextflow.config`.

You can pass `-c` more than once.
Files are applied in order, so later ones override earlier ones.

### Tune how the pipeline executes

Beyond resource requests, config files control **how** each run executes.
For example, which scheduler picks up jobs, what to do when a process fails, and where intermediate files live.
These settings sit alongside the `process` scope you used above, with a few sibling top-level options.

Edit `custom.config` to add automatic retries and redirect the work directory:

```groovy title="custom.config" {1,6-7}
workDir = 'nf-work'

process {
  memory = 3.GB
  cpus = 1
  errorStrategy = 'retry'
  maxRetries    = 2
}
```

Apply it the same way as before:

```bash
nextflow run nf-core/demo -r 1.2.0 -profile test,docker -c custom.config --outdir results_customconfig2
```

You will now see that an `nf-work/` directory specified in the config was generated in the directory where you executed the command, alongside the results directory.

Each setting controls a different aspect of the run.
Common settings are:

- **`process.executor`**: where jobs go — `local` (default), `slurm`, `awsbatch`, `lsf`, `pbs`, `kubernetes`, and so on.
- **`process.queue`**: which queue or partition to submit to on shared clusters.
- **`process.errorStrategy`** with **`process.maxRetries`**: retry transient failures (for example, a node going down or a timeout) instead of failing the whole run.
- **`workDir`**: where Nextflow stages intermediate files — point it at fast scratch storage on HPC to keep your home filesystem clean.

:::info
The syntax

For the full list of config scopes and options (`process`, `executor`, `docker`, `singularity`, `aws`, `azure`, `google`, `report`, `trace`, `timeline`, and more), see the [Nextflow config reference](https://docs.seqera.io/nextflow/reference/config).
nf-core pipelines already enable the standard execution reports under `pipeline_info/`, so you only need to override those if you want custom paths.
:::

### Target specific processes with `withName` and `withLabel`

The blanket `process { cpus = 3 }` example above applies to every process in the pipeline.
In practice you'll usually want to target specific steps of the pipeline.
For example, to bump the resources for an aligner without changing anything else, or pass an extra command-line argument to one tool.
Nextflow gives you two selectors for this:

- **`withName`** matches a specific process by its fully qualified name. `nf-core/demo` runs three processes — [`FASTQC`](https://github.com/nf-core/demo/blob/master/modules/nf-core/fastqc/main.nf), [`SEQTK_TRIM`](https://github.com/nf-core/demo/blob/master/modules/nf-core/seqtk/trim/main.nf), and [`MULTIQC`](https://github.com/nf-core/demo/blob/master/modules/nf-core/multiqc/main.nf) — and each has a name you can see in `.nextflow.log` or in `pipeline_info/execution_trace.txt` after a run.
- **`withLabel`** matches every process that carries a given label. All nf-core pipelines tag their processes with size labels (`process_low`, `process_medium`, `process_high`, `process_high_memory`) — see how they map to resource requests in [`nf-core/demo`'s `conf/base.config`](https://github.com/nf-core/demo/blob/master/conf/base.config). A single `withLabel` block can adjust whole resource tiers at once.

Edit your `custom.config` to include both in the `process block:

```groovy title="custom.config" {8-16}
workDir = 'nf-work'

process {
  memory = 3.GB
  cpus = 1
  errorStrategy = 'retry'
  maxRetries    = 2

  withName: 'NFCORE_DEMO:DEMO:FASTQC' {
    cpus = 1
    memory = 3.GB
  }

  withLabel: 'process_low' {
    cpus = 1
    memory = 1.GB
  }
}
```

Re-run with the `-c custom.config`.

```bash
nextflow run nf-core/demo -r 1.2.0 -profile test,docker -c custom.config --outdir results_customconfig3
```

Compare `results_customconfig2/pipeline_info/execution_report_<datetimestamp>.html` with `results_customconfig3/pipeline_info/execution_report_<datetimestamp>.html`.
FASTQC should report requesting 1 CPU and 3 GB memory instead of 2 and 4 GB.

`SEQTK_TRIM` has also changed, with now at 1 GB of memory and 1 CPU instead of 2 and 4 GB.
Process labels are declared on each `process` definition inside the module.
For example, [`SEQTK_TRIM`'s `main.nf`](https://github.com/nf-core/demo/blob/master/modules/nf-core/seqtk/trim/main.nf) declares `label 'process_low'` on its second line.

:::note
The size labels (`process_low`, `process_medium`, `process_high`, `process_high_memory`) are standard across all nf-core pipelines.
:::

### Pass tool arguments with `ext.args`

`ext.args` is a per-process configuration value that nf-core modules use to inject extra command-line flags into the underlying tool.
It lives inside a `withName` selector (so it follows the same targeting rules as the previous section), and the module picks it up at runtime with `def args = task.ext.args ?: ''`.

By convention, nf-core pipelines collect every per-process option in a dedicated file, [`conf/modules.config`](https://github.com/nf-core/demo/blob/master/conf/modules.config).
It is another source of configuration loaded automatically through the pipeline's `nextflow.config`.
Open `nf-core/demo`'s [`conf/modules.config`](https://github.com/nf-core/demo/blob/45904cb9d12db3d89900e6c479fe604ef71b297b/conf/modules.config#L21-L28) to see two real uses of `ext.args`:

```groovy {2,7}
withName: FASTQC {
    ext.args = '--quiet'
    // ...
}

withName: 'MULTIQC' {
    ext.args = { params.multiqc_title ? "--title \"$params.multiqc_title\"" : '' }
    // ...
}
```

FASTQC always runs with FastQC's `--quiet` mode.
MULTIQC dynamically picks up `--title` if the `multiqc_title` parameter is set.
This is why the `multiqc_title` value from your `params.yaml` flowed through into the MultiQC report.

To pass your own flags without forking the module, override `ext.args` in `custom.config`:

```groovy title="custom.config" {12,19-21}
workDir = 'nf-work'

process {
  memory = 3.GB
  cpus = 1
  errorStrategy = 'retry'
  maxRetries    = 2

  withName: 'NFCORE_DEMO:DEMO:FASTQC' {
    memory = 3.GB
    cpus = 1
    ext.args = '--quiet  --nogroup'
  }

  withLabel: 'process_low' {
    memory = 1.GB
    cpus = 1
  }

  withName: 'MULTIQC' {
      ext.args = { params.multiqc_title ? "--title \"$params.multiqc_title\"" : '' }
  }
}
```

Re-run with the `-c custom.config`.

```bash
nextflow run nf-core/demo -r 1.2.0 -profile test,docker -c custom.config --outdir results_customconfig4
```

FASTQC will pick up the new flags through `task.ext.args` inside its [`main.nf`](https://github.com/nf-core/demo/blob/master/modules/nf-core/fastqc/main.nf).
Compare in your web browser `results_customconfig3/fastqc/SAMPLE1_PE/SAMPLE1_PE_1_fastqc.html` and `results_customconfig4/fastqc/SAMPLE1_PE/SAMPLE1_PE_1_fastqc.html`.
You will see the 'Per base sequence quality' plot has changed, due to addition of the `--nogroup` option.

:::danger
Customising `ext.args` is generally not recommended, and may break the pipeline as the changes have not been tested by the developer.
It is recommended to request official support for a new tool option or argument from the pipeline developer.
Customise `ext.args` in a config as a last resort.
:::

Sibling keys you'll see in the same file:

- **`ext.args2`, `ext.args3`** — second and third argument sets for modules that call more than one tool.
- **`ext.prefix`** — overrides the prefix used for the module's output filenames.

:::tip
`ext.args` is the right way to pass tool-specific flags that aren't exposed as pipeline parameters.
Never edit the module's `main.nf` to add a flag — `ext.args` exists precisely so you can extend the command line without changing the versioned pipeline code.
:::

For `resourceLimits`, executor tuning, and container registry overrides, see [System requirements](../configuration/nextflow-for-your-system.md).
For profile precedence and shared institutional configs, see [Configuration options](../configuration/configuration-options.md).

## Configure with environment variables

<!-- TODO JAMES UP TO HERE -->

`NXF_*` environment variables are the outermost configuration layer.
Nextflow reads them from your shell as it starts, before any config file or parameter is parsed, so they're the right place for settings that don't change between runs.
For example, which Nextflow version to use, where the work directory and container cache live, and whether to operate offline.
They configure **Nextflow itself**, not the pipeline.
You cannot set pipeline inputs or outputs this way.

You set them like any other shell variable: with `export`, inline before a single command, or by adding them to your shell's startup file.

1. Pin a Nextflow version and redirect the work directory for this run:

   ```bash
   export NXF_VER=25.10.4
   export NXF_WORK=$HOME/nf-demo-work
   nextflow run nf-core/demo -r 1.2.0 -profile test,docker --outdir results_customenv
   ```

   The version line at the top of stdout (and in `.nextflow.log`) should now report `25.10.4`, and intermediate files should appear under `$HOME/nf-demo-work/` instead of `./work/`.

2. To set a variable for a single command without exporting it, prefix the invocation:

   ```bash
   NXF_VER=25.10.4 nextflow run nf-core/demo -r 1.2.0 -profile test,docker --outdir results_customenv2
   ```

3. To persist a variable across sessions, add the export to your shell config:

   ```bash title=".bashrc"
   export NXF_VER=25.10.4
   ```

:::tip
The `NXF_*` variables nf-core users reach for most often are:

- `NXF_VER`: Pin a specific Nextflow version
- `NXF_HOME`: Override Nextflow's home directory (defaults to `$HOME/.nextflow`)
- `NXF_WORK`: Override the work directory for intermediate files
- `NXF_OFFLINE`: Run without contacting the network (see [Running pipelines offline](../run-pipelines-offline.md))
- `NXF_SINGULARITY_CACHEDIR` / `NXF_APPTAINER_CACHEDIR`: Reuse downloaded container images across runs

For the full list, see the [Nextflow environment variables reference](https://docs.seqera.io/nextflow/reference/env-vars).
:::

:::note
`NXF_*` variables are distinct from pipeline parameters.
They configure how Nextflow runs, not what the pipeline does.
You cannot set `--input` or `--outdir` with an environment variable.
:::

## Put it all together

Each layer you've explored solves a different problem.
In practice, a single run usually uses several at once.
This step combines them into one invocation so you can see how they interact.

1. Set environment variables for runtime concerns that don't change between runs:

   ```bash
   export NXF_VER=25.10.4
   export NXF_WORK=$HOME/nf-demo-work
   ```

2. Capture pipeline parameters in a `params.yaml` file so the run is reproducible:

   ```yaml title="params.yaml"
   outdir: my_results
   multiqc_title: "nf-core/demo configured run"
   ```

3. Put per-process overrides in `custom.config`:

   ```groovy title="custom.config"
   process {
     withName: 'NFCORE_DEMO:DEMO:FASTQC' {
       cpus = 4
       memory = 8.GB
       ext.args = '--quiet --nogroup'
     }

     withLabel: 'process_low' {
       cpus = 2
       memory = 4.GB
     }
   }
   ```

4. Launch the run, activating the `test` and `docker` profiles and passing the parameter and config files:

   ```bash
   nextflow run nf-core/demo -r 1.2.0 \
     -profile test,docker \
     -params-file params.yaml \
     -c custom.config
   ```

This single command exercises every layer the tutorial introduced.
Each one is resolved independently:

- **Runtime**: `NXF_VER` and `NXF_WORK` are read from the shell before Nextflow parses any config, pinning the version and redirecting intermediate files.
- **Configuration**: config files layer in a fixed order — the pipeline's bundled `nextflow.config` (which pulls in `conf/base.config` and `conf/modules.config`) loads first, the `test` and `docker` profiles override matching keys, and `custom.config` overrides them last. The `withName` and `withLabel` blocks in `custom.config` therefore win for FASTQC and every `process_low` step.
- **Parameters**: values in `params.yaml` set the pipeline parameters. A matching `--flag` on the command line would override them.

To confirm each layer took effect, check:

- `my_results/pipeline_info/execution_report.html` — FASTQC should report 4 CPUs and 8 GB of memory.
- `my_results/multiqc/multiqc_report.html` — the report title should read "nf-core/demo configured run".
- `$HOME/nf-demo-work/` — intermediate files should appear here instead of under `./work/`.

:::tip
A rule of thumb for where each setting belongs:

- **Environment variables** for runtime concerns tied to your machine or session (Nextflow version, cache locations, offline mode).
- **`params.yaml`** for any setting the pipeline exposes as a parameter.
- **Profiles** for reusable bundles you switch on by name (container engine, test data, institutional cluster).
- **`custom.config`** for one-off or team-specific overrides that don't warrant a profile.
  :::

## Next steps

- Learn more about [configuration options](../configuration/configuration-options.md) (i.e., profiles, shared nf-core/configs, and full precedence rules)
- Learn more about [system requirements](../configuration/nextflow-for-your-system.md) (i.e., resource limits, executors, and tool argument overrides)
- Learn more about [Nextflow configuration](https://docs.seqera.io/nextflow/config)
