---
title: Developing pipelines
shortTitle: Developing pipelines
subtitle: Developing pipelines and modules in Codespaces
weight: 3
type: "tutorial"
---

## Running a pipeline in Codespaces

You can run Codespaces from any nf-core pipeline repository on Github.
It also works for the nf-core tools, modules, and website repositories.

Take the nf-core/rnaseq pipeline for example:
It can be test-run from its [repository on Github](https://github.com/nf-core/rnaseq) by clicking the green button labelled "Code" and then selecting "Create codespace" from the Codespaces tab.
Alternatively, visit [github.com/codespaces/new/nf-core/rnaseq](https://github.com/codespaces/new/nf-core/rnaseq) to have more control over the configuration.
Make sure to also select the larger 4 cpu hardware from "Machine Type".

Once Codespaces has loaded, a container with all required tools is run automatically and the following command can run the nf-core workflow:

```bash
nextflow run . \
    -profile test,singularity \
    --outdir my_result
```

:::note
We recommend the use of the singularity profile with `-profile singularity{:bash}` to successfully run the pipeline.
Although docker is available and containers can be run inside the nf-core devcontainers, running nextflow with the docker profile using `-profile docker{:bash}` is currently not yet supported.
:::

This should run the test data through nf-core/rnaseq, using apptainer with your results in the folder: `my_result`.
This may take some time to complete.

The devcontainer method in Codespaces makes it easy to run and test nf-core pipelines quickly, but it lacks the computational power required to run real-size datasets.

## Testing your module with nf-test

If you enter the [Codespaces environment for the modules branch](https://github.com/nf-core/modules), you can run nf-test in order to debug a particular module.

Once you are in the environment, try running an example nf-test for an existing module:

```bash
nf-test test --tag <module_name> --profile singularity
```
