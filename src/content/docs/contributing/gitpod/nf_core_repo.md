---
title: Developing pipelines
subtitle: Creating a new module for the nf-core modules repository.
weight: 5
---

## Trying your first Gitpod environment

You can run Gitpod with any nf-core pipeline repository.

For example, for nf-core RNA-seq pipeline, simply click the green Gitpod button or add the Gitpod prefix before the git URL (instead of `https://github.com/nf-core/rnaseq`, type `https://gitpod.io/#https://github.com/nf-core/rnaseq`).

Once Gitpod has loaded, including the container with all the tools we need, we can go to the terminal and type the following to start the nf-core workflow:

```bash
nextflow run nf-core/rnaseq \
    -profile test,docker \
    --outdir my_result
```

This should run the test data through nf-core rnaseq, using docker with your results in the folder: `my_result`. This will take some time to complete.

Using this Gitpod method makes it easy to run and test nf-core pipelines quickly, but it lacks the parallelization required to run real datasets.

## Testing your module with pytest

With DSL2, much of pipeline writing is spent writing reusable [modules](https://nf-co.re/modules).
If you enter the [Gitpod environment for modules](https://gitpod.io/#https://github.com/nf-core/modules), you can run the pytest function in order to debug a particular module.

You can learn more about `pytest` in [this nf-core bytesize talk](https://nf-co.re/events/2021/bytesize-17-pytest-workflow).

Once you are in the environment, try running an example pytest for an existing module:

```bash
PROFILE=docker pytest --tag <module_name> --symlink --keep-workflow-wd --git-aware
```
