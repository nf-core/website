---
title: 5 - Developing pipelines
subtitle: Creating a new module for the nf-core modules repository.
---

## Trying your first Gitpod environment

You can run Gitpod with the any nf-core pipeline repository.

For example, for nf-core RNA-seq pipeline, simply click the green gitpod button or add the Gitpod prefix before the git URL (of `https://gitpod.io/#` ; to become: `https://gitpod.io/#https://github.com/nf-core/rnaseq`).

Once Gitpod has loaded, including the container with all the tools we need, we can go to the terminal and type the following to start the nf-core workflow:

```bash
nextflow run nf-core/rnaseq \
    -profile test,docker \
    --outdir my_result
```

This should run the test data through nf-core rnaseq, using docker with your results in the folder: `my_result`. This will take some time to complete.

Using this Gitpod method makes it easy to run and test nf-core pipelines quickly, but lacks the parallelization required to run real datasets.

## Testing your module with pytest

With DSL2, much of pipeline writing is spent writing reusable [modules](https://nf-co.re/modules).
If you enter the [gitpod environment for modules](https://gitpod.io/#https://github.com/nf-core/modules), you can run the pytest function in order to debug a particular module.

You can learn more about `pytest` in [this nf-core bytesize talk](https://nf-co.re/events/2021/bytesize-17-pytest-workflow).

Once you are in the environment, try runnning an example pytest for an exisiting module:

```bash
PROFILE=docker pytest --tag <module_name> --symlink --keep-workflow-wd --git-aware
```
