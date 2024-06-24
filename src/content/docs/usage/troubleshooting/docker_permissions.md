---
title: Docker permission errors
subtitle: How to troubleshoot common mistakes and issues
weight: 80
---

## One module container fails due to Docker permissions

The nf-core template `nextflow.config` contains the configuration `docker.runOptions = '-u $(id -u):$(id -g)'{:groovy}` for the profiles `docker` and `arm`.
This is done to emulate the user inside the container.

In some containers, this option may cause permission errors, for example when the Docker container writes to the `$HOME` directory, as the emulated user won't have a `$HOME` directory.

One solution is to override this Docker option with a config file, with `docker.runOptions = ''{:groovy}`. However, this change will affect all the Docker containers, and can't be overridden only in one single process.
To solve this, one option is to replace `docker.runOptions{:groovy}` and use `containerOptions` instead.

Your new `nextflow.config` file should look like this:

```groovy title="nextflow.config"
profiles {
  docker {
    process.containerOptions = '-u $(id -u):$(id -g)'
  }
  arm {
    process.containerOptions = '-u $(id -u):$(id -g)'
  }
}
```

And you can override this value for a particular process selecting it by name, in the `modules.config` file:

```groovy title="modules.config"
process {
  withName: <TOOL> {
        containerOptions = ''
    }
}
```

:::warning
As mentioned in the [Nextflow documentation](https://www.nextflow.io/docs/latest/process.html#containeroptions), the `containerOptions` feature is not supported by the Kubernetes and Google Life Sciences executors.
:::
