---
title: Your first pipeline
subtitle: Rune your first Nextflow and nf-core pipelines.
shortTitle: Your first pipeline
weight: 3
---

## Your first pipeline

Nextflow works best with an active internet connection, as it is able to fetch all pipeline requirements.

:::note
Running Nextflow without an active internet connection requires extra steps to move tooling, pipelines, and other assets offline. See [Running offline](/docs/users/configuration/runningoffline) for more information.
:::

### Hello world!

Run Nextflow’s example [Hello world!](https://github.com/nextflow-io/hello) pipeline to confirm Nextflow is installed and that you have an active internet connection:

```bash
nextflow run hello
```

If everything has been installed properly you should see the pipeline launch and run to completion:

```console
 N E X T F L O W   ~  version 24.04.4
Launching `https://github.com/nextflow-io/hello` [cheesy_sax] DSL2 - revision: afff16a9b4 [master]
executor >  local (4)
[2c/3092e1] sayHello (4) [100%] 4 of 4 ✔
Hello world!
Bonjour world!
Hola world!
Ciao world!
```

### nf-core pipelines

Each nf-core pipeline will have different run requirements. It is best to test it is working using a small dataset that is known to work before running a pipeline with your own data.

nf-core pipelines come packed with a `test` profile that uses small files from the [nf-core/test-datasets](https://github.com/nf-core/test-datasets) repository. The `test` profile is extensively tested and should run without extra configuration. nf-core pipelines come packed with directives for containers and environments that can be flexibly enabled using software dependency profiles (e.g., `docker`, `singularity`, and `conda`). By using a profile for software dependencies an nf-core pipeline should run with all required tools.

:::note
Profiles are not just for software management, but also deeper configuration on different computing infrastructure. Some institutes have profiles to help you run nf-core pipelines on your local infrastructure. See https://nf-co.re/configs to see if yours is listed!
:::

Run the `nf-core/rnaseq` pipeline using the `test` profile and a software dependencies profile (e.g., `docker`) to confirm Nextflow is installed and that you have an active internet connection:

```bash
nextflow run nf-core/rnaseq -profile test,docker --outdir results
```

If everything has been installed properly you should see the pipeline launch and run to completion.

What's next? Check out the rest of the [docs](/docs/) for more information about using, developing and contributing to nf-core pipelines.
