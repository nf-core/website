---
title: Modifying pipelines
subtitle: Configure tool containers and arguments
shortTitle: Modifying pipelines
weight: 3
---

## Modifying tools

Each tool in an nf-core pipeline come preconfigured with a set arguments for an average user.
The arguments are a great place to start and have been tested as a part of the development process.
You normally can change the default settings using parameters using the double dash notation, e.g., `--input`.
However, you may want to modify these to fit your own purposes.

It is **very unlikely** that you will need to edit the pipeline code to configure a tool.

### Tool arguments

You may wish to understand which tool arguments a pipeline uses, update, or add additional arguments not currently supported by a pipeline.

You can sometimes find out what parameters are used by a tool in by checking the longer 'help' description of different pipeline parameters, e.g., by pressing the 'help' button next to [this parameter](https://nf-co.re/funcscan/1.0.1/parameters#annotation_bakta_mincontig) in [nf-core/funcscan](https://nf-co.re/funcscan).

There are two main places that a tool can have a tool argument specified:

- The process `script` block
- The `conf/modules.conf` file

Most arguments (both mandatory or optional) are defined in the `conf/modules.conf` file under the `ext.args` entry. Arguments that are defined in the `conf/modules.conf` file can be flexible modified using custom configuration files.

Arguments specified in `ext.args` are then inserted into the module itself via the `$args` variable in the module's bash code

For example, the `-n` parameter could be added to the `BOWTIE_BUILD` process:

```groovy
process {
    withName: BOWTIE_BUILD {
        ext.args = "-n 0.1"
    }
```

Updated tools may come with major changes and may break a pipeline and/or create missing values in MultiQC version tables.

:::warning
Such changes come with no warranty or support by the the pipeline developers!
:::

### Changing tool versions

You can tell the pipeline to use a different container image within a config file and the `process` scope.

You then need to identify the `process` name and override the Nextflow `container` or `conda` definition using the `withName` process selector.

For example, the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline uses a tool called Pangolin that updates an internal database of COVID-19 lineages quite frequently.

To update the container specification, you can do the following steps:

1. Check the default version used by the pipeline in the module file for the tool under `modules/nf-core/` directory of the pipeline. For example, for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags) for Docker or [Galaxy Project](https://depot.galaxyproject.org/singularity/) for Singularity
   - Note the container version tag is identical for both container systems, but must include the 'build' ID (e.g.`--pyhdfd78af_1`)
3. Create the custom config accordingly:

   - For Docker:

     ```groovy
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.1.17--pyhdfd78af_1'
         }
     }
     ```

   - For Singularity:

     ```groovy
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.1.17--pyhdfd78af_1'
         }
     }
     ```

   - For Conda:

     ```groovy
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.1.17'
         }
     }
     ```

:::warning
Updated tools may come with major changes and may break a pipeline and/or create missing values in MultiQC version tables. Such changes come with no warranty or support by the the pipeline developers.
:::

### Docker registries

nf-core pipelines use `quay.io` as the default docker registry for Docker and Podman images.
When specifying a Docker container, it will pull the image from `quay.io` unless a full URI is specified. For example, if the process container is:

```bash
biocontainers/fastqc:0.11.7--4
```

The image will be pulled from quay.io by default, resulting in a full URI of:

```bash
quay.io/biocontainers/fastqc:0.11.7--4
```

If `docker.registry` is specified, it will be used first. For example, if the config value `docker.registry = 'public.ecr.aws'` is specified the image will be pulled from:

```bash
public.ecr.aws/biocontainers/fastqc:0.11.7--4
```

However, the `docker.registry` setting will be ignored if you specify a full URI:

```bash
docker.io/biocontainers/fastqc:v0.11.9_cv8
```

:::warning
Updated registries may come with unexpected changes and come with no warranty or support by the the pipeline developers.
:::
